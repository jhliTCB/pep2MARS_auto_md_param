#!/public/home/lijunhao/soft/conda_envs/fldev/bin/python
# -*- coding: utf-8 -*-

"""
Simple optimization + atom-name mapping back to PDB template
  - Add Hs (keeps coordinates with addCoords=True)
  - Force all amide/peptide omegas to trans (~ -180°) with MMFF torsion constraints
  - Constrained minimization with RDKit MMFF (or UFF if MMFF fails)

Mapping back to PDB (index-based, no tolerance):
  - Replace ONLY the XYZ columns (30:54) in template PDB ATOM/HETATM lines
  - Coordinates come from the final optimized SDF
  - Atom ordering MUST match between template PDB and SDF.
"""

import os
import copy
import shutil
import tempfile
import subprocess
import argparse

from rdkit import Chem
from rdkit.Chem import AllChem, rdMolTransforms

def _is_carbonyl_carbon(atom):
    """Carbonyl carbon: C with a double-bonded O neighbor."""
    if atom.GetAtomicNum() != 6:
        return False
    mol = atom.GetOwningMol()
    aidx = atom.GetIdx()
    for nb in atom.GetNeighbors():
        if nb.GetAtomicNum() == 8:
            b = mol.GetBondBetweenAtoms(aidx, nb.GetIdx())
            if b is not None and b.GetBondType() == Chem.BondType.DOUBLE:
                return True
    return False

def _pick_acyl_substituent_carbon(c_atom, n_idx):
    """Pick carbon substituent on carbonyl carbon (not O, not amide N)."""
    for nb in c_atom.GetNeighbors():
        if nb.GetIdx() == n_idx:
            continue
        if nb.GetAtomicNum() == 8:
            continue
        if nb.GetAtomicNum() == 6:
            return nb.GetIdx()
    return None


def _pick_n_substituent_carbon(n_atom, c_idx):
    """Pick carbon substituent on amide N (not carbonyl carbon)."""
    for nb in n_atom.GetNeighbors():
        if nb.GetIdx() == c_idx:
            continue
        if nb.GetAtomicNum() == 6:
            return nb.GetIdx()
    return None

def _get_carbonyl_oxygen_idx(c_atom):
    mol = c_atom.GetOwningMol()
    c_idx = c_atom.GetIdx()
    for nb in c_atom.GetNeighbors():
        if nb.GetAtomicNum() == 8:
            b = mol.GetBondBetweenAtoms(c_idx, nb.GetIdx())
            if b is not None and b.GetBondType() == Chem.BondType.DOUBLE:
                return nb.GetIdx()
    return None

def _pick_any_h(atom):
    for nb in atom.GetNeighbors():
        if nb.GetAtomicNum() == 1:
            return nb.GetIdx()
    return None


def _count_heavy_neighbors(atom):
    return sum(1 for nb in atom.GetNeighbors() if nb.GetAtomicNum() > 1)

def _collect_amide_omegas(mol):
    """
    Returns list of tuples:
      (r_acyl, c_idx, n_idx, r_n, o_idx)
    omega torsion = r_acyl - C - N - r_n
    """
    out = []
    for bond in mol.GetBonds():
        if bond.GetBondType() != Chem.BondType.SINGLE:
            continue

        a = bond.GetBeginAtom()
        b = bond.GetEndAtom()

        if a.GetAtomicNum() == 6 and b.GetAtomicNum() == 7 and _is_carbonyl_carbon(a):
            c_atom, n_atom = a, b
        elif a.GetAtomicNum() == 7 and b.GetAtomicNum() == 6 and _is_carbonyl_carbon(b):
            c_atom, n_atom = b, a
        else:
            continue

        c_idx = c_atom.GetIdx()
        n_idx = n_atom.GetIdx()
        r_acyl = _pick_acyl_substituent_carbon(c_atom, n_idx)
        r_n = _pick_n_substituent_carbon(n_atom, c_idx)
        o_idx = _get_carbonyl_oxygen_idx(c_atom)

        if r_acyl is None or r_n is None or o_idx is None:
            continue
        out.append((r_acyl, c_idx, n_idx, r_n, o_idx))
    return out

def _make_mmff(mol, confId=0, mmffVariant="MMFF94s"):
    mp = AllChem.MMFFGetMoleculeProperties(mol, mmffVariant=mmffVariant)
    if mp is None:
        raise ValueError("MMFF parameters not available for this molecule; try UFF or fix atom types.")
    ff = AllChem.MMFFGetMoleculeForceField(mol, mp, confId=confId)
    return ff

def _add_trans_omega_constraints(ff, mol, amides, window_deg=5.0, force_constant=2e5):
    """Constrain every amide omega to trans (~ -180°)."""
    target = -180.0
    minD = target - float(window_deg)
    maxD = target + float(window_deg)
    for (r_acyl, c_idx, n_idx, r_n, _o_idx) in amides:
        ff.MMFFAddTorsionConstraint(
            r_acyl, c_idx, n_idx, r_n,
            False, minD, maxD, float(force_constant)
        )
    ff.Initialize()

def _set_all_omegas(mol, confId, amides, omega_deg=-180.0):
    conf = mol.GetConformer(confId)
    for (r_acyl, c_idx, n_idx, r_n, _o_idx) in amides:
        rdMolTransforms.SetDihedralDeg(conf, r_acyl, c_idx, n_idx, r_n, float(omega_deg))

def _rotate_methyl_about_bond(mol, confId, methyl_c_idx, attach_idx, ref_idx, angle_deg):
    """
    Rotate methyl group by setting dihedral:
        H(methyl) - methylC - attach - ref = angle_deg
    Works only if methyl carbon has explicit H.
    """
    conf = mol.GetConformer(confId)
    methyl_c = mol.GetAtomWithIdx(methyl_c_idx)
    h_idx = _pick_any_h(methyl_c)
    if h_idx is None:
        return False
    rdMolTransforms.SetDihedralDeg(conf, h_idx, methyl_c_idx, attach_idx, ref_idx, float(angle_deg))
    return True

def _min_ring_size_containing_atoms(mol, atom_indices):
    """
    Return the smallest ring size that contains all atoms in atom_indices.
    If none, return None.
    """
    ring_info = mol.GetRingInfo()
    if ring_info is None:
        return None
    rings = ring_info.AtomRings()  # tuple of tuples
    want = set(atom_indices)
    best = None
    for r in rings:
        rs = set(r)
        if want.issubset(rs):
            sz = len(r)
            if best is None or sz < best:
                best = sz
    return best

def _is_small_lactam_amide(mol, c_idx, n_idx, max_ring_size=6):
    """
    Lactam heuristic:
      - The C(=O)-N bond is in a ring
      - The smallest ring containing both C and N has size <= max_ring_size
    """
    bond = mol.GetBondBetweenAtoms(int(c_idx), int(n_idx))
    if bond is None or not bond.IsInRing():
        return False

    min_sz = _min_ring_size_containing_atoms(mol, [c_idx, n_idx])
    return (min_sz is not None) and (min_sz <= int(max_ring_size))

def _filter_amides_for_trans_constraints(mol, amides, lactam_max_ring=6):
    """
    Remove small-ring lactams (<= lactam_max_ring) from constraint list.
    """
    kept = []
    skipped = []
    for (r_acyl, c_idx, n_idx, r_n, o_idx) in amides:
        if _is_small_lactam_amide(mol, c_idx, n_idx, max_ring_size=lactam_max_ring):
            skipped.append((r_acyl, c_idx, n_idx, r_n, o_idx))
        else:
            kept.append((r_acyl, c_idx, n_idx, r_n, o_idx))
    return kept, skipped

def optimize_peptide_sdf_mol(
    mol,
    confId=0,
    omega_window_deg=5.0,
    omega_force_constant=2e5,
    maxIts=2000,
    mmffVariant="MMFF94s",
    pre_relax_caps=True,
    no_trans_constraints=False,
    lactam_max_ring=7,
):
    """
    Constrained minimization keeping all peptide/amide omegas trans.
    """
    #mol = Chem.AddHs(mol, addCoords=True)

    if mol.GetNumConformers() == 0:
        raise ValueError("No conformer/3D coords found in the molecule.")

    amides_all = _collect_amide_omegas(mol)

    amides_to_constrain, amides_skipped = _filter_amides_for_trans_constraints(
        mol, amides_all, lactam_max_ring  )

    if not amides_all:
        #raise ValueError("No amide/peptide C(=O)-N bonds detected.")
        print("No amide/peptide C(=O)-N found, will turn off pre_relax_caps",flush=True)
        pre_relax_caps = False
        apply_constraints = False
    elif not amides_to_constrain:
        print(f"Only small-ring lactam amides found (<= {lactam_max_ring}): will NOT apply trans omega constraints.",
              flush=True)
        pre_relax_caps = False
        apply_constraints = False
    else:
        apply_constraints = (not no_trans_constraints)

    if amides_to_constrain:
        _set_all_omegas(mol, confId, amides_to_constrain, omega_deg=-180.0)
    else:
        pre_relax_caps = False

    if pre_relax_caps:
        ace_methyls = []
        nme_methyls = []
        for (r_acyl, c_idx, n_idx, r_n, _o_idx) in amides_to_constrain:
            ra = mol.GetAtomWithIdx(r_acyl)
            rn = mol.GetAtomWithIdx(r_n)
            if ra.GetAtomicNum() == 6 and _count_heavy_neighbors(ra) == 1:
                ace_methyls.append((r_acyl, c_idx, n_idx))  # rotate around methyl-C(=O), ref=N
            if rn.GetAtomicNum() == 6 and _count_heavy_neighbors(rn) == 1:
                nme_methyls.append((r_n, n_idx, c_idx))     # rotate around methyl-N, ref=C

        ace_methyls = list({t for t in ace_methyls})
        nme_methyls = list({t for t in nme_methyls})

        grid = [0, 60, 120, 180, 240, 300]
        best = copy.deepcopy(mol)
        bestE = None

        ace_opts = ace_methyls if ace_methyls else [(None, None, None)]
        nme_opts = nme_methyls if nme_methyls else [(None, None, None)]

        for a_ang in grid:
            for n_ang in grid:
                trial = copy.deepcopy(mol)

                for (methyl_c, carbonyl_c, ref_n) in ace_opts:
                    if methyl_c is not None:
                        _rotate_methyl_about_bond(trial, confId, methyl_c, carbonyl_c, ref_n, a_ang)

                for (methyl_c, attach_n, ref_c) in nme_opts:
                    if methyl_c is not None:
                        _rotate_methyl_about_bond(trial, confId, methyl_c, attach_n, ref_c, n_ang)

                ff = _make_mmff(trial, confId=confId, mmffVariant=mmffVariant)
                _add_trans_omega_constraints(
                    ff, trial, _collect_amide_omegas(trial),
                    window_deg=omega_window_deg,
                    force_constant=omega_force_constant
                )
                ff.Minimize(maxIts=1000)
                E = ff.CalcEnergy()

                if bestE is None or E < bestE:
                    bestE = E
                    best = trial

        mol = best

    ff = _make_mmff(mol, confId=confId, mmffVariant=mmffVariant)
    if apply_constraints and amides_to_constrain:
        _add_trans_omega_constraints(
        ff, mol, amides_to_constrain,
        window_deg=omega_window_deg,
        force_constant=omega_force_constant
    )
    ff.Minimize(maxIts=int(maxIts))
    return mol


def _read_first_sdf(path: str, no_sanitize) -> Chem.Mol:
    flag = False if no_sanitize else True
    suppl = Chem.SDMolSupplier(path, removeHs=False,sanitize=flag)
    if not suppl or len(suppl) == 0 or suppl[0] is None:
        raise ValueError(f"Failed to read SDF (first molecule): {path}")
    return suppl[0]

def _write_sdf(mol: Chem.Mol, path: str):
    os.makedirs(os.path.dirname(os.path.abspath(path)) or ".", exist_ok=True)
    w = Chem.SDWriter(path)
    w.write(mol)
    w.close()


def mapping_atom_back_sdf_to_pdb(template_pdb: str, opt_sdf: str,
                                 out_pdb: str, sanitize_flag, confId: int = 0):
    """
    Index-based remap:
      - preserves residue/atom names & all non-ATOM records from template_pdb
      - replaces ONLY XYZ fields (cols 30:54) using coordinates from opt_sdf
    No coordinate tolerance is used/needed (mapping is by atom order).

    REQUIREMENT:
      #ATOM/HETATM lines in template_pdb == number of atoms in opt_sdf
      (including H, if present in either file).
    """
    atom_rec = ("ATOM  ", "HETATM")

    pdb_lines = open(template_pdb, "r").readlines()
    atom_line_indices = [k for k, ln in enumerate(pdb_lines) if ln[0:6] in atom_rec]
    if not atom_line_indices:
        raise RuntimeError(f"FATAL: no ATOM/HETATM records in {template_pdb}")

    mol = _read_first_sdf(opt_sdf,sanitize_flag)
    if mol.GetNumConformers() == 0:
        raise RuntimeError(f"FATAL: no conformer/coords found in {opt_sdf}")
    conf = mol.GetConformer(confId)

    n_pdb = len(atom_line_indices)
    n_sdf = mol.GetNumAtoms()
    if n_pdb != n_sdf:
        raise RuntimeError(
            f"FATAL: atom count mismatch: PDB has {n_pdb} ATOM/HETATM lines, "
            f"SDF has {n_sdf} atoms. Index-based mapping cannot proceed.\n"
            f"Tip: make sure both include/exclude H consistently (removeHs=False)."
        )

    out_lines = []
    ai = 0
    for ln in pdb_lines:
        if ln[0:6] not in atom_rec:
            out_lines.append(ln)
            continue

        p = conf.GetAtomPosition(ai)
        newxyz = f"{p.x:8.3f}{p.y:8.3f}{p.z:8.3f}"

        ln2 = ln.rstrip("\n")
        if len(ln2) < 54:
            ln2 = ln2.ljust(54)
        ln2 = ln2[:30] + newxyz + ln2[54:]
        out_lines.append(ln2 + "\n")
        ai += 1

    os.makedirs(os.path.dirname(os.path.abspath(out_pdb)) or ".", exist_ok=True)
    with open(out_pdb, "w") as f:
        f.writelines(out_lines)

    return out_lines

def run_mm_opt(
    input_sdf: str,
    output_sdf: str,
    template_pdb: str = None,
    output_pdb: str = None,
    rdkit_sdf: str = None,
    no_trans_constraints: bool = False,
    no_sanitize: bool = False,
    omega_window_deg: float = 5.0,
    omega_force_constant: float = 2e5,
    maxIts: int = 2000,
    mmffVariant: str = "MMFF94s",
    no_pre_relax_caps: bool = False,
):
    """
    Run RDKit MM optimization and optional index-based coordinate mapping back to PDB.
    """
    # Decide RDKit intermediate filename.
    if rdkit_sdf is None:
        base, _ = os.path.splitext(os.path.abspath(output_sdf))
        rdkit_sdf = base + ".rdkit_opt.tmp.sdf"

    # Skip peptide/amide optimization for pure ligand SDFs when a template PDB is given.
    opt_with_amide = True
    if "pure_lig" in input_sdf and template_pdb:
        opt_with_amide = False

    mol = _read_first_sdf(input_sdf, no_sanitize)

    if opt_with_amide:
        mol_opt = optimize_peptide_sdf_mol(
            mol,
            omega_window_deg=omega_window_deg,
            omega_force_constant=omega_force_constant,
            maxIts=maxIts,
            mmffVariant=mmffVariant,
            pre_relax_caps=(not no_pre_relax_caps),
            no_trans_constraints=no_trans_constraints,
        )
        _write_sdf(mol_opt, rdkit_sdf)
    else:
        rdkit_sdf = input_sdf

    if os.path.abspath(rdkit_sdf) != os.path.abspath(output_sdf):
        shutil.copy2(rdkit_sdf, output_sdf)

    final_sdf = output_sdf

    # Optional mapping back to PDB.
    if output_pdb:
        if not template_pdb:
            raise ValueError(
                "If output_pdb is set, you must provide template_pdb."
            )
        mapping_atom_back_sdf_to_pdb(
            template_pdb,
            final_sdf,
            output_pdb,
            no_sanitize,
        )

    print(f"[OK] RDKit intermediate: {rdkit_sdf}")
    print(f"[OK] Final SDF:         {final_sdf}")
    if output_pdb:
        print(f"[OK] Output PDB:        {output_pdb}")

    # Return paths so an importing script can use them programmatically.
    return {
        "rdkit_sdf": rdkit_sdf,
        "final_sdf": final_sdf,
        "output_pdb": output_pdb,
    }

def main():
    ap = argparse.ArgumentParser(
        description="Two-step optimization (RDKit trans-omega + optional xTB) and index-based mapping back to PDB."
    )
    ap.add_argument("-isdf", "--input_sdf", required=True, help="Input SDF with 3D coords, required")
    ap.add_argument("-osdf", "--output_sdf", required=True, help="Final optimized SDF (RDKit-only or xTB result), required")

    ap.add_argument("--rdkit_sdf", default=None, help="Intermediate RDKit SDF (default: <output>.two_step_opt.sdf)")
    ap.add_argument("--no_trans_constraints",action="store_true",default=False,
                    help="Do not constrain all amide bonds to have trans conf, default=False"
                        "use this flag if your have amide inside a small ring!")
    ap.add_argument("--no_sanitize", default=False, action="store_true",
                    help="Do sanitize the SDF when reading it, default=True")
    ap.add_argument("--omega_window_deg", type=float, default=5.0,help="peptide bond omega sample window, default=5.0")
    ap.add_argument("--omega_force_constant", type=float, default=2e5,help="peptide bond constrain constance, default=2e5")
    ap.add_argument("--maxIts", type=int, default=2000,help="max iteractions for rdkit FF optimizer, default=2000")
    ap.add_argument("--mmffVariant", default="MMFF94s")
    ap.add_argument("--no_pre_relax_caps", action="store_true", help="Disable ACE/NME methyl pre-rotation")
    # mapping back
    ap.add_argument("-ipdb", "--template_pdb", default=None, help="Template PDB with residue/atom naming")
    ap.add_argument("-opdb", "--output_pdb", default=None, help="Output PDB with optimized coordinates\n")
    args = ap.parse_args()

    return run_mm_opt(
        input_sdf=args.input_sdf,
        output_sdf=args.output_sdf,
        template_pdb=args.template_pdb,
        output_pdb=args.output_pdb,
        rdkit_sdf=args.rdkit_sdf,
        no_trans_constraints=args.no_trans_constraints,
        no_sanitize=args.no_sanitize,
        omega_window_deg=args.omega_window_deg,
        omega_force_constant=args.omega_force_constant,
        maxIts=args.maxIts,
        mmffVariant=args.mmffVariant,
        no_pre_relax_caps=args.no_pre_relax_caps,
    )

if __name__ == "__main__":
    main()
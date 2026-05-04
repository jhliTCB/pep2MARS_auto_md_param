#!/public/home/lijunhao/soft/conda_envs/fldev/bin/python
# protocol to prepare capped blocks for NCAA peptides
# TODO: Incorparate DNA/RNA capping procedure
  #### use neighboring CONH as caps
import json
import os
import re
import shutil
import subprocess
import sys
import time
import argparse
import zipfile, string
from collections import defaultdict, Counter
from utils.blocks_builder import (unitClassifier,known_residues,nucRes)
from utils.mm_opt import run_mm_opt
from pymol import cmd, editor
from pymol import CmdException
from rdkit import Chem
import MDAnalysis as mda
from MDAnalysis.topology.tables import vdwradii

topDir = os.path.dirname(os.path.abspath(__file__))
def zipDir(dirpath, outFullName):
    """
    :param dirpath: target directory path
    :param outFullName: path+xxxx.zip
    :return: None
    """
    zip = zipfile.ZipFile(outFullName, "w", zipfile.ZIP_DEFLATED)
    for path, dirnames, filenames in os.walk(dirpath):
        # remove full path name
        fpath = path.replace(dirpath, '')
        for filename in filenames:
            zip.write(os.path.join(path, filename), os.path.join(fpath, filename))
    zip.close()

def rename_residues_in_pymol(rename_map):
    '''rename_map format:
          { (old_resn, chain, resi): (new_resn, chain, resi), ... }
    '''
    if len(rename_map) == 0: return rename_map
    for (old_resn, chain, resi), (new_resn, _, _) in rename_map.items():
        if new_resn != old_resn:
            cmd.alter(f"r. {old_resn} & c. {chain} & i. {resi}", f"resn='{new_resn}'")
    return rename_map

def split_nonprotein(pdb, output_dir="./", renamed_pdb=None):
    """
    - split ncaa or ligand into standalone PDB files
    - residues with same names, the less atom residue will be renamed, e.g. ZDC -> ZDc
    - rename ZDC -> ZDc in renamed_pdb
    - return nonProt list and all_keys list
    depands on global values of three2one, nucRes, ionsRes, watRes
    """
    with open(pdb, 'r') as f:
        lines = f.readlines()

    # 1) collects ncaa and ligand,  key: resn:chain:resi
    residue_data = defaultdict(list)   # list for all ncaa and ligands
    residue_atoms = defaultdict(set)
    for line in lines:
        if line.startswith(("ATOM  ", "HETATM")):
            residue_name = line[17:20].strip()
            if residue_name in known_residues:
                continue
            chain_name   = line[21:22].strip()
            residue_num  = line[22:26].strip()
            atom_name    = line[12:16].strip()
            res_key      = f"{residue_name}:{chain_name}:{residue_num}"
            residue_data[res_key].append(line)
            residue_atoms[residue_name].add(atom_name)

    # 2) if there are same-name residues differ in atom number
    res_instances = defaultdict(list)  # resn -> [(res_key, atom_count), ...]
    for res_key, block in residue_data.items():
        resn, ch, ri = res_key.split(':')
        atoms = {L[12:16].strip() for L in block if L.startswith(("ATOM", "HETATM"))}
        res_instances[resn].append((res_key, len(atoms)))

    # 3) Find out less-atom residue, rename it: e.g. ZDC -> ZDc
    rename_map = {}                    # res_key -> new_resn（e.g. ZDc）
    for resn, insts in res_instances.items():
        if len(insts) <= 1:
            continue
        base_count = max(c for _, c in insts)
        if len(resn) == 3:
            lower3 = resn[:2] + resn[2].lower()
        else:
            # non 3-letter residue, keep it the same name
            continue
        for res_key, cnt in insts:
            if cnt < base_count:
                rename_map[res_key] = lower3

    os.makedirs(output_dir, exist_ok=True)

    # 4) write out each nacc/ligand pdb
    nonProt, allNPs = [], []
    kf = open('ncaa_ligs.keys', 'w') # for pymol selection
    for res_key, block in residue_data.items():
        resn, ch, ri = [s.strip() for s in res_key.split(':')]
        out_resn = rename_map.get(res_key, resn)
        if out_resn not in nonProt:
            nonProt.append(out_resn)
        output_file = os.path.join(output_dir, f"{out_resn}.pdb")
        with open(output_file, 'w') as f:
            for L in block:
                if L.startswith(("ATOM", "HETATM")):
                    # replace resn in clolumns 18-20 
                    f.write(L[:17] + out_resn.ljust(3) + L[20:])
                else:
                    f.write(L)
            f.write("END\n")
        allNPs.append([out_resn,ch,ri])
        kf.write(f"{out_resn},{ch},{ri}\n")
    kf.close()

    # 5) get renamed_pdb for the duplicated resn for differed residues
    if renamed_pdb is None:
        base = os.path.basename(pdb)
        renamed_pdb = os.path.join(output_dir, f".renamed.{pdb}")

    # constructing a set for comparison: { (resn, chain, resi) } -> out_resn
    inst_renames = {}
    for res_key, new_resn in rename_map.items():
        r, c, i = [s.strip() for s in res_key.split(':')]
        inst_renames[(r, c, i)] = new_resn
    
    with open(renamed_pdb, 'w') as g:
        for L in lines:
            if L.startswith(("ATOM  ", "HETATM")):
                r = L[17:20].strip()
                c = L[21:22].strip()
                i = L[22:26].strip()
                newr = inst_renames.get((r, c, i))
                if newr:
                    g.write(L[:17] + newr.ljust(3) + L[20:])
                else:
                    g.write(L)
            else:
                g.write(L)

    # 6) get the list of all NCAAs and LIGs
    with open(os.path.join(output_dir, "ligs.ncaas.resn.lst"), "w") as fp:
        fp.write('\n'.join(nonProt))

    return nonProt, allNPs

def add_remark(data, pdbF):
    '''
    example data:
    [
    ( (('GLy', 'J', '8'), 'N'), (('LNB', 'J', '5'), 'C3') ),
    ( (('HIp', 'J', '6'), 'N'), (('LNB', 'J', '5'), 'C1') ),
    ( (('LNB', 'J', '5'), 'N1'), (('THr', 'J', '4'), 'C') )]
    '''
    remarks = []
    with open(pdbF, 'r') as f:
        pdbFL = f.read().splitlines()
    for edge in data:
        (set1, set2) = edge
        resn1, atom1 = set1[0][0], set1[1]
        resn2, atom2 = set2[0][0], set2[1]
        line = f"REMARK: bond cpx.{resn1}.{atom1} cpx.{resn2}.{atom2}"
        if line not in remarks: remarks.append(line)
    new_pdb = remarks + pdbFL
    with open(pdbF, 'w') as f:
        f.writelines('\n'.join(new_pdb))
    return remarks

def find_linking(selection="r. MOL & c. A & i. 777"):
    '''
    selection of a residue/molecule, uses the pymol's bonding and connection table
    '''
    link_dict = {}
    for atom in cmd.get_model(selection).atom:
        itself0 = f'(name "{atom.name}" & r. {atom.resn} & i. {atom.resi})'
        # exclude hydrogen for case of close hydrogen atom in HB
        extend1 = f'{selection} & name "{atom.name}" extend 1 & ! {itself0} & ! elem H'
        cmd.select("tmp_atm",extend1)
        itself1 = f"r. {atom.resn} & c. {atom.chain} & i. {atom.resi}"
        for x in cmd.get_model("tmp_atm").atom:
            sele2 = f"r. {x.resn} & c. {x.chain} & i. {x.resi}"
            if sele2 != itself1:
                link_dict[f'{selection} & name "{atom.name}"'] = f'{sele2} & name "{x.name}"'
                break
    return link_dict

def editor_run(sele,cap):
    '''cap == ace or nme'''
    '''pitfall: when editor command failed, it will create non-necessary states'''
    cmd.select("extend_one",f"{sele} extend 1")
    # termini should have 3 atoms, including itself
    # if it was capped with HXT, make sure if you need it
    if cmd.count_atoms("extend_one") != 3:
        print(f"the connection number for {sele} is not 2",flush=True)
        cmd.delete("extend_one")
        return 0
    cmd.delete("extend_one")
    try:
        editor.attach_amino_acid(sele, cap, ss=1)
    except CmdException:
        print(f"Error for adding {cap} to {sele}: {CmdException}",flush=True)
        return 0
    else:
        print(f"{cap} cap added to {sele} successfully",flush=True)
        return 1

def attach_caps_for_res(obj, parent_obj, resN, resC):
    """add ACE/NME caps to the N/C atoms of the given obj
    could be odd number of caps"""
    capped_num = 0
    if len(resN) == len(resC):
        for capN, capC in zip(resN,resC):
            NList = [obj] + list(capN)
            CList = [obj] + list(capC)
            if len(NList) > 1:
                selN = "{} and r. {} and c. {} and i. {} and name N".format(*NList)
                if cmd.count_atoms(selN) == 1:
                    capped_num += editor_run(selN, 'ace')
            if len(CList) > 1:
                selC = "{} and r. {} and c. {} and i. {} and name C".format(*CList)
                if cmd.count_atoms(selC) == 1:
                    capped_num += editor_run(selC,'nme')
    else:
        if len(resN) == 1 and len(resC) == 0:
            resList = [obj] + list(resN[0]) + ['N']
            cap_type = 'ace'
        if len(resC) == 1 and len(resN) == 0:
            resList = [obj] + list(resC[0]) + ['C']
            cap_type = 'nme'
        sel = "{} and r. {} and c. {} and i. {} and name {}".format(*resList)
        capped_num += editor_run(sel, cap_type)

    return capped_num

def get_sele_string(unit,keys):
    '''get pymol selection key for a given unit/block'''
    seleStr = ""
    for key in keys:
        if key in ['a','b']:
            kl = list(unit[key])
            kl.append(' ') if key == "b" else kl.append('|')
            seleStr += "(r. {} & c. {} & i. {}) {}".format(*kl)
        else:
            for i in range(len(unit[key])):
                kl = list(unit[key][i])
                kl.append(' ') if i == len(unit[key])-1 else kl.append('|')
                seleStr += "(r. {} & c. {} & i. {}) {}".format(*kl)
    return seleStr

def NC_swapping(atom):
    '''always N first, C last'''
    _resN, _resC = [], []
    if atom.name == 'N':
        _resN = [tuple([atom.resn, atom.chain, atom.resi])]
        _resC = []
    elif atom.name == 'C':
        _resN = []
        _resC = [tuple([atom.resn, atom.chain, atom.resi])]
    return _resN, _resC

def create_with_caps(sele,obj,capping,color_id=0):
    ''' sele: the block selection string
    obj: the parent object, to check if there is linking not belong to sele's (resn+resi+chain)
    '''
    colors = [5,154,6,9,144,11,10,5262,12,36,5271,124]*5
    _sele = cmd.get_model(sele).atom
    res0 = []
    for x in _sele:
        if (x.resn,x.chain,x.resi) not in res0:
            res0.append([x.resn,x.chain,x.resi])
    # possible better way: get the keys of "C" and "N" atom names!
    extend1 = cmd.get_model(f"({sele} extend 1 & {obj}) & ! {sele}").atom
    if len(extend1) == 0:
        cmd.create(f"capped_{sele}", sele)
        attach_caps_for_res(obj=f"capped_{sele}", parent_obj=obj,
                            resN=capping, resC=capping)
        return "", []
    else:
        extended_keys = [f"r. {extend1[i].resn} & c. {extend1[i].chain} & i. {extend1[i].resi}"
                         for i in range(len(extend1))]
        extended_name = [x.name for x in extend1]
        counter = Counter(extended_keys)
        duplicates = {item for item, count in counter.items() if count > 1}
        x, x_name = [], [] # x is the non duplicated ress list of the extend1 residues, y is the duplicated ress list
        for item, ext in zip(extended_keys, extended_name):
            if item not in duplicates:
                x.append(item)
                x_name.append(ext)
        y = [item for item, count in counter.items() if count > 1]

        extended_sele = []
        for i in range(len(x)):
            keep = cmd.get_model(f"{x[i]} & {obj}").atom
            #_key = f"r. {x[i].resn} & c. {x[i].chain} & i. {x[i].resi}"
            _key = x[i]
            if x_name[i] == "N": # connecting to block's "C", alters to NME
                extended_sele.append(f"({_key} & name N+H+CA+CB)")
            elif x_name[i] == "C": # connecting to block's "N", alters to ACE
                #cmd.select(f"tmp_ace_{i}", f"{_key} & {obj} & name C+O+CA")
                extended_sele.append(f"({_key} & name C+O+CA+CB)")
        capped_sele = f"({sele} + {'+'.join(extended_sele)}) & {obj}"
        cmd.select(f"extended_{sele}",capped_sele)
        cmd.create(f"capped_{sele}",capped_sele)
        for i in range(len(extended_sele)):
            '''Good to know that empty selections did not show errors 
                in cmd.select, cmd.valence, and cmd.remove!'''
            cmd.valence(1, f"capped_{sele} & {extended_sele[i]} & name CA",
                        f"capped_{sele} & {extended_sele[i]} & name CB")
            cmd.remove(f"capped_{sele} & {extended_sele[i]} & name CB")
            # extended_sele[i] is "(r. PRO & c. A & i. 2 & name C+O+CA+CB)"
            just_the_key = f"{'&'.join(extended_sele[i].split('&')[0:-1])})"
            if 'N+H+CA' in extended_sele[i]:
                #it is possible that CA does not exist
                to_be_handle = ['NME','N','C'] # resn, capped place, CH3's PDB name
            if 'C+O+CA' in extended_sele[i]:
                to_be_handle = ['ACE','C','CH3'] # resn, capped place, CH3's PDB name
                # some PDB's C=O valence is 1, set it to 2
                cmd.valence(2, f"capped_{sele} & {extended_sele[i]} & name C",
                            f"capped_{sele} & {extended_sele[i]} & name O")
            ## handle together
            tmp_sel = cmd.get_model(f"capped_{sele} & {just_the_key}").atom
            dist_sel1 = f"capped_{sele} & {just_the_key} & name {to_be_handle[1]}"
            dist_sel2 = f"capped_{sele} & {just_the_key} & name CA"
            if "CA" not in [xx.name for xx in tmp_sel]:
                cmd.edit(f"capped_{sele} & {just_the_key} & name {to_be_handle[1]}")
                cmd.attach("C", 1, 3, to_be_handle[2])
                cmd.edit()
            elif cmd.distance(dist_sel1, dist_sel2) > 1.4: # type: ignore
                print(f"distance of\n{dist_sel1}\n{dist_sel2} > 1.4 Angstrom\n, Removing {dist_sel2}")
                cmd.remove(f"capped_{sele} & {just_the_key} & name CA")
                cmd.edit(f"capped_{sele} & {just_the_key} & name {to_be_handle[1]}")
                cmd.attach("C", 1, 3, to_be_handle[2])
                cmd.edit()
            cmd.alter(f"capped_{sele} & {just_the_key}", f"resn='{to_be_handle[0]}'")
        cmd.h_add(f"capped_{sele} & r. ACE+NME")
        cmd.alter(f"capped_{sele} & r. NME & elem N extend 1 & elem H", f"name='HNN'")
        cmd.alter(f"capped_{sele} & r. ACE & elem O extend 1 & elem C", f"name='C'")
        cmd.alter(f"capped_{sele} & r. ACE & elem H extend 1 & elem C", f"name='CA'")

        for i in range(len(y)):
            overlaps = cmd.get_model(f"{y[i]} extend 1 & {sele} & {obj}").atom
            for o in overlaps:
                _resN, _resC = NC_swapping(o)
                attach_caps_for_res(obj=f"capped_{sele}", parent_obj=obj,
                                        resN=_resN,resC=_resC)

        for NC in cmd.get_model(f"({sele} & name C+N) extend 1 & {obj}").atom:
            express = f"((r. {NC.resn} & c. {NC.chain} & i. {NC.resi} & name {NC.name} & {obj}) extend 1) & ! {sele}"
            if len(cmd.get_model(express).atom) == 0:
                _resN, _resC = NC_swapping(NC)
                attach_caps_for_res(obj=f"capped_{sele}", parent_obj=obj,
                                    resN=_resN,resC=_resC)
    cmd.util.cba(colors[color_id], f"capped_{sele}", _self=cmd)
    return f"capped_{sele}", extended_sele

def get_blocks(root_obj,pdbf, capping_method, nonProt,allNPs):
    '''
    Previously: allNPs.append([out_resn,ch,ri])
    link_dict's lenght:
    0 = standard ligand, export 1 resn as it is
    1 = ligand with a residue, export 2 resn with ACE+NME
    2 = a residue without ligand, export 1 resn with ACE+NME
        when linking != C/N, could be a ligand's part, export 3 resn with 2 ACE+NME!
    3 = a residue for cyclic connection, the cases need special cares:
        export 2 or 3 resn with 2 ACE+NME
    4 : needs to be avoided when preparing the structure, however, if there are common
        residues non-backbone connecting a pair and a triple nodes, quartet-residue block occurs
    '''
    cmd.load(pdbf, root_obj)
    cmd.set("retain_order", 1)
    cmd.set("pdb_conect_all", 1)
    link_list, blocks_list = [], []
    ligand_id = 0
    for res in nonProt:
        for item in allNPs:
            if item[0] == res:
                sele = "r. {} & c. {} & i. {}".format(*item)
                link_dict = find_linking(sele)
                print(f"the link_dict for {res} is:\n{link_dict}",flush=True)
                link_list.append(link_dict)
                if len(link_dict) == 0:
                    print(f"no linking atom found for {res}, assuming it is a pure ligand",flush=True)
                    # {res}.pdb should be there after split_nonprotein
                    cmd.select(f"pure_lig_{ligand_id}_{res}",sele)
                    cmd.save(f"pure_lig_{ligand_id}_{res}.pdb",f"pure_lig_{ligand_id}_{res}")
                    cmd.save(f"pure_lig_{ligand_id}_{res}.sdf", f"pure_lig_{ligand_id}_{res}")
                    blocks_list.append(f"pure_lig_{ligand_id}_{res}.pdb")
                    
    print(f"link_list = {link_list}",flush=True)
    clear = unitClassifier()
    units = clear.classify(link_list)
    rename_residues_in_pymol(units['rename_map'])
    color_id = 0
    for block in ['quartets','triples','pairs','singles']:
        print(f"{block} in units is:\n {units[block]}\n",flush=True)
        is_nucRes = False
        if block in ['triples','quartets']:
            sele_keys = ['nodes']
        elif block == 'pairs':
            sele_keys = ['a','b']
        if len(units[block]) > 0 and block != 'singles':
            for i in range(len(units[block])):
                outPdb = f"capped_{block[:-1]}_unit{i}.pdb"
                outSdf = f"capped_{block[:-1]}_unit{i}.sdf"
                seleStr = get_sele_string(units[block][i],sele_keys)
                selName = f"{block[:-1]}_unit{i}"
                capping = []
                for e in units[block][i]['edges']:
                    for ee in e:
                        if ee[0][0] in nucRes:
                            is_nucRes = True
                        if ee[1] in ['N','C']:
                            capping.append(ee[0])
                if is_nucRes:
                    print(f"nuclear acid residues found: {units[block][i]['edges']}, skipping...")
                    print("check this tutorial: https://amberhub.chpc.utah.edu/create-new-nucleotide-for-amber/")
                    continue
                cmd.select(selName,seleStr)
                capped_name, extended_sele = create_with_caps(selName, root_obj, capping, color_id)
                color_id += 1
                cmd.save(outPdb, capped_name)
                cmd.save(outSdf, capped_name)
                add_remark(units[block][i]['edges'], outPdb)
                with open(f"ter_{outPdb}", 'w') as fp:
                    fp.write(''.join(add_het_ter(outPdb)))
                os.system(f"cp -vf ter_{outPdb} {outPdb} ")
                blocks_list.append(outPdb)

        if len(units[block]) > 0 and block == 'singles':
            for unit in units['singles']:
                unit = list(unit)
                if unit[0] in nucRes:
                    print(f"nuclear acid residues found: {units[block][i]['edges']}, skipping...")
                    print("check this tutorial: https://amberhub.chpc.utah.edu/create-new-nucleotide-for-amber/")
                    continue
                cmd.select(f"single_{unit[0]}", "r. {} & c. {} & i. {}".format(*unit))
                capped_name, extended_sele = create_with_caps(f"single_{unit[0]}", root_obj, unit[0], color_id)
                color_id += 1
                cmd.save(f"{capped_name}.pdb", capped_name)
                cmd.save(f"{capped_name}.sdf", capped_name)
                blocks_list.append(f"{capped_name}.pdb")
    orienting_center = "capped_*"
    # simplified logics now
    #TODO: non-standard capping is still not well supported yet!
    if capping_method.lower() in ['non-standard','non_standard','nonstandard']:
        nsc = 'nonStd_cAAPEd'
        os.makedirs(nsc, exist_ok=True)
        blocks_list = [x for x in blocks_list if x.startswith('pure_lig')]
        std_capping = cmd.get_names(selection='capped_*')
        for nsc_block in std_capping:
            new_AA_list = []
            for _res in cmd.get_model(nsc_block).atom:
                if _res.resn not in ['ACE','NME'] and _res.resn not in new_AA_list:
                    new_AA_list.append(_res.resn)
            for _res in new_AA_list:
                cmd.select(f"NonStd_{_res}",f"{nsc_block} & r. {_res}")
                cmd.select(f"NonStd_{_res}_ext2", f"(NonStd_{_res} extend 2) & {nsc_block}")
                cmd.create(f"{nsc}_{_res}_ext2", f"NonStd_{_res}_ext2")
                t_sele = f"(({nsc}_{_res}_ext2 & r. {_res}) extend 2) & ! elem H & ! " + \
                         f"(({nsc}_{_res}_ext2 & r. {_res}) extend 1)"
                cmd.select("terminis",t_sele)
                for atm in cmd.get_model("terminis").atom:
                    cmd.edit(f"r. {atm.resn} & i. {atm.resi} & name {atm.name} & {nsc}_{_res}_ext2")
                    cmd.h_fill()
                    cmd.edit()
                cmd.save(f"{nsc}/cAPPEd_{_res}.pdb", f"{nsc}_{_res}_ext2")
                cmd.save(f"{nsc}/cAPPEd_{_res}.sdf", f"{nsc}_{_res}_ext2")
                blocks_list.append(f"{nsc}/cAPPEd_{_res}.pdb")
        cmd.disable(orienting_center)
        orienting_center = f"{nsc}_*"
    # save a pymol session when everything is done:
    cmd.show("lines", "all")
    cmd.show("sticks", "capped*")
    cmd.create(f"{root_obj}_noh",root_obj)
    cmd.select("deH", f"elem H & {root_obj}_noh")
    cmd.remove("deH")
    cmd.delete("deH")
    cmd.set("stick_radius",0.15) # type: ignore
    cmd.set("cartoon_oval_length", 0.8) # type: ignore
    cmd.set("cartoon_oval_width", 0.12) # type: ignore
    cmd.set("cartoon_rect_length", 0.8) # type: ignore
    cmd.set("cartoon_rect_width", 0.12) # type: ignore
    cmd.set("cartoon_loop_radius", 0.12) # type: ignore
    cmd.set("stick_h_scale", 1.0) # type: ignore
    cmd.set("cartoon_gap_cutoff", 0)
    cmd.bg_color("white")
    cmd.set("ray_opaque_background", 1)
    cmd.set("use_shaders", 1)
    cmd.save(f"{root_obj}_nh_leap.pdb", f"{root_obj}_noh")
    cmd.hide("cartoon", "capped*")
    cmd.delete('dist*')
    cmd.disable("full_st*")
    cmd.select("CAPs","r. ACE+NME")
    cmd.hide("sticks","CAPs")
    cmd.hide("spheres")
    cmd.set("line_width",2)
    cmd.orient(orienting_center)
    #cmd.save(f"{root_obj}_nh_leap.mol2",f"{root_obj}_noh")
    with open(f"{root_obj}_nh_leap_ter.pdb",'w') as fp:
        fp.write(''.join(add_het_ter(f"{root_obj}_nh_leap.pdb")))

    cmd.save("capping_results.pse")
    cmd.reinitialize()
    return blocks_list

def guessing_charge(blocks_list,guess_method,charge_file):
    '''
    pdbf2 is the connection corrected original PDB files for selections
    blocks_list is a pdb file list of blocks:
        capped_triple_unit0.pdb, capped_pair_unit0.pdb, capped_single_0A1.pdb, etc...
    if standard capping is used, charge_dict is like:
        {'capped_triple_unit0.pdb': 0, 'capped_pair_unit0.pdb': 0, 'capped_single_0A1.pdb': 0}
    else if non-standard capping is used, charge_dict is like:
        {'cAPPEd_89N.pdb': 0, 'cAPPEd_Z50.pdb': 0, 'cAPPEd_0A1.pdb': 0, 'cAPPEd_DLY.pdb': 1}
    '''
    charge_dict = {}
    if os.path.isfile(charge_file):
        print(f"{charge_file} file found from previous calculation, will read charges from it",flush=True)
        charge_dict = json.load(open(charge_file,'r',encoding="utf-8"))
        print(f"charge of each block:\n{charge_dict}\n",flush=True)
    else:
        for i in range(len(blocks_list)):
            zip_flag = True if i == len(blocks_list) - 1 else False
            pdb = blocks_list[i]
            charge_dict[pdb] = guess_chg_run(pdb,zip_flag,guess_method)
        json.dump(charge_dict, open(charge_file,"w"), indent=4)

    return charge_dict

def guess_chg_run(pdb,zip_flag,guess_method='mm'):
    '''
    choices: Gaussian guess from a range or with correct connections for rdkit to derive formal charge!
    '''
    if zipfile.is_zipfile('chgG.zip'):
        zip = zipfile.ZipFile('chgG.zip')
        zip.extractall('./chgG')
    else:
        os.makedirs('chgG',exist_ok=True)
    os.chdir('chgG')
    basename = pdb[:-4]
    base = basename.split('/')[-1]
    # make sure that the sdf's connection is correct!
    #if 'single' in pdb:     no_trans_const = '--no_trans_constraints --no_pre_relax_caps'
    no_trans_const = ' '
    run_mm_opt(
    input_sdf=f"../{basename}.sdf",
    output_sdf=f"{base}.sdf",
    template_pdb=f"../{pdb}",
    output_pdb=f"{base}.pdb",
    maxIts=3000,
    no_trans_constraints=no_trans_const,
    )
    
    if 'pure_lig' in f'{base}.pdb':
        os.system(f"cp -vf ../{pdb} .")
        os.system(f"cp -vf ../{basename}.sdf .")
    if guess_method in ['qm','QM','Gaussian']:
        min_list = qm_chg_run(f'{base}.pdb') # ('jobid',[chg, logF, ener])
        guessChg = min_list[1][0]
        if zip_flag: os.system("/usr/bin/rm -f *.chk slurm*")
    elif guess_method in ['mm','MM','rdkit']:
        guessChg = mm_chg_run(f"{base}.sdf")
    os.chdir('..')

    if zip_flag:
        zipDir("./chgG","./chgG.zip")
        shutil.rmtree("./chgG")

    return guessChg

def mm_chg_run(sdf_path):
    '''
    May not good to use for some special cases
    '''
    suppl = Chem.SDMolSupplier(sdf_path, removeHs=False, sanitize=True)
    if len(suppl) == 0:
        raise ValueError(f"No records found in SDF: {sdf_path}")
    mol = suppl[0]
    if mol is None:
        #raise ValueError(f"Could not parse molecule at index 0 from {sdf_path}")
        print(f"Could not parse molecule at index 0 from {sdf_path}, Turn off Sanitize!")
        suppl = Chem.SDMolSupplier(sdf_path, removeHs=False, sanitize=False)
        mol = fix_phosphates(suppl[0])
    # Closed-shell check (no unpaired electrons)
    radical_e = sum(a.GetNumRadicalElectrons() for a in mol.GetAtoms())
    if radical_e != 0:
        raise ValueError(f"Not closed-shell: found {radical_e} radical electrons")
    # Net formal charge
    return sum(a.GetFormalCharge() for a in mol.GetAtoms())

def fix_phosphates(mol):
    '''
    Phosphates systems, bonding information is prone to error, need to fix it
    '''
    mol = Chem.RWMol(mol)

    for p in [a for a in mol.GetAtoms() if a.GetSymbol() == "P"]:
        o_neis = [n for n in p.GetNeighbors() if n.GetSymbol() == "O"]
        if len(o_neis) < 3:
            continue

        # terminal oxygens: O bonded only to P (degree 1) OR bonded to P and H (degree 2 with H)
        terminal = []
        bridging = []
        for o in o_neis:
            # count heavy neighbors (non-H)
            heavy_nbrs = [n for n in o.GetNeighbors() if n.GetSymbol() != "H"]
            if len(heavy_nbrs) == 1:
                terminal.append(o)
            else:
                bridging.append(o)

        if not terminal:
            continue

        # choose one terminal oxygen to become P=O (prefer neutral O)
        terminal_sorted = sorted(
            terminal,
            key=lambda o: (o.GetFormalCharge() != 0, o.GetIdx())
        )
        o_dbl = terminal_sorted[0]

        bond = mol.GetBondBetweenAtoms(p.GetIdx(), o_dbl.GetIdx())
        if bond is not None:
            bond.SetBondType(Chem.rdchem.BondType.DOUBLE)
            o_dbl.SetFormalCharge(0)

        # remaining terminal oxygens: single bond and negative charge if not already
        for o in terminal_sorted[1:]:
            b = mol.GetBondBetweenAtoms(p.GetIdx(), o.GetIdx())
            if b is not None:
                b.SetBondType(Chem.rdchem.BondType.SINGLE)
            if o.GetFormalCharge() == 0:
                o.SetFormalCharge(-1)

    mol = mol.GetMol()
    Chem.SanitizeMol(mol)
    return mol

def qm_chg_run(pdb):
    '''
    Gaussian Workflow to guess/scan formal charges
    By default using the slurm submit wrapper in the root dir
    '''
    atom_num = 0
    pre, mul, eners, jobid, coord = pdb[0:-3], 1, {}, [], []
    for line in open(f"{pdb}", 'r'):
        if line[0:6] in ['ATOM  ', 'HETATM']:
            elem = ''.join(re.split(r'[^A-Za-z]', line.split()[-1]))
            xyzs = line[27:55]
            H = [elem] + [float(x) for x in xyzs.split()]
            coord.append(" {:5s} {:10.5f} {:10.5f} {:10.5f}\n".format(*H))
            atom_num += 1
    i = 0
    if atom_num < 60:
        cpu_num = "12"
        mem_num = "24GB"
    elif atom_num < 100:
        cpu_num = "24"
        mem_num = "48GB"
    else:
        cpu_num = "48"
        mem_num = "96GB"

    for chg in range(-8, 9):
        ### writing gau input starts
        gau_inp = f"{pre.strip('.')}_chg{chg}"
        with open(f"{gau_inp}.com", 'w') as fp:
            fp.write(f'''%chk={gau_inp}.chk
%mem={mem_num}
%nproc={cpu_num}
#T hf/6-31g(d,p)

 HF is better than b3lyp for systems like lysine/arginine

{chg} {mul}
{''.join(coord)}


''')  ### Writing gau input ends
        ter = ' termination '  # Normal/Error termination
        if os.path.isfile(f"{gau_inp}.log") and ter in subprocess.run(
                ['grep', ter, f"{gau_inp}.log"], capture_output=True, text=True).stdout:
            i += 1
            jid = 123 + i
        else:
            cmd = ['/public/home/lijunhao/soft/submit_gau.sh', '-i', f"{gau_inp}.com"]
            gaurun = subprocess.run(cmd, capture_output=True, text=True)
            jid = gaurun.stdout.strip().split()[-1]
        jobid.append(str(jid))
        eners[str(jid)] = [chg, f"{gau_inp}.log", 1.000]  # charge, fileName, initial_high_energy
    print(f"jobid list: {jobid}", flush=True)
    print(f"initail eners dict: {eners}", flush=True)
    counting = 0
    normal_ended_logs = 0
    while len(jobid) > 0:
        time.sleep(10)
        counting += 1
        # if counting >= sleep_t and normal_ended_logs > 0:
        #    # maximum of xxx minutes for long pending jobs
        #    # still risky, what if the first normal terminated log did not represent the actual charge!!!
        #    print(f'at least {sleep_t*10} seconds passed with 1 normal terminated log file for {pdb},'
        #          f'\n moving to next pdb file!',flush=True)
        #    break
        for jid in jobid:
            # grep_cmd = ["/usr/bin/grep"," termination ",eners[jid][1]]
            # ter_line = subprocess.run(grep_cmd,capture_output=True,text=True)
            time.sleep(1)
            log_lines = []
            # sacct_run = f"sacct -X --parsable2 --noheader --format State -j {jid}"
            squ_run = f"squeue -j {jid}"
            sacct_out = subprocess.run(squ_run.split(), capture_output=True, text=True)
            # if sacct_out.stdout.strip() not in ['RUNNING','PENDING']:
            if sacct_out.returncode != 0:
                time.sleep(1)
                jobid.remove(jid)
                if os.path.isfile(eners[jid][1]):
                    with open(eners[jid][1], "r") as logF:
                        log_lines = logF.readlines()
                else:
                    print(f"jobid {jid} is not running/pending, but its {eners[jid][1]} does not exist?", flush=True)
            else:
                continue

            if 'Normal termination' in str(log_lines):
                scf_lines = [x for x in log_lines if "SCF Done" in x]
                normal_ended_logs += 1
                if len(scf_lines) == 0:
                    print(f"Normal termination detected, why no SCF done for {jid}?", flush=True)
                ener = float(scf_lines[-1].split()[4])
                eners[jid][2] = ener  # eners[jobid]= [chg,logF,ener]
                # print(f"checking eners for removed jid = {jid},log = {eners[jid][1]}: {eners} ",flush=True)
            elif "Error termination" in str(log_lines):
                print(f"removing jobid {jid} due to inappropriate charge", flush=True)
                # print(f"checking eners for removed jid = {jid}, log = {eners[jid][1]}: {eners} ", flush=True)

    print(f"eners dict after all guessing: {eners},flush=True")
    min_list = min(eners.items(), key=lambda x: x[1][2])
    return min_list

def add_het_ter(pdb):
    '''
    In 5NES system, ter is need to distinguish the two ligand resn
    however, this function should be deprecated for further improvement
    '''
    with open(pdb,'r') as fp:
        lines = fp.readlines()
    output = []
    prekey = None
    terkey = None
    hetatm = False
    for i, line in enumerate(lines):
        if line.startswith('HETATM'):
            hetatm = True
            chainN = line[21:22]
            resNam = line[17:20].strip()
            resNum = line[22:26].strip()
            curkey = (chainN,resNum,resNam)
            if prekey and curkey != prekey and terkey != "TER":
                output.append("TER\n")
            prekey = curkey
            if line.startswith("TER"):
                terkey = "TER"
            else:
                terkey = None
        if line.startswith("ATOM  ") and hetatm:
            hetatm = False
            output.append("TER\n")
        output.append(line)
    return output

def mda_pre_process(pdb):
    '''
    MDAnalysis's guess_bonds module is used here to get right connection of PDB file, 
    openbabel can be used, the rules for connection records for pymol to
    correctly display the bond orders of complex protein-peptide systems is still not easy
    '''
    VDWRADII: dict[str, float] = {
        symbol.capitalize(): radius for symbol, radius in vdwradii.items()
    }
    u = mda.Universe(pdb)
    #u.atoms.guess_bonds(vdwradii=VDWRADII, fudge_factor=0.6)
    u.atoms.guess_bonds(vdwradii=VDWRADII) # type: ignore
    mol = u.atoms.convert_to("RDKIT") # type: ignore
    Chem.MolToPDBFile(mol=mol,filename="full_st_connected.pdb",flavor=4)
    return "full_st_connected.pdb"

def run_capping(pdbF, just_capping, guess_method,capping_method,charge_file, bonding_method='mda'):
    # TODO: BCC charge protocol
    pdbf = f".renamed.{pdbF}"
    nonProt, allNPs = split_nonprotein(pdbF, "./", pdbf)
    pdbf2 = 'full_st_connected.pdb'
    if bonding_method == 'mda':
        full_st_conected = mda_pre_process(pdbf)
    # the prepare_caps is more likely just working for systems with NCAAs
    _nonProt, _allNPs = find_SC_connections(pdbf2)

    nonProt += _nonProt
    allNPs += _allNPs

    blocks_list = get_blocks("full_st", pdbf2, capping_method, nonProt, allNPs)
    # if ct.lower() == 'resp':  # no need to get block charges when using BCC
    # still using blocks, non-stardard capping only applies inside the block's non-NME/ACE residues!
    if not just_capping:
        return guessing_charge(blocks_list, guess_method, charge_file)
    else:
        return blocks_list

def judge_BB(name,resn):
    # small wrapper to judge a backbone, nucleotide system is considered
    BB_atoms = {"N", "CA", "C", "O", "P", "O1P", "O2P", "O5'",
                "C5'", "C4'", "O4'", "C3'", "O3'", "C2'", "C1'"}
    if name in BB_atoms and resn in known_residues:
        return True
    return False

def find_SC_connections(pdbF,distance_cutoff=1.9):
    '''
    Using BioPython for quiker distance calculations
    '''
    from io import StringIO
    from Bio.PDB import PDBParser # type: ignore
    import itertools
    pdb_text = open(pdbF).read()
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("pdb_structure", StringIO(pdb_text))
    residues = [res for res in structure.get_residues() if res.id[0] == " "] # type: ignore
    SC_SC, SC_BB, newNPs, nonProt = [], [], [], []
    for res1, res2 in itertools.combinations(residues, 2):
        for atom1 in res1.get_atoms():
            if atom1.element == "H": continue
            for atom2 in res2.get_atoms():
                if atom2.element == "H": continue
                dist = atom1 - atom2
                if dist <= distance_cutoff:
                    a1_is_BB = judge_BB(atom1.get_name(),res1.resname)
                    a2_is_BB = judge_BB(atom2.get_name(),res2.resname)
                    if not a1_is_BB and not a2_is_BB:
                        newNPs += [[x.get_resname(),x.get_parent().id,x.id[1]] for x in [res1,res2] # type: ignore
                                   if [x.get_resname(),x.get_parent().id,x.id[1]] not in newNPs] # type: ignore
                        nonProt += [x.get_resname() for x in [res1,res2] if x.get_resname() not in nonProt]
                        SC_SC.append((res1, res2))
                    elif (not a1_is_BB and a2_is_BB) or (a1_is_BB and not a2_is_BB):
                        newNPs += [[x.get_resname(),x.get_parent().id,x.id[1]] for x in [res1,res2] # type: ignore
                                   if [x.get_resname(),x.get_parent().id,x.id[1]] not in newNPs] # type: ignore
                        SC_BB.append((res1, res2))
                        nonProt += [x.get_resname() for x in [res1, res2] if x.get_resname() not in nonProt]
    print("=== SC-SC connections ===")
    for r1, r2 in set(SC_SC):
        print(f"{r1.get_resname()} {r1.id[1]} - {r2.get_resname()} {r2.id[1]}")
    print("\n=== BB-BB connections ===")
    for r1, r2 in set(SC_BB):
        print(f"{r1.get_resname()} {r1.id[1]} - {r2.get_resname()} {r2.id[1]}")
    return nonProt, newNPs

def parse_args(argv=None):
    ''' CLI is here: more options can be added'''
    parser = argparse.ArgumentParser(
        prog="prepare_pep_caps",
        description="Prepare peptide caps for the splitted blocks\n"
                    " from a PDB file. With addional function of formal charge detection"
    )
    # Positional: pdb file
    parser.add_argument(
        "pdb_file",
        help="Input PDB file (must end with .pdb)."
    )
    parser.add_argument(
        "--just_capping",
        action="store_true", 
        default=False,
        help="Optional, Use 'just_capping' to disable charge guessing. When this is the case,\n"
             "please make sure --charge_file specified file exists and in correct format"
    )
    parser.add_argument(
        "--charge_file",
        type=str,
        default="blocks_charges.json",
        help="Charge file path, default is blocks_charges.json"
    )
    parser.add_argument(
        "--guess-method",
        choices=["mm", "qm"],
        default="mm",
        help="Method used to guess formal charges of capped blocks (mm or qm). Default: mm."
    )
    parser.add_argument(
        "--charge-method",
        choices=["resp", "bcc"],
        default="resp",
        help="Method used to calculate partial charges. Default: resp"
    )
    parser.add_argument(
        "--capping-method",
        choices=["standard", "nonStandard", "non-standard"],
        default="standard",
        help="Capping method. Default: standard."
    )
    args = parser.parse_args(argv)
    if not args.pdb_file.endswith(".pdb"):
        parser.error("pdb_file must end with .pdb")
    if args.just_capping and not os.path.isfile(args.charge_file):
        parser.error("--charge_file ('blocks_charges.json') is required when --just_capping is set")
    if args.charge_method == 'bcc' and args.capping_method == 'standard':
        parser.error("it is difficult to use standard caps for BCC charge, not recommended!")
    return args

if __name__ == "__main__":
    args = parse_args()
    pdbF = args.pdb_file
    bonding_method = 'mda' # only mda is supported for now, other methods can be added in the future
    guess_method = args.guess_method
    capping_method = args.capping_method
    #charge_method = args.charge_method
    charge_file = args.charge_file  # now available if you need it
    resp_dict = run_capping(pdbF, args.just_capping, guess_method,
                            capping_method,charge_file,bonding_method)
    print(f"resp_dict in prepare_pep_cap's main is {resp_dict}", flush=True)
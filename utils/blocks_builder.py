import re
import string
from collections import defaultdict, deque

ionsRes = ["Cl","Na","K","CL","NA","Na+","Cl-","K+","CA","MG"]
watRes = ["HOH","WAT"]
nucRes = ["DA","DG","DC","DT","A","G","A3","A5","AN","C","C3","U",
          "U3","U5","G3","G5","GN","UN"]
for_rename = list(string.ascii_uppercase)+[str(x) for x in range(0,10)]
three2one = {
    'ALA':'A', 'CYS':'C', 'ASP':'D', 'GLU':'E', 'PHE':'F',
    'GLY':'G', 'HIS':'H', 'ILE':'I', 'LYS':'K', 'LEU':'L',
    'MET':'M', 'ASN':'N', 'PRO':'P', 'GLN':'Q', 'ARG':'R',
    'SER':'S', 'THR':'T', 'VAL':'V', 'TRP':'W', 'TYR':'Y',
    'HID':'H', 'HIE':'H', 'HIP':'H', 'ASH':'D', 'GLH':'E',
    'CYX':'C',
    # capps
    'ACE':'X', 'NME':'X', 'NHE':'X',
    # phosphorylated residues recognized by Amber phosaa14SB
    'S1P':'S', 'SEP':'S', 'T1P':'T', 'TPO':'T', 'Y1P':'Y',
    'PTR':'Y', 'H1D':'H', 'H2D':'H', 'H1E':'H', 'H2E':'H',
}
known_residues = list(three2one.keys())+ionsRes+watRes+nucRes

class unitClassifier:
    """
    Classify quartets / triples / pairs as CONNECTED COMPONENTS in the eligible (non-backbone) graph,
    then post-process singles by residue name, and finally rename residues inside units to avoid
    duplicated residue names (and to normalize standard residues from three2one).

    Eligible "non-backbone connection":
      - NOT a peptide backbone C–N link, AND
      - at least one atom name is not in MAINCHAIN {"N","C"}.
    """
    #TODO: the data structure should be unified, needs to deprecate "a","b" keys in pairs
    SEL_RE = re.compile(
    r"""r\.\s*(?P<resn>\S+)\s*&\s*c\.\s*(?P<chain>\S*)\s*&\s*i\.\s*(?P<resi>\S+)\s*&\s*name\s+"?(?P<atom>[^"]+)"?""",
    re.IGNORECASE
    )
    MAINCHAIN = {"N", "C"}

    @staticmethod
    def parse_sel(sel: str):
        m = unitClassifier.SEL_RE.fullmatch(sel.strip())
        if not m:
            raise ValueError(f"Unrecognized selection: {sel!r}")
        return (m.group("resn"), m.group("chain"), m.group("resi"), m.group("atom"))

    @staticmethod
    def is_backbone_link(a_atom: str, b_atom: str) -> bool:
        # peptide bond backbone strictly when C–N (either direction)
        return {a_atom.upper(), b_atom.upper()} == {"C", "N"}

    @staticmethod
    def _norm(resn: str) -> str:
        return (resn or "").upper()

    @staticmethod
    def _resi_key(resi: str):
        # stable numeric sort when possible
        try:
            return (0, int(resi))
        except Exception:
            return (1, resi)

    @staticmethod
    def _canon_edge(u, u_atom, v, v_atom):
        # canonicalize atom-edge so reverse duplicates collapse
        if u <= v:
            return (u, u_atom, v, v_atom)
        return (v, v_atom, u, u_atom)

    def __init__(self):
        self.all_nodes = set()              # set of (resn, chain, resi)
        self.atom_edges = []                # list of dicts: u,u_atom,v,v_atom,eligible
        self.adj = defaultdict(set)         # eligible adjacency (undirected)
        self.rename_map = {}                # populated by classify()

    # -----------------------
    # Build eligible graph
    # -----------------------
    def build(self, linking_dicts):
        for d in linking_dicts:
            for k, v in d.items():
                a = self.parse_sel(k)
                b = self.parse_sel(v)

                u, u_atom = (a[0], a[1], a[2]), a[3]
                w, w_atom = (b[0], b[1], b[2]), b[3]
                self.all_nodes.update([u, w])

                backbone = self.is_backbone_link(u_atom, w_atom)
                eligible = (not backbone) and (
                    (u_atom.upper() not in self.MAINCHAIN)
                    or (w_atom.upper() not in self.MAINCHAIN)
                )

                self.atom_edges.append(
                    {
                        "u": u,
                        "u_atom": u_atom,
                        "v": w,
                        "v_atom": w_atom,
                        "eligible": eligible,
                    }
                )

                if eligible:
                    self.adj[u].add(w)
                    self.adj[w].add(u)

    def connected_components(self):
        seen = set()
        comps = []

        for n in sorted(self.all_nodes):
            if n in seen:
                continue

            # isolated in eligible graph
            if n not in self.adj or not self.adj[n]:
                seen.add(n)
                comps.append({n})
                continue

            q = deque([n])
            seen.add(n)
            comp = {n}

            while q:
                x = q.popleft()
                for y in self.adj.get(x, ()):
                    if y in seen:
                        continue
                    seen.add(y)
                    comp.add(y)
                    q.append(y)

            comps.append(comp)

        return comps

    def edges_within(self, nodes_set):
        """Deduplicated eligible atom-edges fully inside nodes_set."""
        out = []
        seen = set()

        for e in self.atom_edges:
            if not e["eligible"]:
                continue
            if e["u"] in nodes_set and e["v"] in nodes_set:
                key = self._canon_edge(e["u"], e["u_atom"], e["v"], e["v_atom"])
                if key in seen:
                    continue
                seen.add(key)
                out.append(((key[0], key[1]), (key[2], key[3])))

        out.sort(key=lambda x: (x[0][0], x[0][1], x[1][0], x[1][1]))
        return out

    # -----------------------
    # Globals helpers
    # -----------------------
    def _get_three2one_keys(self):
        """
        three2one is expected to be a GLOBAL dict.
        We only need its .keys() for filtering/standard-AA detection.
        """
        t = globals().get("three2one", None)
        if t is None:
            return set()
        try:
            return {self._norm(k) for k in t.keys()}
        except Exception:
            return set()

    def _get_for_rename(self):
        """
        for_rename is expected to be a GLOBAL list like:
          list(string.ascii_uppercase) + [str(x) for x in range(10)]
        If missing, use that default.
        """
        fr = globals().get("for_rename", None)
        if fr is None:
            fr = list(string.ascii_uppercase) + [str(x) for x in range(10)]
        return list(fr)

    def _postprocess_singles_by_resname(self, raw_singles, quartets, triples, pairs):
        """
        Bases on residue name only:
          - de-duplicate residue names in singles
          - remove residue names that occur in three2one.keys()
          - remove residue names that occur in pairs/triples/quartets
        """
        unit_resnames = set()

        for q in quartets:
            for node in q.get("nodes", []):
                unit_resnames.add(self._norm(node[0]))

        for t in triples:
            for node in t.get("nodes", []):
                unit_resnames.add(self._norm(node[0]))

        for p in pairs:
            unit_resnames.add(self._norm(p["a"][0]))
            unit_resnames.add(self._norm(p["b"][0]))

        banned = unit_resnames | self._get_three2one_keys()

        seen = set()
        singles = []
        for n in raw_singles:
            r = self._norm(n[0])
            if r in banned:
                continue
            if r in seen:
                continue
            seen.add(r)
            singles.append(n)

        return singles

    # -----------------------
    def _rename_units_by_resname(self, quartets, triples, pairs):
        """
        Renaming rules:

        Let base = RESN uppercased. Let for_rename = [A..Z,0..9].
        Let start_idx be the index of base[-1] in for_rename (or 0 if missing).

        1) If base in three2one.keys() (standard AA):
           - ALWAYS rename (even if only one occurrence):
             occurrence i => base[:-1] + lower(for_rename[start_idx + i])

           Examples:
             GLY -> GLy
             HIP -> HIp
             THR -> THr
             CYS occurrences -> CYs, CYt, CYu, ...

        2) If base NOT in three2one.keys() (nonstandard / modified):
           - first occurrence keeps original spelling (node[0] as-is)
           - subsequent occurrences rename using:
             occurrence i (i>=1) => base[:-1] + lower(for_rename[start_idx + (i-1)])

           Example:
             DCY occurrences -> DCY, DCy, DCz, DC0, ...
             DBB occurrences -> DBB, DBb, DBc, DBd, ...

        Occurrence order is by (chain, resi numeric).
        Renaming is applied consistently across *all* nodes in pairs/triples/quartets,
        and edges are updated accordingly.

        Returns (new_quartets, new_triples, new_pairs, rename_map)
        rename_map format:
          { (old_resn, chain, resi): (new_resn, chain, resi), ... }
        """
        three2one_keys = self._get_three2one_keys()
        fr = self._get_for_rename()
        rename_map = {}
        if not fr:
            return quartets, triples, pairs, rename_map

        # collect all unit nodes (unique by node tuple)
        nodes = []
        for q in quartets:
            nodes.extend(q.get("nodes", []))
        for t in triples:
            nodes.extend(t.get("nodes", []))
        for p in pairs:
            nodes.append(p["a"])
            nodes.append(p["b"])

        uniq_nodes = sorted(
            set(nodes),
            key=lambda n: (self._norm(n[0]), n[1], self._resi_key(n[2]), n[2]),
        )

        by_resn = defaultdict(list)
        for n in uniq_nodes:
            by_resn[self._norm(n[0])].append(n)

        # build node->node mapping, then expose as (resn,chain,resi)->(...)
        node_map = {}

        for base, group in by_resn.items():
            group_sorted = sorted(group, key=lambda n: (n[1], self._resi_key(n[2]), n[2]))
            is_std = base in three2one_keys

            last = base[-1] if base else ""
            try:
                start_idx = fr.index(last)
            except ValueError:
                start_idx = 0

            for i, node in enumerate(group_sorted):
                old_resn, chain, resi = node

                # keep original if degenerate base
                if len(base) < 2:
                    new_resn = old_resn
                else:
                    if is_std:
                        # standard: ALWAYS rename, even first occurrence
                        offset = i
                    else:
                        # nonstandard: keep first, rename 2nd+
                        if i == 0:
                            offset = None
                        else:
                            offset = i - 1

                    if offset is None:
                        new_resn = old_resn
                    else:
                        sym = fr[(start_idx + offset) % len(fr)]
                        sym = sym.lower() if str(sym).isalpha() else str(sym)
                        new_resn = base[:-1] + sym

                new_node = (new_resn, chain, resi)
                node_map[node] = new_node
                rename_map[node] = new_node  # same shape (resn,chain,resi)->(resn,chain,resi)

        def map_node(n):
            return node_map.get(n, n)

        def map_edge(edge):
            (n1, a1), (n2, a2) = edge
            return ((map_node(n1), a1), (map_node(n2), a2))

        new_quartets = []
        for q in quartets:
            new_quartets.append(
                {
                    "nodes": [map_node(n) for n in q.get("nodes", [])],
                    "edges": [map_edge(e) for e in q.get("edges", [])],
                }
            )

        new_triples = []
        for t in triples:
            new_triples.append(
                {
                    "nodes": [map_node(n) for n in t.get("nodes", [])],
                    "edges": [map_edge(e) for e in t.get("edges", [])],
                }
            )

        new_pairs = []
        for p in pairs:
            new_pairs.append(
                {
                    "a": map_node(p["a"]),
                    "b": map_node(p["b"]),
                    "edges": [map_edge(e) for e in p.get("edges", [])],
                }
            )

        return new_quartets, new_triples, new_pairs, rename_map

    # -----------------------
    # Main API
    # -----------------------
    def classify(self, linking_dicts):
        self.__init__()
        self.build(linking_dicts)

        comps = self.connected_components()

        quartets, triples, pairs = [], [], []
        used_nodes = set()

        for comp in comps:
            sz = len(comp)
            # for one linker with 4 side-chain connections with other residues
            # it will be better to use non-block capping methods!
            if sz == 4:
                quartets.append({"nodes": sorted(comp), "edges": self.edges_within(comp)})
                used_nodes |= comp
            elif sz == 3:
                triples.append({"nodes": sorted(comp), "edges": self.edges_within(comp)})
                used_nodes |= comp
            elif sz == 2:
                u, v = sorted(comp)
                pairs.append({"a": u, "b": v, "edges": self.edges_within(comp)})
                used_nodes |= comp

        # singles: nodes not used by any 2/3/4-mer unit, then postprocess by residue name
        raw_singles = [n for n in sorted(self.all_nodes) if n not in used_nodes]
        singles = self._postprocess_singles_by_resname(raw_singles, quartets, triples, pairs)

        # rename residues inside units + get mapping
        quartets, triples, pairs, rename_map = self._rename_units_by_resname(quartets, triples, pairs)
        self.rename_map = rename_map  # also store on the instance

        units = {"quartets": quartets, "triples": triples, "pairs": pairs, "singles": singles}
        units["rename_map"] = rename_map  # convenient for external PyMOL alter()
        return units

"""
Stage 3 - is the pocket actually druggable? Size, burial, hydrophobicity.

Three measurements on the model (all geometric, all reproducible):

1. Cavity volume. A 0.7 A grid is laid over a 22 A cube at the pocket
   centre. A point counts as cavity if it is 3.2-6.5 A from the nearest
   protein heavy atom (outside van der Waals contact but still in the first
   shell) and buried, where buried means >= 65% of 42 evenly spread rays
   from that point hit protein within 12 A. The connected component
   containing the co-crystal ligand is the pocket; its volume is the point
   count times 0.343 A^3.

2. Burial of the reference ligand: Shrake-Rupley SASA of AR-A014418 alone
   minus its SASA in the complex, as a percentage.

3. Hydrophobicity of the lining: fraction of lining heavy atoms that are
   apolar carbon, and the residue-type split.

Output: 03_pocket/druggability.json / .md
"""
import json, os
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, ".."))
MODEL = os.path.join(ROOT, "02_structure", "models", "gsk3b_model_aligned.pdb")
REF = os.path.join(ROOT, "05_docking", "ref_ligand_TMU.pdb")
POCK = os.path.join(ROOT, "03_pocket")

HYDROPHOBIC_RES = {"ALA", "VAL", "LEU", "ILE", "MET", "PHE", "TRP", "PRO", "GLY", "CYS"}
GRID = 0.7
CUBE = 22.0
D_MIN, D_MAX = 3.2, 6.5
BURIAL_CUT = 0.65
N_RAYS = 42


def fibonacci_sphere(n):
    i = np.arange(n) + 0.5
    phi = np.arccos(1 - 2 * i / n)
    theta = np.pi * (1 + 5 ** 0.5) * i
    return np.stack([np.cos(theta) * np.sin(phi),
                     np.sin(theta) * np.sin(phi), np.cos(phi)], axis=1)


def main():
    from Bio.PDB import PDBParser, ShrakeRupley
    p = PDBParser(QUIET=True)
    st = p.get_structure("m", MODEL)
    prot = np.array([a.coord for a in st.get_atoms() if a.element != "H"])
    lig = np.array([[float(l[30:38]), float(l[38:46]), float(l[46:54])]
                    for l in open(REF) if l.startswith("HETATM")])
    centre = lig.mean(axis=0)

    # ---- 1. cavity volume ------------------------------------------------
    n = int(CUBE / GRID)
    ax = (np.arange(n) - n / 2) * GRID
    gx, gy, gz = np.meshgrid(ax, ax, ax, indexing="ij")
    pts = np.stack([gx.ravel(), gy.ravel(), gz.ravel()], axis=1) + centre

    from scipy.spatial import cKDTree
    tree = cKDTree(prot)
    dmin, _ = tree.query(pts, k=1)
    keep = (dmin > D_MIN) & (dmin < D_MAX)
    cand = pts[keep]

    rays = fibonacci_sphere(N_RAYS)
    steps = np.arange(1.5, 12.0, 0.7)
    buried = np.zeros(len(cand), dtype=float)
    for r in rays:
        probe = cand[:, None, :] + r[None, None, :] * steps[None, :, None]
        d, _ = tree.query(probe.reshape(-1, 3), k=1,
                          distance_upper_bound=2.6)
        hit = np.isfinite(d).reshape(len(cand), len(steps)).any(axis=1)
        buried += hit
    buried /= len(rays)
    cav = cand[buried >= BURIAL_CUT]

    # connected component containing the ligand
    ct = cKDTree(cav)
    pairs = ct.query_pairs(GRID * 1.45, output_type="ndarray")
    import scipy.sparse as sp
    from scipy.sparse.csgraph import connected_components
    g = sp.coo_matrix((np.ones(len(pairs)), (pairs[:, 0], pairs[:, 1])),
                      shape=(len(cav), len(cav)))
    ncomp, lab = connected_components(g, directed=False)
    seed = ct.query(centre, k=1)[1]
    pocket_pts = cav[lab == lab[seed]]
    volume = len(pocket_pts) * GRID ** 3

    # ---- 2. burial of the reference ligand -------------------------------
    sr = ShrakeRupley()
    from Bio.PDB.Structure import Structure
    from Bio.PDB.Model import Model
    from Bio.PDB.Chain import Chain
    from Bio.PDB.Residue import Residue
    from Bio.PDB.Atom import Atom

    def mk_lig():
        s = Structure("l"); m = Model(0); c = Chain("Z")
        res = Residue((" ", 1, " "), "LIG", "")
        for i, xyz in enumerate(lig):
            a = Atom(f"C{i}", xyz, 0.0, 1.0, " ", f"C{i}", i, element="C")
            res.add(a)
        c.add(res); m.add(c); s.add(m)
        return s

    free = mk_lig(); sr.compute(free, level="R")
    sasa_free = list(free.get_residues())[0].sasa

    both = p.get_structure("m2", MODEL)
    ch = list(both[0])[0]
    res = Residue((" ", 9001, " "), "LIG", "")
    for i, xyz in enumerate(lig):
        res.add(Atom(f"C{i}", xyz, 0.0, 1.0, " ", f"C{i}", 10000 + i, element="C"))
    ch.add(res)
    sr.compute(both[0], level="R")
    sasa_bound = res.sasa
    burial = 100.0 * (1 - sasa_bound / sasa_free)

    # ---- 3. hydrophobicity of the lining ---------------------------------
    lining, lin_atoms = {}, []
    for a in st.get_atoms():
        if a.element == "H":
            continue
        if np.linalg.norm(lig - a.coord, axis=1).min() < 5.0:
            r = a.get_parent()
            lining[(r.id[1] + 34, r.get_resname())] = True
            lin_atoms.append(a.element)
    carbon_frac = 100.0 * sum(1 for e in lin_atoms if e == "C") / len(lin_atoms)
    hydro_res = sum(1 for (_, rn) in lining if rn in HYDROPHOBIC_RES)

    out = {"cavity_volume_A3": round(volume, 1),
           "cavity_grid_points": int(len(pocket_pts)),
           "grid_spacing_A": GRID, "burial_cutoff": BURIAL_CUT,
           "reference_ligand_buried_pct": round(burial, 1),
           "lining_residues": len(lining),
           "lining_hydrophobic_residues": hydro_res,
           "lining_hydrophobic_residue_pct": round(100 * hydro_res / len(lining), 1),
           "lining_carbon_atom_pct": round(carbon_frac, 1),
           "lining_residue_list": sorted(f"{rn.title()}{num}" for num, rn in lining)}
    json.dump(out, open(os.path.join(POCK, "druggability.json"), "w"), indent=1)
    with open(os.path.join(POCK, "druggability.md"), "w", encoding="utf8") as fh:
        fh.write("# ATP-pocket druggability (geometric)\n\n")
        for k, v in out.items():
            fh.write(f"- **{k}**: {v}\n")
    for k, v in out.items():
        print(f"  {k}: {v}")


if __name__ == "__main__":
    main()

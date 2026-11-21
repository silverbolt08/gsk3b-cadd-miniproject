"""
Stage 3b - define the binding site and write the Vina search box.

Works on the aligned model (crystal frame, see 05a_align_model.py), so the
same box applies to the model and to the 1Q5K crystal receptor.

The box is anchored on two independent pieces of evidence:
  (a) P2Rank pocket 1 on the model - the top-ranked pocket, whose residues
      are the canonical ATP site;
  (b) the crystallographic position of AR-A014418 (ligand TMU) in 1Q5K.
The centre is the midpoint of the two centroids. Each edge is the extent of
the pocket-1 residues together with the reference ligand plus a 5 A margin,
rounded up and clipped to 22-26 A: large enough for every library member to
rotate freely, small enough to keep the search inside the cleft.

Outputs: 03_pocket/box.json, 03_pocket/pocket_residues.csv,
         05_docking/ref_ligand_TMU.pdb
"""
import os, json, csv
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
MOD = os.path.join(HERE, "..", "02_structure", "models", "gsk3b_model_aligned.pdb")
TPL = os.path.join(HERE, "..", "02_structure", "templates")
POCK = os.path.join(HERE, "..", "03_pocket")
DOCK = os.path.join(HERE, "..", "05_docking")
OFFSET = 34
MARGIN = 5.0


def reference_ligand():
    """TMU = AR-A014418 as it sits in 1Q5K chain A (crystal frame)."""
    from Bio.PDB import PDBParser
    c = PDBParser(QUIET=True).get_structure("c", os.path.join(TPL, "1Q5K.pdb"))[0]
    res = [r for r in c["A"] if r.get_resname() == "TMU"][0]
    atoms = list(res)
    os.makedirs(DOCK, exist_ok=True)
    out = os.path.join(DOCK, "ref_ligand_TMU.pdb")
    from Bio.PDB import PDBIO, Select

    class OnlyTMU(Select):
        def accept_residue(self, r):
            return r.get_resname() == "TMU" and r.get_parent().id == "A"

    io = PDBIO()
    io.set_structure(c.get_parent())
    io.save(out, OnlyTMU())
    print(f"  reference ligand TMU ({len(atoms)} heavy atoms) -> {out}")
    return np.array([a.coord for a in atoms])


def p2rank_pocket1():
    f = os.path.join(POCK, "p2rank", "gsk3b_model_aligned.pdb_predictions.csv")
    rows = list(csv.reader(open(f)))
    hdr = [h.strip() for h in rows[0]]
    r = dict(zip(hdr, [c.strip() for c in rows[1]]))
    centre = np.array([float(r["center_x"]), float(r["center_y"]), float(r["center_z"])])
    nums = sorted(int(x.split("_")[1]) for x in r["residue_ids"].split() if x)
    return centre, nums, float(r["score"]), float(r["probability"])


def pocket_atoms(resnums_model):
    from Bio.PDB import PDBParser
    ch = list(PDBParser(QUIET=True).get_structure("m", MOD)[0])[0]
    pts = []
    for r in ch:
        if r.id[0] == " " and r.id[1] in resnums_model:
            pts.extend([a.coord for a in r if a.element != "H"])
    return np.array(pts)


def main():
    lig_xyz = reference_ligand()
    lig_centre = lig_xyz.mean(axis=0)
    p_centre, p_resnums, p_score, p_prob = p2rank_pocket1()
    print(f"  P2Rank pocket 1: score {p_score}, probability {p_prob}, "
          f"{len(p_resnums)} residues")
    print(f"  P2Rank centre   {p_centre.round(2)}")
    print(f"  TMU centroid    {lig_centre.round(2)}  "
          f"(separation {np.linalg.norm(p_centre - lig_centre):.2f} A)")

    centre = (p_centre + lig_centre) / 2.0
    pts = pocket_atoms(p_resnums)
    lo = np.minimum(pts.min(axis=0), lig_xyz.min(axis=0))
    hi = np.maximum(pts.max(axis=0), lig_xyz.max(axis=0))
    size = np.clip(np.ceil((hi - lo) + 2 * MARGIN), 22.0, 26.0)

    uni = [n + OFFSET for n in p_resnums]
    box = {"frame": "1Q5K crystal frame",
           "center": [round(float(x), 3) for x in centre],
           "size": [float(x) for x in size],
           "p2rank_center": [round(float(x), 3) for x in p_centre],
           "tmu_centroid": [round(float(x), 3) for x in lig_centre],
           "p2rank_pocket1_score": p_score, "p2rank_pocket1_probability": p_prob,
           "pocket_residues_model": p_resnums, "pocket_residues_uniprot": uni,
           "margin_A": MARGIN, "offset_model_to_uniprot": OFFSET}
    json.dump(box, open(os.path.join(POCK, "box.json"), "w"), indent=1)

    from Bio.PDB import PDBParser
    ch = list(PDBParser(QUIET=True).get_structure("m", MOD)[0])[0]
    with open(os.path.join(POCK, "pocket_residues.csv"), "w", newline="") as fh:
        w = csv.writer(fh); w.writerow(["uniprot_resnum", "model_resnum", "resname"])
        for n in p_resnums:
            rn = [r.get_resname() for r in ch if r.id[0] == " " and r.id[1] == n]
            w.writerow([n + OFFSET, n, rn[0] if rn else "?"])

    print(f"\n  box centre {box['center']}  size {box['size']} A")
    print(f"  pocket residues (UniProt): {uni}")
    print("  wrote", os.path.join(POCK, "box.json"))


if __name__ == "__main__":
    main()

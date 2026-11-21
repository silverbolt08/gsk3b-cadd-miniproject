"""
Protocol validation - redock the co-crystallised inhibitor.

AR-A014418 (HET code TMU) is resolved in 1Q5K chain A. We dock the same
molecule, prepared from SMILES by our own pipeline, back into
  (a) the 1Q5K crystal receptor, and
  (b) the MODELLER model,
and measure the symmetry-corrected heavy-atom RMSD of each predicted pose
against the crystallographic pose. Under 2.0 A is the usual pass mark; it
tells us the box, the receptor preparation and the scoring function can
reproduce a known answer before we trust them on unknowns.

Run after: python 06_dock.py --receptor crystal --only ACT10
           python 06_dock.py --receptor model   --only ACT10

Outputs: 05_docking/redock_validation.json / .md
"""
import os, json, subprocess, sys
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, ".."))
REF_PDB = os.path.join(ROOT, "05_docking", "ref_ligand_TMU.pdb")
# AR-A014418: N-(4-methoxybenzyl)-N'-(5-nitro-1,3-thiazol-2-yl)urea
TMU_SMILES = "COc1ccc(CNC(=O)Nc2ncc([N+](=O)[O-])s2)cc1"
LIG_ID = "ACT10"


def reference_mol():
    from rdkit import Chem
    from rdkit.Chem import AllChem
    ref = Chem.MolFromPDBFile(REF_PDB, removeHs=True, sanitize=False)
    tpl = Chem.MolFromSmiles(TMU_SMILES)
    ref = AllChem.AssignBondOrdersFromTemplate(tpl, ref)
    return ref


def poses(pdbqt):
    from meeko import PDBQTMolecule, RDKitMolCreate
    pm = PDBQTMolecule.from_file(pdbqt, skip_typing=True)
    mols = RDKitMolCreate.from_pdbqt_mol(pm, only_cluster_leads=False)
    return mols[0]        # one RDKit mol holding every pose as a conformer


def rmsds(ref, docked):
    from rdkit import Chem
    from rdkit.Chem import rdMolAlign
    out = []
    for i in range(docked.GetNumConformers()):
        probe = Chem.Mol(docked, confId=i) if False else Chem.Mol(docked)
        probe.RemoveAllConformers()
        probe.AddConformer(docked.GetConformer(i), assignId=True)
        probe = Chem.RemoveHs(probe)
        out.append(rdMolAlign.CalcRMS(probe, ref))
    return out


def scores(pdbqt):
    s = []
    for line in open(pdbqt):
        if line.startswith("REMARK VINA RESULT"):
            s.append(float(line.split()[3]))
    return s


def main():
    ref = reference_mol()
    result = {"ligand": "AR-A014418 (TMU in 1Q5K)", "smiles": TMU_SMILES,
              "reference_atoms": ref.GetNumAtoms(), "receptors": {}}
    for tag in ("crystal", "model"):
        p = os.path.join(ROOT, "05_docking", tag, "poses", LIG_ID + ".pdbqt")
        if not os.path.exists(p):
            print(f"!! {p} missing - run 06_dock.py --receptor {tag} --only {LIG_ID}")
            continue
        d = poses(p)
        r = rmsds(ref, d)
        sc = scores(p)
        best_i = int(np.argmin(r))
        result["receptors"][tag] = {
            "top_pose_score": sc[0], "top_pose_rmsd": round(r[0], 2),
            "best_rmsd_any_mode": round(min(r), 2), "best_rmsd_mode": best_i + 1,
            "score_of_best_rmsd_mode": sc[best_i],
            "all_scores": sc, "all_rmsd": [round(x, 2) for x in r]}
        print(f"{tag:8s} top pose {sc[0]:6.2f} kcal/mol  RMSD {r[0]:.2f} A   "
              f"(best of {len(r)} modes: {min(r):.2f} A, mode {best_i+1})")

    json.dump(result, open(os.path.join(ROOT, "05_docking", "redock_validation.json"), "w"), indent=1)
    with open(os.path.join(ROOT, "05_docking", "redock_validation.md"), "w", encoding="utf8") as fh:
        fh.write("# Redocking validation - AR-A014418 into GSK-3beta\n\n")
        fh.write("| receptor | top-pose score (kcal/mol) | top-pose RMSD (A) | best RMSD over 9 modes (A) |\n")
        fh.write("|---|---|---|---|\n")
        for tag, v in result["receptors"].items():
            fh.write(f"| {tag} | {v['top_pose_score']:.2f} | {v['top_pose_rmsd']:.2f} | "
                     f"{v['best_rmsd_any_mode']:.2f} (mode {v['best_rmsd_mode']}) |\n")
    print("wrote", os.path.join(ROOT, "05_docking", "redock_validation.md"))


if __name__ == "__main__":
    main()

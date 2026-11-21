"""
Stage 4 - ligand preparation.

For every molecule (11 positive controls + 128 library compounds):
  1. Open Babel sets the protonation state at pH 7.4
     (OBMol.AddHydrogens(polaronly=False, correctForPH=True, pH=7.4))
  2. RDKit generates one 3D conformer with ETKDGv3 (seed 0xf00d) and
     minimises it with MMFF94s (500 steps); if MMFF typing fails we fall
     back to UFF.
  3. Meeko (mk_prepare_ligand) writes the AutoDock PDBQT: Gasteiger charges,
     rotatable bonds detected automatically, amide bonds kept rigid.

Outputs: 04_ligands/sdf/<id>.sdf , 04_ligands/pdbqt/<id>.pdbqt
         04_ligands/ligand_manifest.csv
"""
import csv, os, sys, warnings

warnings.filterwarnings("ignore")
HERE = os.path.dirname(os.path.abspath(__file__))
LIGDIR = os.path.join(HERE, "..", "04_ligands")
SDF = os.path.join(LIGDIR, "sdf")
PDBQT = os.path.join(LIGDIR, "pdbqt")
SEED = 0xF00D

from rdkit import Chem, RDLogger
from rdkit.Chem import AllChem, Descriptors, Crippen, rdMolDescriptors
RDLogger.DisableLog("rdApp.*")
from openbabel import pybel
from meeko import MoleculePreparation, PDBQTWriterLegacy


def protonate_ph74(smiles):
    """Return the SMILES of the dominant microspecies at pH 7.4 (Open Babel)."""
    m = pybel.readstring("smi", smiles)
    m.OBMol.AddHydrogens(False, True, 7.4)
    out = m.write("smi").split()[0]
    return out


def embed3d(smiles):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None, "unparseable"
    mol = Chem.AddHs(mol)
    ps = AllChem.ETKDGv3()
    ps.randomSeed = SEED
    ps.useSmallRingTorsions = True
    if AllChem.EmbedMolecule(mol, ps) != 0:
        ps.useRandomCoords = True
        if AllChem.EmbedMolecule(mol, ps) != 0:
            return None, "embed_failed"
    try:
        if AllChem.MMFFHasAllMoleculeParams(mol):
            AllChem.MMFFOptimizeMolecule(mol, mmffVariant="MMFF94s", maxIters=500)
            ff = "MMFF94s"
        else:
            AllChem.UFFOptimizeMolecule(mol, maxIters=500)
            ff = "UFF"
    except Exception:
        ff = "none"
    return mol, ff


def to_pdbqt(mol, path):
    prep = MoleculePreparation(rigid_macrocycles=False)
    setups = prep.prepare(mol)
    s, ok, err = PDBQTWriterLegacy.write_string(setups[0])
    if not ok:
        return False, err
    open(path, "w").write(s)
    return True, ""


def props(mol):
    m = Chem.RemoveHs(mol)
    return dict(mw=round(Descriptors.MolWt(m), 2),
                clogp=round(Crippen.MolLogP(m), 2),
                tpsa=round(rdMolDescriptors.CalcTPSA(m), 1),
                hbd=rdMolDescriptors.CalcNumHBD(m),
                hba=rdMolDescriptors.CalcNumHBA(m),
                rotb=rdMolDescriptors.CalcNumRotatableBonds(m),
                arom_rings=rdMolDescriptors.CalcNumAromaticRings(m),
                heavy=m.GetNumHeavyAtoms(),
                formula=rdMolDescriptors.CalcMolFormula(m))


def load_inputs():
    items = []
    ap = os.path.join(LIGDIR, "actives", "known_actives.csv")
    for i, r in enumerate(csv.DictReader(open(ap, encoding="utf8")), 1):
        items.append({"id": f"ACT{i:02d}", "name": r["name"], "smiles": r["smiles"],
                      "set": "known_active", "chembl_id": r["chembl_id"],
                      "potency_nM": r.get("potency_nM", ""),
                      "note": r.get("mechanism_note", "")})
    lp = os.path.join(LIGDIR, "library", "screening_library.csv")
    for r in csv.DictReader(open(lp, encoding="utf8")):
        items.append({"id": r["ligand_id"], "name": r["name"], "smiles": r["smiles"],
                      "set": r["block"], "chembl_id": r["chembl_id"],
                      "potency_nM": "", "note": ""})
    return items


def main():
    os.makedirs(SDF, exist_ok=True); os.makedirs(PDBQT, exist_ok=True)
    rows, fails = [], []
    for it in load_inputs():
        try:
            smi74 = protonate_ph74(it["smiles"])
        except Exception:
            smi74 = it["smiles"]
        mol, ff = embed3d(smi74)
        if mol is None:
            mol, ff = embed3d(it["smiles"])          # retry with neutral form
            smi74 = it["smiles"]
        if mol is None:
            fails.append((it["id"], it["name"], ff)); continue
        mol.SetProp("_Name", f"{it['id']}_{it['name']}")
        sdf_path = os.path.join(SDF, it["id"] + ".sdf")
        w = Chem.SDWriter(sdf_path); w.write(mol); w.close()
        pq = os.path.join(PDBQT, it["id"] + ".pdbqt")
        ok, err = to_pdbqt(mol, pq)
        if not ok:
            fails.append((it["id"], it["name"], "pdbqt:" + str(err))); continue
        p = props(mol)
        charge = Chem.GetFormalCharge(mol)
        rows.append({**{k: it[k] for k in ("id", "name", "set", "chembl_id", "potency_nM", "note")},
                     "smiles_input": it["smiles"], "smiles_ph74": smi74,
                     "formal_charge_ph74": charge, "forcefield": ff, **p})
        print(f"  {it['id']:7s} {it['name'][:28]:28s} {p['formula']:14s} q={charge:+d} {ff}")

    man = os.path.join(LIGDIR, "ligand_manifest.csv")
    with open(man, "w", newline="", encoding="utf8") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys())); w.writeheader(); w.writerows(rows)
    print(f"\nprepared {len(rows)} ligands -> {man}")
    if fails:
        print("FAILED:")
        for f in fails:
            print("  ", f)


if __name__ == "__main__":
    main()

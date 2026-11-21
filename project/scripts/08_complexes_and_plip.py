"""
Stage 6a - build the receptor-ligand complexes and profile the interactions
with PLIP.

For every compound asked for (default: the 11 known actives plus the top 20
library hits) the best-scoring Vina pose is written into a complex PDB with
the prepared receptor, then PLIP is run on that complex. The parsed
interactions become
  - a per-interaction table (06_analysis/plip/interactions.csv)
  - a residue x compound interaction fingerprint
    (06_analysis/plip/fingerprint.csv)
so a hit can be compared residue by residue with the known actives.

Usage: python 08_complexes_and_plip.py [--receptor model] [--top 20]
"""
import argparse, csv, json, os, sys
from collections import defaultdict

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, ".."))
ANA = os.path.join(ROOT, "06_analysis")
CPX = os.path.join(ANA, "complexes")
PLIPDIR = os.path.join(ANA, "plip")
OFFSET = 34   # model numbering -> UniProt; set to 0 for the crystal receptor

INTERACTION_KINDS = ["hydrophobic", "hbond", "waterbridge", "saltbridge",
                     "pistacking", "pication", "halogen"]


def pose_to_pdb_block(pdbqt, resname="LIG", chain="Z", resnum=1, start_serial=90000):
    """best (first) pose from a Vina PDBQT -> PDB HETATM lines with hydrogens"""
    from meeko import PDBQTMolecule, RDKitMolCreate
    from rdkit import Chem
    pm = PDBQTMolecule.from_file(pdbqt, skip_typing=True)
    mol = RDKitMolCreate.from_pdbqt_mol(pm, only_cluster_leads=False)[0]
    conf = mol.GetConformer(0)
    lines, n = [], start_serial
    counts = defaultdict(int)
    for a in mol.GetAtoms():
        n += 1
        el = a.GetSymbol()
        counts[el] += 1
        name = f"{el}{counts[el]}"[:4]
        p = conf.GetAtomPosition(a.GetIdx())
        lines.append(f"HETATM{n % 100000:5d} {name:<4s} {resname:<3s} {chain}{resnum:4d}    "
                     f"{p.x:8.3f}{p.y:8.3f}{p.z:8.3f}  1.00  0.00          {el:>2s}\n")
    return lines, Chem.MolToSmiles(Chem.RemoveHs(mol))


def build_complex(receptor_pdb, pose_pdbqt, out_pdb):
    rec = [l for l in open(receptor_pdb) if l.startswith(("ATOM", "HETATM"))]
    lig, smi = pose_to_pdb_block(pose_pdbqt)
    with open(out_pdb, "w") as fh:
        fh.writelines(rec)
        fh.write("TER\n")
        fh.writelines(lig)
        fh.write("END\n")
    return out_pdb, smi


def _patch_openbabel_inchi():
    """The Open Babel wheel used here is built without the InChI plugin, and
    PLIP asks every ligand for an InChI key purely for its report header.
    Return an empty string for that one format instead of letting it raise."""
    from openbabel import pybel
    if getattr(pybel.Molecule, "_inchi_patched", False):
        return
    original = pybel.Molecule.write

    def write(self, format="smi", *args, **kwargs):
        fmt = kwargs.get("format", format)
        if fmt in ("inchikey", "inchi"):
            return ""
        return original(self, format, *args, **kwargs)

    pybel.Molecule.write = write
    pybel.Molecule._inchi_patched = True


def run_plip(complex_pdb):
    _patch_openbabel_inchi()
    from plip.structure.preparation import PDBComplex
    c = PDBComplex()
    c.load_pdb(complex_pdb)
    c.analyze()
    for key, site in c.interaction_sets.items():
        if key.startswith("LIG"):
            return site
    return None


def collect(site):
    """flatten a PLIP binding site into (kind, resnum, resname, detail) rows
    plus the geometry needed to draw each interaction in a figure"""
    rows, geom = [], []

    def add(kind, resnr, restype, detail, p1=None, p2=None):
        rows.append((kind, resnr, restype, detail))
        if p1 is not None and p2 is not None:
            geom.append({"kind": kind, "res": f"{restype}{resnr + OFFSET}",
                         "p1": [round(float(x), 3) for x in p1],
                         "p2": [round(float(x), 3) for x in p2],
                         "detail": detail})

    for hb in site.hbonds_pdon + site.hbonds_ldon:
        add("hbond", hb.resnr, hb.restype,
            f"d={hb.distance_ad:.2f} angle={hb.angle:.0f} "
            f"{'protein-donor' if hb.protisdon else 'ligand-donor'}",
            hb.d.coords, hb.a.coords)
    for hc in site.hydrophobic_contacts:
        add("hydrophobic", hc.resnr, hc.restype, f"d={hc.distance:.2f}",
            hc.bsatom.coords, hc.ligatom.coords)
    for ps in site.pistacking:
        add("pistacking", ps.resnr, ps.restype, f"d={ps.distance:.2f} type={ps.type}",
            ps.proteinring.center, ps.ligandring.center)
    for sb in site.saltbridge_lneg + site.saltbridge_pneg:
        add("saltbridge", sb.resnr, sb.restype, f"d={sb.distance:.2f}",
            sb.positive.center, sb.negative.center)
    for pc in site.pication_laro + site.pication_paro:
        add("pication", pc.resnr, pc.restype, f"d={pc.distance:.2f}",
            pc.ring.center, pc.charge.center)
    for wb in site.water_bridges:
        add("waterbridge", wb.resnr, wb.restype,
            f"dAW={wb.distance_aw:.2f} dDW={wb.distance_dw:.2f}",
            wb.d.coords, wb.a.coords)
    for hal in site.halogen_bonds:
        add("halogen", hal.resnr, hal.restype, f"d={hal.distance:.2f}",
            hal.don.x.coords, hal.acc.o.coords)
    return rows, geom


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--receptor", default="model", choices=["model", "crystal"])
    ap.add_argument("--top", type=int, default=20)
    a = ap.parse_args()
    global OFFSET
    if a.receptor == "crystal":
        OFFSET = 0            # 1Q5K already uses UniProt numbering

    rec_pdb = os.path.join(ROOT, "05_docking", "receptor",
                           ("gsk3b_model_H.pdb" if a.receptor == "model" else "1q5kA_H.pdb"))
    scores = list(csv.DictReader(open(os.path.join(ROOT, "05_docking", a.receptor,
                                                   "scores.csv"), encoding="utf8")))
    scores.sort(key=lambda r: float(r["best_score"]))
    actives = [r for r in scores if r["set"] == "known_active"]
    hits = [r for r in scores if r["set"] != "known_active"][:a.top]
    todo = actives + hits
    print(f"profiling {len(todo)} compounds ({len(actives)} actives + top {len(hits)} library)")

    os.makedirs(CPX, exist_ok=True); os.makedirs(PLIPDIR, exist_ok=True)
    inter_rows, fp, geoms = [], {}, {}
    for r in todo:
        lid = r["id"]
        pose = os.path.join(ROOT, "05_docking", a.receptor, "poses", lid + ".pdbqt")
        cpx = os.path.join(CPX, f"{lid}_{a.receptor}_complex.pdb")
        build_complex(rec_pdb, pose, cpx)
        try:
            site = run_plip(cpx)
        except Exception as e:
            print(f"  !! PLIP failed for {lid}: {e}"); continue
        if site is None:
            print(f"  !! no binding site detected for {lid}"); continue
        rows, geom = collect(site)
        geoms[lid] = {"name": r["name"], "set": r["set"],
                      "score": float(r["best_score"]), "interactions": geom}
        fp[lid] = set()
        for kind, resnr, restype, detail in rows:
            inter_rows.append({"id": lid, "name": r["name"], "set": r["set"],
                               "score": r["best_score"], "kind": kind,
                               "resnum": resnr, "resname": restype, "detail": detail})
            fp[lid].add((kind, resnr, restype))
        print(f"  {lid:7s} {r['name'][:26]:26s} {float(r['best_score']):6.2f}  "
              f"{len(rows):2d} interactions")

    with open(os.path.join(PLIPDIR, "interactions.csv"), "w", newline="", encoding="utf8") as fh:
        w = csv.DictWriter(fh, fieldnames=list(inter_rows[0].keys()))
        w.writeheader(); w.writerows(inter_rows)

    # residue x compound fingerprint (any interaction type)
    residues = sorted({(rn, rt) for s in fp.values() for _, rn, rt in s})
    meta = {r["id"]: r for r in todo}
    with open(os.path.join(PLIPDIR, "fingerprint.csv"), "w", newline="", encoding="utf8") as fh:
        w = csv.writer(fh)
        w.writerow(["id", "name", "set", "score"] +
                   [f"{rt}{rn + OFFSET}" for rn, rt in residues])
        for lid, s in fp.items():
            present = {(rn, rt) for _, rn, rt in s}
            w.writerow([lid, meta[lid]["name"], meta[lid]["set"], meta[lid]["best_score"]] +
                       [1 if r in present else 0 for r in residues])
    json.dump(geoms, open(os.path.join(PLIPDIR, "interaction_geometry.json"), "w"), indent=1)
    print("\nwrote", os.path.join(PLIPDIR, "interactions.csv"),
          ", fingerprint.csv and interaction_geometry.json")


if __name__ == "__main__":
    main()

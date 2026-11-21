"""
Stage 3c - receptor preparation (both receptors, same frame, same protocol).

  model    = MODELLER model, aligned to the 1Q5K frame
  crystal  = 1Q5K chain A  (cross-check receptor)

Ten crystallographic waters that sit within 6 A of the co-crystal ligand are
kept in both receptors. The redocking test (07_redock_validation.py) shows
why: without them Vina never samples the crystallographic orientation of
AR-A014418 (best of 20 modes 4.2 A); with them the crystal pose is sampled
to 1.3 A. They are transplanted into the model unchanged - the model and the
crystal are in one frame, and the ATP-site backbone differs by 0.35 A RMSD,
so the water positions are transferable.

Protocol for each receptor:
  1. PDB2PQR 3.x : AMBER force field, PROPKA titration states at pH 7.4,
     hydrogens added and the H-bond network optimised.
  2. Meeko mk_prepare_receptor : rigid PDBQT, Gasteiger charges, merged
     non-polar hydrogens.

Outputs: 05_docking/receptor/{gsk3b_model,1q5kA}.pdbqt (+ _H.pdb, .pqr, .log)
"""
import os, subprocess, sys
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, ".."))
TPL = os.path.join(ROOT, "02_structure", "templates")
MOD = os.path.join(ROOT, "02_structure", "models")
REC = os.path.join(ROOT, "05_docking", "receptor")
MK = os.path.join(os.path.expanduser("~"), "AppData", "Roaming", "Python",
                  f"Python{sys.version_info.major}{sys.version_info.minor}",
                  "Scripts", "mk_prepare_receptor.exe")
WATER_CUTOFF = 6.0


def conserved_waters():
    from Bio.PDB import PDBParser
    c = PDBParser(QUIET=True).get_structure("c", os.path.join(TPL, "1Q5K.pdb"))[0]
    tmu = [r for r in c["A"] if r.get_resname() == "TMU"][0]
    tx = np.array([a.coord for a in tmu])
    out = []
    for r in c["A"]:
        if r.get_resname() != "HOH":
            continue
        if min(np.linalg.norm(a.coord - tx, axis=1).min() for a in r) < WATER_CUTOFF:
            out.append(r)
    return out


def write_with_waters(src_pdb, waters, out_pdb, chain="A"):
    """append the waters to a protein-only PDB as HETATM records"""
    lines = [l for l in open(src_pdb) if l.startswith(("ATOM", "TER"))]
    n = max(int(l[6:11]) for l in lines if l.startswith("ATOM"))
    w = []
    for i, res in enumerate(waters, 1):
        for a in res:
            n += 1
            x, y, z = a.coord
            w.append(f"HETATM{n:5d}  O   HOH {chain}{900+i:4d}    "
                     f"{x:8.3f}{y:8.3f}{z:8.3f}  1.00  0.00           O\n")
    open(out_pdb, "w").writelines(lines + w + ["END\n"])
    return out_pdb


def prepare(src_pdb, tag):
    h = os.path.join(REC, f"{tag}_H.pdb")
    pqr = os.path.join(REC, f"{tag}.pqr")
    subprocess.run([sys.executable, "-m", "pdb2pqr", "--ff=AMBER", "--with-ph=7.4",
                    "--titration-state-method=propka", "--keep-chain",
                    f"--pdb-output={h}", src_pdb, pqr], check=True,
                   capture_output=True, text=True)
    r = subprocess.run([MK, "--read_pdb", h, "-o", os.path.join(REC, tag), "-p"],
                       capture_output=True, text=True)
    pdbqt = os.path.join(REC, tag + ".pdbqt")
    if not os.path.exists(pdbqt):
        print(r.stdout[-2000:], r.stderr[-2000:]); raise SystemExit(f"{tag}: no pdbqt")
    n_h2o = sum(1 for l in open(pdbqt) if "HOH" in l)
    n_at = sum(1 for l in open(pdbqt) if l.startswith(("ATOM", "HETATM")))
    print(f"  {tag:14s} {n_at:5d} atoms, {n_h2o} water atoms -> {pdbqt}")


def main():
    os.makedirs(REC, exist_ok=True)
    wat = conserved_waters()
    print(f"conserved waters within {WATER_CUTOFF} A of the co-crystal ligand: {len(wat)}")
    src_model = write_with_waters(os.path.join(MOD, "gsk3b_model_aligned.pdb"), wat,
                                  os.path.join(MOD, "gsk3b_model_aligned_wat.pdb"))
    src_xtal = write_with_waters(os.path.join(TPL, "1Q5K_A.pdb"), wat,
                                 os.path.join(TPL, "1Q5K_A_wat.pdb"))
    prepare(src_model, "gsk3b_model")
    prepare(src_xtal, "1q5kA")


if __name__ == "__main__":
    main()

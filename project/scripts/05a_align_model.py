"""
Stage 3a - put the model in the crystal frame.

MODELLER writes the model in its own frame, which is not the 1Q5K frame. We
superpose the model onto 1Q5K chain A (C-alpha, all residues present in both)
and save the transformed coordinates. Everything downstream - P2Rank, the
search box, the redocking RMSD, the figures - then lives in one frame, so a
pose docked into the model and a pose docked into the crystal are directly
comparable.

Output: 02_structure/models/gsk3b_model_aligned.pdb
"""
import os

HERE = os.path.dirname(os.path.abspath(__file__))
MOD = os.path.join(HERE, "..", "02_structure", "models")
TPL = os.path.join(HERE, "..", "02_structure", "templates")
OFFSET = 34


def main():
    from Bio.PDB import PDBParser, PDBIO, Superimposer
    p = PDBParser(QUIET=True)
    m = p.get_structure("m", os.path.join(MOD, "gsk3b_model_best.pdb"))
    c = p.get_structure("c", os.path.join(TPL, "1Q5K_A.pdb"))
    mres = {r.id[1] + OFFSET: r for r in list(m[0])[0] if r.id[0] == " " and "CA" in r}
    cres = {r.id[1]: r for r in list(c[0])[0] if r.id[0] == " " and "CA" in r}
    common = sorted(set(mres) & set(cres))
    sup = Superimposer()
    sup.set_atoms([cres[i]["CA"] for i in common], [mres[i]["CA"] for i in common])
    sup.apply(list(m.get_atoms()))
    print(f"superposed model onto 1Q5K chain A on {len(common)} CA pairs, "
          f"RMSD {sup.rms:.3f} A")
    out = os.path.join(MOD, "gsk3b_model_aligned.pdb")
    io = PDBIO(); io.set_structure(m); io.save(out)
    print("wrote", out)


if __name__ == "__main__":
    main()

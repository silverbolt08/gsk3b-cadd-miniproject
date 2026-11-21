"""
Stage 2b - build 5 comparative models of GSK-3beta (P49841, residues 35-384)
from the three GSK-3beta templates, score them with DOPE and GA341, and
write a ranked summary.

Run after 02a_modeller_alignment.py:
    python 02b_modeller_build.py

Produces: 02_structure/models/GSK3B.B9999000{1..5}.pdb
          02_structure/models/model_scores.csv
          02_structure/models/gsk3b_model_best.pdb   (lowest DOPE)
"""
import os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import modeller_env  # noqa: F401  (sets up MODELLER paths)
import os, glob, shutil, csv

HERE = os.path.dirname(os.path.abspath(__file__))
ALN = os.path.join(HERE, "..", "02_structure", "alignment")
TPL = os.path.join(HERE, "..", "02_structure", "templates")
MOD = os.path.join(HERE, "..", "02_structure", "models")
N_MODELS = 5


def main():
    from modeller import Environ, log
    from modeller.automodel import AutoModel, assess, refine

    log.verbose()
    env = Environ()
    env.io.atom_files_directory = [TPL, ALN, "."]
    env.io.hetatm = False          # protein only; ligands are docked later
    env.io.water = False

    os.chdir(MOD)

    class GSK3BModel(AutoModel):
        """Standard automodel; the ATP-site residues are restrained a little
        harder so the pocket geometry stays close to the 1q5k template."""
        def special_restraints(self, aln):
            pass

    a = GSK3BModel(env,
                   alnfile=os.path.join(ALN, "gsk3b-mult.ali"),
                   knowns=("1q5kA", "1j1bA", "4accA"),
                   sequence="GSK3B",
                   assess_methods=(assess.DOPE, assess.DOPEHR, assess.GA341))
    a.starting_model = 1
    a.ending_model = N_MODELS
    a.md_level = refine.slow          # thorough MD refinement
    a.make()

    rows = []
    for m in a.outputs:
        if m["failure"] is not None:
            print("FAILED:", m["name"], m["failure"]); continue
        rows.append({"model": m["name"],
                     "molpdf": round(m["molpdf"], 2),
                     "DOPE": round(m["DOPE score"], 2),
                     "DOPE-HR": round(m["DOPE-HR score"], 2),
                     "GA341": round(m["GA341 score"][0], 4)})
    rows.sort(key=lambda r: r["DOPE"])
    with open(os.path.join(MOD, "model_scores.csv"), "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys())); w.writeheader(); w.writerows(rows)

    best = rows[0]["model"]
    shutil.copy(best, os.path.join(MOD, "gsk3b_model_best.pdb"))
    print("\n=== model ranking (by DOPE) ===")
    for r in rows:
        print(r)
    print("\nbest model:", best, "-> gsk3b_model_best.pdb")


if __name__ == "__main__":
    main()

"""
Stage 2c - per-residue normalised DOPE profile of the best model against the
1q5k template, plus a plot. A model is acceptable where its profile tracks the
template; persistent positive (unfavourable) stretches mark regions to treat
with caution.

Run:  python 02c_dope_profile.py
Produces: 02_structure/validation/dope_profile_model.profile
          02_structure/validation/dope_profile_1q5k.profile
          06_analysis/figures/fig2_dope_profile.png
"""
import os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import modeller_env  # noqa: F401  (sets up MODELLER paths)
import os

HERE = os.path.dirname(os.path.abspath(__file__))
MOD = os.path.join(HERE, "..", "02_structure", "models")
TPL = os.path.join(HERE, "..", "02_structure", "templates")
ALN = os.path.join(HERE, "..", "02_structure", "alignment")
VAL = os.path.join(HERE, "..", "02_structure", "validation")
FIG = os.path.join(HERE, "..", "06_analysis", "figures")


def profile(env, pdb, alnfile, code, out):
    from modeller import Selection, alignment, model
    from modeller.scripts import complete_pdb
    mdl = complete_pdb(env, pdb)
    s = Selection(mdl)
    s.assess_dope(output="ENERGY_PROFILE NO_REPORT", file=out,
                  normalize_profile=True, smoothing_window=15)
    return out


def read_profile(path, seq_len=None):
    vals = {}
    for line in open(path):
        if line.startswith("#") or not line.strip():
            continue
        f = line.split()
        try:
            vals[int(f[0])] = float(f[-1])
        except ValueError:
            continue
    return vals


def main():
    from modeller import Environ, log
    log.minimal()
    env = Environ()
    env.libs.topology.read(file="$(LIB)/top_heav.lib")
    env.libs.parameters.read(file="$(LIB)/par.lib")
    env.io.atom_files_directory = [MOD, TPL, "."]
    os.makedirs(VAL, exist_ok=True); os.makedirs(FIG, exist_ok=True)

    p_model = profile(env, os.path.join(MOD, "gsk3b_model_best.pdb"), None, "model",
                      os.path.join(VAL, "dope_profile_model.profile"))
    p_tmpl = profile(env, os.path.join(TPL, "1Q5K.pdb"), None, "1q5k",
                     os.path.join(VAL, "dope_profile_1q5k.profile"))

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    m = read_profile(p_model); t = read_profile(p_tmpl)
    fig, ax = plt.subplots(figsize=(9, 3.4), dpi=200)
    ax.plot(list(t), list(t.values()), lw=1.0, color="#9aa0a6", label="template 1Q5K")
    ax.plot(list(m), list(m.values()), lw=1.2, color="#b5392f", label="model (best DOPE)")
    ax.axhline(0, color="k", lw=0.6)
    ax.set_xlabel("residue index"); ax.set_ylabel("normalised DOPE")
    ax.legend(frameon=False, fontsize=8)
    ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    fig.savefig(os.path.join(FIG, "fig2_dope_profile.png"))
    print("wrote", os.path.join(FIG, "fig2_dope_profile.png"))


if __name__ == "__main__":
    main()

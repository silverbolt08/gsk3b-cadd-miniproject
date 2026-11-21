"""
Stage 2d - validate the comparative model.

1. Extract chain A (protein only) from the 1Q5K crystal -> 1Q5K_A.pdb
2. Per-residue normalised DOPE profile for the model and for 1Q5K chain A,
   plotted on the UniProt numbering so they can be read against each other.
3. C-alpha RMSD of the model to the crystal over all residues present in both,
   and separately over the ATP-site residues.
4. Ramachandran analysis using the Top8000 kernel density shipped with
   RamachanDraw. Favoured / allowed contours are taken the MolProbity way:
   the iso-density level enclosing 98% (favoured) and 99.95% (allowed) of the
   reference density mass.

Outputs: 02_structure/validation/validation_summary.json / .md
         06_analysis/figures/fig2_dope_profile.png
         06_analysis/figures/fig3_ramachandran.png
"""
import os, sys, json
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import modeller_env  # noqa: F401
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
MOD = os.path.join(HERE, "..", "02_structure", "models")
TPL = os.path.join(HERE, "..", "02_structure", "templates")
VAL = os.path.join(HERE, "..", "02_structure", "validation")
FIG = os.path.join(HERE, "..", "06_analysis", "figures")
MODEL = os.path.join(MOD, "gsk3b_model_best.pdb")
OFFSET = 34          # model residue 1 == UniProt residue 35

# ATP-site residues (UniProt numbering) from the GSK-3beta literature
ATP_SITE = [62, 63, 64, 65, 66, 67, 68, 69, 70, 85, 97, 110, 132, 133, 134,
            135, 136, 138, 141, 185, 186, 188, 199, 200]


def extract_chain_a():
    from Bio.PDB import PDBParser, PDBIO, Select
    out = os.path.join(TPL, "1Q5K_A.pdb")
    p = PDBParser(QUIET=True)
    st = p.get_structure("1q5k", os.path.join(TPL, "1Q5K.pdb"))

    class Sel(Select):
        def accept_chain(self, c):  return c.id == "A"
        def accept_residue(self, r): return r.id[0] == " "
        def accept_atom(self, a):   return a.element != "H" and not a.is_disordered() or a.get_altloc() in (" ", "A")

    io = PDBIO(); io.set_structure(st); io.save(out, Sel())
    return out


def dope_profiles(crystal_a):
    from modeller import Environ, Selection, log
    from modeller.scripts import complete_pdb
    log.none()
    env = Environ()
    env.libs.topology.read(file="$(LIB)/top_heav.lib")
    env.libs.parameters.read(file="$(LIB)/par.lib")
    out = {}
    for tag, pdb in (("model", MODEL), ("1q5k", crystal_a)):
        mdl = complete_pdb(env, pdb)
        f = os.path.join(VAL, f"dope_profile_{tag}.profile")
        Selection(mdl).assess_dope(output="ENERGY_PROFILE NO_REPORT", file=f,
                                   normalize_profile=True, smoothing_window=15)
        out[tag] = f
    return out


def read_profile(path, resnums):
    """profile rows are in file order; map them onto the given residue numbers"""
    vals = []
    for line in open(path):
        if line.startswith("#") or not line.strip():
            continue
        f = line.split()
        try:
            vals.append(float(f[-1]))
        except ValueError:
            continue
    n = min(len(vals), len(resnums))
    return resnums[:n], vals[:n]


def resnums_of(pdb, offset=0):
    from Bio.PDB import PDBParser
    st = PDBParser(QUIET=True).get_structure("x", pdb)
    ch = list(st[0])[0]
    return [r.id[1] + offset for r in ch if r.id[0] == " "]


def ca_rmsd(crystal_a):
    from Bio.PDB import PDBParser, Superimposer
    p = PDBParser(QUIET=True)
    m = p.get_structure("m", MODEL)[0]
    c = p.get_structure("c", crystal_a)[0]
    mch, cch = list(m)[0], list(c)[0]
    mres = {r.id[1] + OFFSET: r for r in mch if r.id[0] == " " and "CA" in r}
    cres = {r.id[1]: r for r in cch if r.id[0] == " " and "CA" in r}
    common = sorted(set(mres) & set(cres))
    fixed = [cres[i]["CA"] for i in common]
    moving = [mres[i]["CA"] for i in common]
    sup = Superimposer(); sup.set_atoms(fixed, moving)
    rms_all = sup.rms
    # apply the transform, then measure the pocket on its own
    sup.apply([a for r in mch for a in r])
    pocket = [i for i in common if i in ATP_SITE]
    d = [np.linalg.norm(mres[i]["CA"].coord - cres[i]["CA"].coord) for i in pocket]
    # all-heavy-atom RMSD over pocket side chains
    heavy = []
    for i in pocket:
        for at in cres[i]:
            if at.element == "H":
                continue
            if at.get_id() in mres[i]:
                heavy.append(np.linalg.norm(mres[i][at.get_id()].coord - at.coord))
    return dict(n_common=len(common), ca_rmsd_all=round(float(rms_all), 3),
                ca_rmsd_pocket=round(float(np.sqrt(np.mean(np.square(d)))), 3),
                heavy_rmsd_pocket=round(float(np.sqrt(np.mean(np.square(heavy)))), 3),
                n_pocket=len(pocket))


def ramachandran(pdb, label, ax=None):
    from ramachandraw.utils import get_phi_psi
    from importlib import resources
    z = np.fromfile(str(resources.files("ramachandraw") / "kde.dat")).reshape(100, 100)
    flat = np.sort(z.ravel())[::-1]
    cum = np.cumsum(flat) / flat.sum()
    lvl_fav = flat[np.searchsorted(cum, 0.98)]
    lvl_all = flat[np.searchsorted(cum, 0.9995)]
    angles = get_phi_psi(pdb_filepath=pdb)
    edges = np.linspace(-180, 180, 101)
    fav = allowed = outlier = 0
    outliers = []
    for key, (phi, psi) in angles.items():
        if phi is None or psi is None:
            continue
        i = min(np.searchsorted(edges, phi) - 1, 99)
        j = min(np.searchsorted(edges, psi) - 1, 99)
        d = z[max(i, 0), max(j, 0)]
        if d >= lvl_fav:
            fav += 1
        elif d >= lvl_all:
            allowed += 1
        else:
            outlier += 1; outliers.append(key)
    tot = fav + allowed + outlier
    return dict(label=label, n=tot,
                favoured_pct=round(100 * fav / tot, 1),
                allowed_pct=round(100 * allowed / tot, 1),
                outlier_pct=round(100 * outlier / tot, 1),
                outliers=outliers), angles


def main():
    os.makedirs(VAL, exist_ok=True); os.makedirs(FIG, exist_ok=True)
    crystal_a = extract_chain_a()
    print("chain A extracted:", crystal_a)

    prof = dope_profiles(crystal_a)
    mnums = resnums_of(MODEL, OFFSET)
    cnums = resnums_of(crystal_a, 0)
    mx, my = read_profile(prof["model"], mnums)
    cx, cy = read_profile(prof["1q5k"], cnums)

    import matplotlib; matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(figsize=(9, 3.2), dpi=200)
    ax.plot(cx, cy, lw=1.0, color="#9aa0a6", label="template 1Q5K chain A")
    ax.plot(mx, my, lw=1.2, color="#b5392f", label="model 5 (best DOPE)")
    for r in (120, 121, 287, 288, 289):
        ax.axvline(r, color="#4a7fb5", lw=3, alpha=0.18)
    ax.axhline(0, color="k", lw=0.5)
    ax.set_xlabel("residue (UniProt P49841 numbering)")
    ax.set_ylabel("normalised DOPE")
    ax.set_title("Per-residue DOPE profile (blue bands = loops rebuilt by MODELLER)", fontsize=9)
    ax.legend(frameon=False, fontsize=8)
    ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout(); fig.savefig(os.path.join(FIG, "fig2_dope_profile.png")); plt.close(fig)

    rms = ca_rmsd(crystal_a)
    print("RMSD:", rms)

    rama_m, ang_m = ramachandran(MODEL, "model")
    rama_c, ang_c = ramachandran(crystal_a, "1Q5K_A")
    print("Ramachandran model:", {k: v for k, v in rama_m.items() if k != "outliers"})
    print("Ramachandran 1q5k :", {k: v for k, v in rama_c.items() if k != "outliers"})

    # Ramachandran figure
    from importlib import resources
    z = np.fromfile(str(resources.files("ramachandraw") / "kde.dat")).reshape(100, 100)
    fig, ax = plt.subplots(figsize=(4.6, 4.4), dpi=200)
    ax.imshow(np.log10(np.rot90(z)), cmap="Blues", extent=(-180, 180, -180, 180), alpha=0.9)
    ax.contour(np.rot90(np.fliplr(z)), colors="k", linewidths=0.4,
               levels=[10 ** i for i in range(-7, 0)], extent=[-180, 180, -180, 180], alpha=0.4)
    xs = [v[0] for v in ang_m.values() if v[0] is not None]
    ys = [v[1] for v in ang_m.values() if v[1] is not None]
    ax.scatter(xs, ys, s=3, c="#b5392f", lw=0)
    ax.set(aspect="equal", xlim=(-180, 180), ylim=(-180, 180),
           xlabel="phi (deg)", ylabel="psi (deg)",
           title=f"Model: {rama_m['favoured_pct']}% favoured, "
                 f"{rama_m['outlier_pct']}% outliers", )
    ax.title.set_fontsize(9)
    fig.tight_layout(); fig.savefig(os.path.join(FIG, "fig3_ramachandran.png")); plt.close(fig)

    summary = {"model": os.path.basename(MODEL), "rmsd": rms,
               "ramachandran_model": rama_m, "ramachandran_crystal": rama_c}
    json.dump(summary, open(os.path.join(VAL, "validation_summary.json"), "w"), indent=1)

    with open(os.path.join(VAL, "validation_summary.md"), "w", encoding="utf8") as fh:
        fh.write("# Model validation\n\n")
        fh.write(f"- residues compared with 1Q5K chain A: {rms['n_common']}\n")
        fh.write(f"- C-alpha RMSD (all common residues): **{rms['ca_rmsd_all']} A**\n")
        fh.write(f"- C-alpha RMSD (ATP-site residues, n={rms['n_pocket']}): "
                 f"**{rms['ca_rmsd_pocket']} A**\n")
        fh.write(f"- all-heavy-atom RMSD over ATP-site residues: "
                 f"**{rms['heavy_rmsd_pocket']} A**\n")
        fh.write(f"- Ramachandran (model): {rama_m['favoured_pct']}% favoured, "
                 f"{rama_m['allowed_pct']}% allowed, {rama_m['outlier_pct']}% outliers\n")
        fh.write(f"- Ramachandran (1Q5K chain A, for reference): "
                 f"{rama_c['favoured_pct']}% favoured, {rama_c['allowed_pct']}% allowed, "
                 f"{rama_c['outlier_pct']}% outliers\n")
    print("wrote", os.path.join(VAL, "validation_summary.md"))


if __name__ == "__main__":
    main()

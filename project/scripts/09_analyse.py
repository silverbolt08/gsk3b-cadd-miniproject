"""
Stage 5/6 - rank the screen, measure how well it separates the known actives
from the background library, and draw the figures.

Metrics
  ROC-AUC            over the 139 docked compounds, actives = positives.
                     The library compounds are a background set, not
                     confirmed inactives, so this is an optimistic-but-honest
                     retrospective check, not a benchmark.
  EF1% / EF5% / EF10% enrichment factor at the usual cut-offs.
  BEDROC(alpha=20)   early-recognition weighted metric.
  model vs crystal   Pearson/Spearman correlation of the two score sets, i.e.
                     how much the homology model changes the ranking.

Outputs: 06_analysis/ranking_model.csv, enrichment.json,
         figures fig4..fig7
"""
import csv, json, math, os
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, ".."))
ANA = os.path.join(ROOT, "06_analysis")
FIG = os.path.join(ANA, "figures")
C_ACT, C_DRUG, C_NP = "#b5392f", "#4a7fb5", "#2f8f6f"



def cns_mpo(mw, clogp, tpsa, hbd):
    """Four of the six CNS-MPO desirability terms (Wager et al., ACS Chem
    Neurosci 2010): MW, clogP, TPSA and HBD. clogD and pKa are left out -
    we have no measured pKa - so the score runs 0-4, not 0-6."""
    def hump(x, lo_good, hi_good):          # 1 inside the window, ramps to 0
        if x <= lo_good: return 1.0
        if x >= hi_good: return 0.0
        return (hi_good - x) / (hi_good - lo_good)

    def window(x, a, b, c, d):              # 0 below a, 1 between b and c, 0 above d
        if x <= a or x >= d: return 0.0
        if b <= x <= c: return 1.0
        if x < b: return (x - a) / (b - a)
        return (d - x) / (d - c)

    return round(hump(mw, 360.0, 500.0) + hump(clogp, 3.0, 5.0)
                 + window(tpsa, 20.0, 40.0, 90.0, 120.0) + hump(hbd, 0.5, 3.5), 2)


def load(tag):
    p = os.path.join(ROOT, "05_docking", tag, "scores.csv")
    rows = list(csv.DictReader(open(p, encoding="utf8")))
    for r in rows:
        r["best_score"] = float(r["best_score"])
        r["ligand_efficiency"] = float(r["ligand_efficiency"])
        r["mw"] = float(r["mw"]); r["clogp"] = float(r["clogp"])
        r["tpsa"] = float(r["tpsa"]); r["is_active"] = r["set"] == "known_active"
        r["cns_mpo4"] = cns_mpo(r["mw"], r["clogp"], r["tpsa"], int(r["hbd"]))
    rows.sort(key=lambda r: r["best_score"])
    for i, r in enumerate(rows, 1):
        r["rank"] = i
    return rows


def roc_auc(rows):
    y = np.array([r["is_active"] for r in rows], dtype=float)
    n_pos, n_neg = y.sum(), len(y) - y.sum()
    # rows are already sorted best-first
    tp = np.cumsum(y); fp = np.cumsum(1 - y)
    tpr = tp / n_pos; fpr = fp / n_neg
    auc = np.trapezoid(np.concatenate([[0], tpr]), np.concatenate([[0], fpr]))
    return float(auc), fpr, tpr


def enrichment(rows, frac):
    n = max(1, int(round(len(rows) * frac)))
    top = rows[:n]
    hits = sum(r["is_active"] for r in top)
    total_act = sum(r["is_active"] for r in rows)
    return round((hits / n) / (total_act / len(rows)), 2), hits, n


def bedroc(rows, alpha=20.0):
    n = len(rows); act = [i for i, r in enumerate(rows, 1) if r["is_active"]]
    N, Na = n, len(act)
    ra = Na / N
    s = sum(math.exp(-alpha * i / N) for i in act)
    rie = s / (Na / N * (1 - math.exp(-alpha)) / (math.exp(alpha / N) - 1))
    return round(rie * ra * math.sinh(alpha / 2) / (math.cosh(alpha / 2) - math.cosh(alpha / 2 - alpha * ra))
                 + 1 / (1 - math.exp(alpha * (1 - ra))), 3)


def main():
    os.makedirs(FIG, exist_ok=True)
    model = load("model")
    try:
        crystal = load("crystal")
        if len(crystal) < 0.9 * len(model):
            print(f"crystal run incomplete ({len(crystal)}/{len(model)}) - "
                  "model-only analysis")
            crystal = []
    except Exception:
        crystal = []
    cm = {r["id"]: r for r in crystal}

    for r in model:
        r["crystal_score"] = cm[r["id"]]["best_score"] if r["id"] in cm else ""
        r["crystal_rank"] = cm[r["id"]]["rank"] if r["id"] in cm else ""

    cols = ["rank", "id", "name", "set", "chembl_id", "potency_nM", "best_score",
            "ligand_efficiency", "cns_mpo4", "crystal_score", "crystal_rank", "mw",
            "clogp", "tpsa", "hbd", "hba", "rotb", "heavy_atoms", "mode_scores"]
    with open(os.path.join(ANA, "ranking_model.csv"), "w", newline="", encoding="utf8") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, extrasaction="ignore")
        w.writeheader(); w.writerows(model)

    auc, fpr, tpr = roc_auc(model)
    if crystal:
        auc_c, fpr_c, tpr_c = roc_auc(crystal)
    else:
        auc_c, fpr_c, tpr_c = float("nan"), None, None
    stats = {"n_compounds": len(model),
             "n_actives": sum(r["is_active"] for r in model),
             "roc_auc_model": round(auc, 3),
             "roc_auc_crystal": round(auc_c, 3) if crystal else None,
             "bedroc20_model": bedroc(model),
             "bedroc20_crystal": bedroc(crystal) if crystal else None}
    for f in (0.01, 0.05, 0.10):
        ef, hits, n = enrichment(model, f)
        stats[f"EF{int(f*100)}pct_model"] = {"EF": ef, "actives_found": hits, "n_top": n}

    from scipy.stats import pearsonr, spearmanr
    common = [r for r in model if r["crystal_score"] != ""]
    if common:
        x = np.array([r["best_score"] for r in common])
        y = np.array([r["crystal_score"] for r in common])
        stats["model_vs_crystal_pearson"] = round(float(pearsonr(x, y)[0]), 3)
        stats["model_vs_crystal_spearman"] = round(float(spearmanr(x, y)[0]), 3)
        stats["mean_score_shift_model_minus_crystal"] = round(float((x - y).mean()), 2)
    else:
        x = y = np.array([])
        stats["model_vs_crystal_pearson"] = None
        stats["model_vs_crystal_spearman"] = None
        stats["mean_score_shift_model_minus_crystal"] = None
    # does the score track measured potency among the actives?
    import math
    act = [r for r in model if r["is_active"] and r["potency_nM"]]
    px = [-math.log10(float(r["potency_nM"]) * 1e-9) for r in act]
    py = [r["best_score"] for r in act]
    stats["actives_n_with_potency"] = len(act)
    stats["score_vs_pchembl_pearson"] = round(float(pearsonr(px, py)[0]), 3)
    stats["score_vs_pchembl_spearman"] = round(float(spearmanr(px, py)[0]), 3)

    # structural-alert screen on the compounds we might shortlist
    try:
        from rdkit import Chem, RDLogger
        from rdkit.Chem import FilterCatalog
        RDLogger.DisableLog("rdApp.*")
        man = {r["id"]: r for r in csv.DictReader(
            open(os.path.join(ROOT, "04_ligands", "ligand_manifest.csv"), encoding="utf8"))}
        params = FilterCatalog.FilterCatalogParams()
        for cat in (FilterCatalog.FilterCatalogParams.FilterCatalogs.PAINS,
                    FilterCatalog.FilterCatalogParams.FilterCatalogs.BRENK):
            params.AddCatalog(cat)
        fc = FilterCatalog.FilterCatalog(params)
        alerts = []
        for r in model[:40]:
            smi = man[r["id"]]["smiles_ph74"]
            m = Chem.MolFromSmiles(smi)
            hits = [e.GetDescription() for e in fc.GetMatches(m)] if m else ["unparsed"]
            alerts.append({"rank": r["rank"], "id": r["id"], "name": r["name"],
                           "set": r["set"], "score": r["best_score"],
                           "n_alerts": len(hits), "alerts": "; ".join(hits[:4])})
        with open(os.path.join(ANA, "structural_alerts.csv"), "w", newline="",
                  encoding="utf8") as fh:
            w = csv.DictWriter(fh, fieldnames=list(alerts[0].keys()))
            w.writeheader(); w.writerows(alerts)
        stats["alerts_in_top40"] = sum(1 for a in alerts if a["n_alerts"])
    except Exception as e:
        print("alert screen failed:", e)

    json.dump(stats, open(os.path.join(ANA, "enrichment.json"), "w"), indent=1)
    print(json.dumps(stats, indent=1))

    import matplotlib; matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    # fig4 - score distribution
    fig, ax = plt.subplots(figsize=(7.2, 3.4), dpi=200)
    lo = min(r["best_score"] for r in model)
    bins = np.arange(np.floor(lo * 2) / 2, -3.4, 0.4)
    for s, c, lab in (("A_approved_drug", C_DRUG, "approved drugs"),
                      ("B_natural_product", C_NP, "natural products")):
        v = [r["best_score"] for r in model if r["set"] == s]
        ax.hist(v, bins=bins, alpha=0.62, color=c, label=f"{lab} ({len(v)})",
                edgecolor="white", lw=0.4)
    for r in model:
        if r["is_active"]:
            ax.axvline(r["best_score"], color=C_ACT, lw=1.0, alpha=0.9)
    ax.plot([], [], color=C_ACT, lw=1.0,
            label=f"known actives ({sum(r['is_active'] for r in model)})")
    ax.set_xlabel("Vina score (kcal/mol)"); ax.set_ylabel("compounds")
    ax.set_title("Docking score distribution, GSK-3$\\beta$ homology model", fontsize=10)
    ax.legend(frameon=False, fontsize=8); ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout(); fig.savefig(os.path.join(FIG, "fig4_score_distribution.png")); plt.close(fig)

    # fig5 - ROC
    fig, ax = plt.subplots(figsize=(4.0, 3.8), dpi=200)
    ax.plot(fpr, tpr, color=C_ACT, lw=1.6, label=f"model  AUC {auc:.2f}")
    if crystal:
        ax.plot(fpr_c, tpr_c, color=C_DRUG, lw=1.2, ls="--",
                label=f"1Q5K crystal  AUC {auc_c:.2f}")
    ax.plot([0, 1], [0, 1], color="#999", lw=0.8, ls=":")
    ax.set_xlabel("false positive rate"); ax.set_ylabel("true positive rate")
    ax.set_title("Recovery of the 11 known actives", fontsize=9)
    ax.legend(frameon=False, fontsize=8, loc="lower right")
    ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout(); fig.savefig(os.path.join(FIG, "fig5_roc.png")); plt.close(fig)

    # fig6 - model vs crystal scores
    fig, ax = plt.subplots(figsize=(4.2, 4.0), dpi=200)
    for s, c, lab in (("A_approved_drug", C_DRUG, "approved drug"),
                      ("B_natural_product", C_NP, "natural product"),
                      ("known_active", C_ACT, "known active")):
        xs = [r["best_score"] for r in common if r["set"] == s]
        ys = [r["crystal_score"] for r in common if r["set"] == s]
        ax.scatter(xs, ys, s=14, c=c, alpha=0.75, lw=0, label=lab)
    if not common:
        plt.close(fig)
    lim = ([min(x.min(), y.min()) - 0.4, max(x.max(), y.max()) + 0.4]
           if common else [-12, -4])
    ax.plot(lim, lim, color="#999", lw=0.8, ls=":")
    ax.set_xlim(lim); ax.set_ylim(lim)
    ax.set_xlabel("score on homology model (kcal/mol)")
    ax.set_ylabel("score on 1Q5K crystal (kcal/mol)")
    ax.set_title(f"r = {stats['model_vs_crystal_pearson']}, "
                 f"$\\rho$ = {stats['model_vs_crystal_spearman']}", fontsize=9)
    ax.legend(frameon=False, fontsize=7.5, loc="upper left")
    ax.spines[["top", "right"]].set_visible(False)
    if common:
        fig.tight_layout(); fig.savefig(os.path.join(FIG, "fig6_model_vs_crystal.png"))
    plt.close(fig)

    # fig7 - score vs ligand efficiency for the top 30
    top = model[:30]
    fig, ax = plt.subplots(figsize=(7.4, 3.8), dpi=200)
    for r in top:
        c = C_ACT if r["is_active"] else (C_DRUG if r["set"] == "A_approved_drug" else C_NP)
        ax.scatter(r["best_score"], r["ligand_efficiency"], s=26, c=c, lw=0)
        ax.annotate(r["name"][:16], (r["best_score"], r["ligand_efficiency"]),
                    fontsize=5.5, xytext=(2, 2), textcoords="offset points")
    ax.set_xlabel("Vina score (kcal/mol)")
    ax.set_ylabel("ligand efficiency (kcal/mol per heavy atom)")
    ax.set_title("Top 30 of the screen", fontsize=10)
    ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout(); fig.savefig(os.path.join(FIG, "fig7_top30_le.png")); plt.close(fig)

    print("\nTop 20 on the model:")
    for r in model[:20]:
        tag = "ACTIVE" if r["is_active"] else ("drug" if r["set"] == "A_approved_drug" else "NP")
        print(f"  {r['rank']:3d} {r['id']:7s} {r['name'][:28]:28s} {r['best_score']:6.2f}  "
              f"LE {r['ligand_efficiency']:.3f}  [{tag}]")


if __name__ == "__main__":
    main()

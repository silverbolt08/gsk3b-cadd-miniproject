"""
Stage 6 - interaction fingerprint heat map.

Rows are compounds (known actives first, then the top library hits), columns
are the ATP-site residues PLIP reports contacts with. A filled cell means at
least one interaction of any type; H-bonds are marked separately so the
hinge pattern (Asp133 / Val135) is readable at a glance.

Output: 06_analysis/figures/fig8_fingerprint.png
"""
import csv, os
from collections import defaultdict
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, ".."))
PLIPDIR = os.path.join(ROOT, "06_analysis", "plip")
FIG = os.path.join(ROOT, "06_analysis", "figures")
OFFSET = 34
HINGE = {133, 135, 134}


def main():
    rows = list(csv.DictReader(open(os.path.join(PLIPDIR, "interactions.csv"), encoding="utf8")))
    comps, order = {}, []
    any_i, hb_i = defaultdict(set), defaultdict(set)
    for r in rows:
        lid = r["id"]
        if lid not in comps:
            comps[lid] = (r["name"], r["set"], float(r["score"]))
            order.append(lid)
        res = int(r["resnum"]) + OFFSET
        any_i[lid].add((res, r["resname"]))
        if r["kind"] in ("hbond", "waterbridge", "saltbridge"):
            hb_i[lid].add((res, r["resname"]))

    residues = sorted({x for s in any_i.values() for x in s})
    order.sort(key=lambda l: (comps[l][1] != "known_active", comps[l][2]))

    M = np.zeros((len(order), len(residues)))
    for i, lid in enumerate(order):
        for j, res in enumerate(residues):
            if res in hb_i[lid]:
                M[i, j] = 2
            elif res in any_i[lid]:
                M[i, j] = 1

    import matplotlib; matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.colors import ListedColormap
    cmap = ListedColormap(["#f2f2f4", "#b9cde4", "#2f6fd0"])
    h = max(3.2, 0.21 * len(order) + 1.6)
    fig, ax = plt.subplots(figsize=(0.30 * len(residues) + 3.2, h), dpi=200)
    ax.imshow(M, cmap=cmap, aspect="auto", vmin=0, vmax=2)
    ax.set_xticks(range(len(residues)))
    ax.set_xticklabels([f"{rt.title()}{rn}" for rn, rt in residues],
                       rotation=90, fontsize=7)
    for t, (rn, _) in zip(ax.get_xticklabels(), residues):
        if rn in HINGE:
            t.set_color("#b5392f"); t.set_fontweight("bold")
    ax.set_yticks(range(len(order)))
    ax.set_yticklabels([f"{comps[l][0][:22]}  ({comps[l][2]:.1f})" for l in order], fontsize=7)
    for t, l in zip(ax.get_yticklabels(), order):
        if comps[l][1] == "known_active":
            t.set_color("#b5392f"); t.set_fontweight("bold")
    ax.set_xticks(np.arange(-.5, len(residues), 1), minor=True)
    ax.set_yticks(np.arange(-.5, len(order), 1), minor=True)
    ax.grid(which="minor", color="white", lw=0.8)
    ax.tick_params(which="minor", length=0)
    n_act = sum(1 for l in order if comps[l][1] == "known_active")
    ax.axhline(n_act - 0.5, color="#b5392f", lw=1.2)
    ax.set_title("PLIP interaction fingerprint - known actives (red, above the line) "
                 "vs top library hits\ndark = H-bond / water bridge / salt bridge, "
                 "light = hydrophobic or $\\pi$ contact; hinge residues in red",
                 fontsize=8.5, loc="left")
    fig.tight_layout()
    out = os.path.join(FIG, "fig8_fingerprint.png")
    fig.savefig(out); print("wrote", out)

    # how often is the hinge used?
    import json
    from collections import Counter
    print("\nhinge (Asp133/Val135) H-bond usage:")
    summary = {"n_actives": 0, "actives_hinge": 0, "n_hits": 0, "hits_hinge": 0,
               "per_compound": {}}
    freq = Counter()
    for l in order:
        hinge_hb = {rn for rn, _ in hb_i[l]} & {133, 135}
        is_act = comps[l][1] == "known_active"
        summary["n_actives" if is_act else "n_hits"] += 1
        if hinge_hb:
            summary["actives_hinge" if is_act else "hits_hinge"] += 1
        summary["per_compound"][l] = {"name": comps[l][0], "set": comps[l][1],
                                      "score": comps[l][2],
                                      "hinge": sorted(hinge_hb)}
        for rn, rt in any_i[l]:
            freq[f"{rt.title()}{rn}"] += 1
        print(f"  {comps[l][0][:26]:26s} {comps[l][1][:14]:14s} "
              f"{comps[l][2]:6.2f}  {'yes ' + str(sorted(hinge_hb)) if hinge_hb else 'no'}")
    summary["top_residues"] = [r for r, _ in freq.most_common(6)]
    json.dump(summary, open(os.path.join(ROOT, "06_analysis", "hinge_summary.json"), "w"),
              indent=1)
    print("\nmost-contacted residues:", summary["top_residues"])


if __name__ == "__main__":
    main()

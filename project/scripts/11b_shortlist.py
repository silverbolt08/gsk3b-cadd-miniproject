"""
Stage 6 - apply the selection criteria mechanically, so the shortlist is
reproducible rather than a matter of taste.

A library compound is eligible if it:
  1. scores at least as well as the median positive control;
  2. makes at least one H-bond / water bridge to the hinge (Asp133 or Val135),
     the anchor the known actives use;
  3. has ligand efficiency >= 0.30 kcal/mol per heavy atom;
  4. has CNS-MPO (4 terms) >= 3.0, i.e. a plausible shot at the brain.

Compounds carrying PAINS/Brenk structural alerts are not excluded
automatically - the alerts are printed next to each candidate so the choice
can be made with them in view.

Output: 06_analysis/shortlist.csv  (every library compound, with pass/fail
        on each criterion)
"""
import csv, json, os, statistics

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, ".."))
ANA = os.path.join(ROOT, "06_analysis")

LE_MIN = 0.30
MPO_MIN = 3.0


def main():
    ranking = list(csv.DictReader(open(os.path.join(ANA, "ranking_model.csv"), encoding="utf8")))
    hinge = json.load(open(os.path.join(ANA, "hinge_summary.json")))["per_compound"]
    alerts = {}
    ap = os.path.join(ANA, "structural_alerts.csv")
    if os.path.exists(ap):
        alerts = {r["id"]: r for r in csv.DictReader(open(ap, encoding="utf8"))}

    ctrl = [float(r["best_score"]) for r in ranking if r["set"] == "known_active"]
    med_ctrl = statistics.median(ctrl)
    print(f"median positive-control score: {med_ctrl:.2f} kcal/mol "
          f"(range {min(ctrl):.2f} to {max(ctrl):.2f})")

    rows = []
    for r in ranking:
        if r["set"] == "known_active":
            continue
        h = hinge.get(r["id"], {}).get("hinge", [])
        sc = float(r["best_score"]); le = float(r["ligand_efficiency"])
        mpo = float(r["cns_mpo4"])
        c1 = sc <= med_ctrl
        c2 = bool(h)
        c3 = le >= LE_MIN
        c4 = mpo >= MPO_MIN
        rows.append({"rank": int(r["rank"]), "id": r["id"], "name": r["name"],
                     "set": r["set"], "score": sc, "LE": le, "cns_mpo4": mpo,
                     "hinge_hbond": ",".join(str(x) for x in h) or "-",
                     "c1_score": int(c1), "c2_hinge": int(c2), "c3_LE": int(c3),
                     "c4_cns": int(c4), "n_pass": sum([c1, c2, c3, c4]),
                     "profiled": r["id"] in hinge,
                     "alerts": alerts.get(r["id"], {}).get("alerts", "")[:70]})

    rows.sort(key=lambda x: (-x["n_pass"], x["score"]))
    with open(os.path.join(ANA, "shortlist.csv"), "w", newline="", encoding="utf8") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys())); w.writeheader(); w.writerows(rows)

    print("\ncompounds passing all four criteria:")
    any_pass = False
    for x in rows:
        if x["n_pass"] == 4:
            any_pass = True
            print(f"  rank {x['rank']:3d}  {x['name'][:26]:26s} {x['score']:6.2f}  "
                  f"LE {x['LE']:.3f}  MPO {x['cns_mpo4']:.1f}  hinge {x['hinge_hbond']:8s}"
                  f"  {x['alerts'][:46]}")
    if not any_pass:
        print("  none")
    print("\nnext best (3 of 4), among the PLIP-profiled hits:")
    for x in rows:
        if x["n_pass"] == 3 and x["profiled"]:
            fails = [k for k in ("c1_score", "c2_hinge", "c3_LE", "c4_cns") if not x[k]]
            print(f"  rank {x['rank']:3d}  {x['name'][:26]:26s} {x['score']:6.2f}  "
                  f"LE {x['LE']:.3f}  MPO {x['cns_mpo4']:.1f}  hinge {x['hinge_hbond']:8s}"
                  f"  fails {fails}  {x['alerts'][:34]}")
    print("\nwrote", os.path.join(ANA, "shortlist.csv"))


if __name__ == "__main__":
    main()

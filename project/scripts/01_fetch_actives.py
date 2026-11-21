"""
Stage 1 / Stage 4 : known GSK-3beta actives (positive controls) from ChEMBL.

Step 1  pull every binding-assay activity recorded against human GSK-3beta
        (target CHEMBL262) that carries a pChEMBL value, i.e. a curated
        IC50/Ki/Kd on a -log10(M) scale.
Step 2  aggregate per molecule (median pChEMBL, number of independent
        measurements).
Step 3  keep the named tool compounds / clinical candidates we want as
        positive controls, each identified by its ChEMBL synonym.

Outputs
    04_ligands/actives/gsk3b_chembl_activities.csv   full aggregated table
    04_ligands/actives/known_actives.csv             the curated control set
"""
import json, time, urllib.request, urllib.parse, statistics, csv, os, collections

BASE = "https://www.ebi.ac.uk/chembl/api/data"
TARGET = "CHEMBL262"
HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "..", "04_ligands", "actives")

# Positive controls we want, matched case-insensitively against ChEMBL
# synonyms / preferred names. Mechanism notes come from the primary papers.
CONTROLS = {
    "AR-A014418":            "ATP-competitive thiazole urea; co-crystallised with GSK-3b in PDB 1Q5K (Bhat et al. JBC 2003)",
    "SB-216763":             "ATP-competitive maleimide (Coghlan et al. Chem Biol 2000)",
    "SB-415286":             "ATP-competitive maleimide (Coghlan et al. Chem Biol 2000)",
    "LADUVIGLUSIB":          "CHIR-99021, highly selective ATP-competitive aminopyrimidine",
    "TIDEGLUSIB":            "non-ATP-competitive thiadiazolidinedione; phase II in Alzheimer's (ARGO)",
    "LY-2090314":            "ATP-competitive maleimide, clinical candidate",
    "ALSTERPAULLONE":        "ATP-competitive paullone",
    "KENPAULLONE":           "ATP-competitive paullone",
    "INDIRUBIN-3'-MONOXIME": "ATP-competitive bisindole",
    "TDZD-8":                "non-ATP-competitive thiadiazolidinedione tool compound",
    "STAUROSPORINE":         "pan-kinase reference inhibitor",
}


def get(url, tries=4):
    req = urllib.request.Request(url, headers={"Accept": "application/json",
                                               "User-Agent": "cadd-miniproject/1.0"})
    for a in range(tries):
        try:
            with urllib.request.urlopen(req, timeout=120) as r:
                return json.load(r)
        except Exception:
            if a == tries - 1:
                raise
            time.sleep(3 * (a + 1))


def fetch_activities():
    rows, offset, limit = [], 0, 1000
    while True:
        url = (f"{BASE}/activity.json?target_chembl_id={TARGET}"
               f"&pchembl_value__isnull=false&assay_type=B"
               f"&limit={limit}&offset={offset}")
        d = get(url)
        acts = d.get("activities", [])
        rows.extend(acts)
        print(f"  fetched {len(rows)} activities", end="\r")
        if len(acts) < limit or offset > 20000:
            break
        offset += limit
        time.sleep(0.2)
    print()
    return rows


def main():
    os.makedirs(OUT, exist_ok=True)
    print("Fetching GSK-3beta (CHEMBL262) binding activities from ChEMBL ...")
    acts = fetch_activities()

    per_mol = collections.defaultdict(list)
    types = collections.defaultdict(collections.Counter)
    for a in acts:
        try:
            v = float(a["pchembl_value"])
        except (TypeError, ValueError, KeyError):
            continue
        mid = a["molecule_chembl_id"]
        per_mol[mid].append(v)
        types[mid][a.get("standard_type", "?")] += 1

    print(f"{len(acts)} activities on {len(per_mol)} distinct molecules")

    agg = []
    for mid, vals in per_mol.items():
        agg.append({"chembl_id": mid,
                    "pchembl_median": round(statistics.median(vals), 2),
                    "pchembl_max": round(max(vals), 2),
                    "n": len(vals),
                    "main_type": types[mid].most_common(1)[0][0]})
    agg.sort(key=lambda r: (-r["pchembl_median"], -r["n"]))
    with open(os.path.join(OUT, "gsk3b_chembl_activities.csv"), "w", newline="", encoding="utf8") as fh:
        w = csv.DictWriter(fh, fieldnames=list(agg[0].keys())); w.writeheader(); w.writerows(agg)

    # resolve the named controls
    print("\nResolving positive controls by synonym ...")
    wanted = {k.upper(): v for k, v in CONTROLS.items()}
    found = {}
    for name in list(wanted):
        q = urllib.parse.quote(name)
        for url in (f"{BASE}/molecule.json?molecule_synonyms__molecule_synonym__iexact={q}&limit=5",
                    f"{BASE}/molecule.json?pref_name__iexact={q}&limit=5",
                    f"{BASE}/molecule/search.json?q={q}&limit=5"):
            try:
                d = get(url, tries=2)
            except Exception:
                continue
            mols = d.get("molecules", [])
            if mols:
                m = mols[0]
                st = m.get("molecule_structures") or {}
                if st.get("canonical_smiles"):
                    found[name] = m
                    break
            time.sleep(0.2)
        if name not in found:
            print(f"  !! could not resolve {name}")

    rows = []
    for name, m in found.items():
        mid = m["molecule_chembl_id"]
        st = m["molecule_structures"]; pr = m.get("molecule_properties") or {}
        rec = next((r for r in agg if r["chembl_id"] == mid), None)
        pch = rec["pchembl_median"] if rec else None
        rows.append({
            "name": name.title() if not name.startswith(("AR-", "SB-", "LY-", "TDZD", "BIO")) else name,
            "chembl_id": mid,
            "pref_name": m.get("pref_name") or "",
            "smiles": st["canonical_smiles"],
            "pchembl_median": pch if pch else "",
            "potency_nM": round(10 ** (9 - pch), 1) if pch else "",
            "n_measurements": rec["n"] if rec else 0,
            "assay_type": rec["main_type"] if rec else "",
            "mw": pr.get("full_mwt"), "alogp": pr.get("alogp"), "hbd": pr.get("hbd"),
            "hba": pr.get("hba"), "tpsa": pr.get("psa"), "rotb": pr.get("rtb"),
            "max_phase": m.get("max_phase"),
            "mechanism_note": CONTROLS[name],
        })
        print(f"  {name:24s} {mid:15s} pChEMBL {pch}  n={rec['n'] if rec else 0}")

    rows.sort(key=lambda r: -(r["pchembl_median"] or 0))
    path = os.path.join(OUT, "known_actives.csv")
    with open(path, "w", newline="", encoding="utf8") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys())); w.writeheader(); w.writerows(rows)
    print("\nwrote", path, len(rows), "positive controls")


if __name__ == "__main__":
    main()

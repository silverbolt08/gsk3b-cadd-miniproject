"""
Stage 4 - build the screening library.

Hypothesis: a GSK-3beta inhibitor for Alzheimer's has to reach the brain, and
the fastest route to that is a molecule that is already in humans. The library
is therefore two hypothesis-driven blocks:

  Block A (80)  approved small-molecule drugs (ChEMBL max_phase = 4) filtered
                with CNS-MPO-like physicochemical criteria, then reduced to 80
                by MaxMin diversity picking on Morgan fingerprints so the block
                spans chemotype space instead of 80 near-duplicates.

  Block B (~40) natural products with published activity in Alzheimer models or
                on GSK-3 / tau, pulled from ChEMBL by name. These are the
                "nutraceutical" claims the AD literature is full of - docking
                them against the same pocket is a fair test of those claims.

Positive controls (11 known GSK-3beta actives) come from 01_fetch_actives.py
and are docked alongside the library.

Outputs: 04_ligands/library/screening_library.csv
         04_ligands/library/library_composition.csv
"""
import json, os, time, urllib.request, urllib.parse, csv, random

BASE = "https://www.ebi.ac.uk/chembl/api/data"
HERE = os.path.dirname(os.path.abspath(__file__))
LIB = os.path.join(HERE, "..", "04_ligands", "library")
ACT = os.path.join(HERE, "..", "04_ligands", "actives")
N_DRUGS = 80
SEED = 2024

# CNS-MPO-like filter (Wager et al. ACS Chem Neurosci 2010, simplified)
CNS_FILTER = dict(mw_max=450.0, mw_min=150.0, alogp_min=0.5, alogp_max=5.0,
                  hbd_max=3, tpsa_max=90.0, rotb_max=8)

NATURAL_PRODUCTS = [
    "CURCUMIN", "RESVERATROL", "QUERCETIN", "LUTEOLIN", "APIGENIN", "FISETIN",
    "GENISTEIN", "BERBERINE", "HUPERZINE A", "GALANTAMINE", "MYRICETIN",
    "KAEMPFEROL", "NARINGENIN", "HESPERETIN", "BAICALEIN", "CHRYSIN",
    "CATECHIN", "EPICATECHIN", "SILIBININ", "ROSMARINIC ACID", "FERULIC ACID",
    "CAFFEIC ACID", "TANSHINONE IIA", "CELASTROL", "URSOLIC ACID",
    "BETULINIC ACID", "ASIATIC ACID", "PIPERINE", "SULFORAPHANE",
    "THYMOQUINONE", "HARMINE", "HONOKIOL", "MAGNOLOL", "EMODIN",
    "OLEUROPEIN", "HYDROXYTYROSOL", "CROCETIN", "GINKGOLIDE B", "BILOBALIDE",
    "WITHAFERIN A", "SCOPOLETIN", "ICARIIN", "PAEONIFLORIN", "GASTRODIN",
    "VANILLIN", "PUERARIN", "DAIDZEIN", "NOBILETIN", "TANGERETIN", "LIMONIN",
]


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


def fetch_approved():
    """All approved small molecules that pass the CNS-like property filter."""
    f = CNS_FILTER
    q = (f"{BASE}/molecule.json?max_phase=4&molecule_type=Small%20molecule"
         f"&molecule_properties__full_mwt__lte={f['mw_max']}"
         f"&molecule_properties__full_mwt__gte={f['mw_min']}"
         f"&molecule_properties__hbd__lte={f['hbd_max']}"
         f"&molecule_properties__psa__lte={f['tpsa_max']}"
         f"&molecule_properties__rtb__lte={f['rotb_max']}"
         f"&limit=1000")
    out, url = [], q + "&offset=0"
    offset = 0
    while True:
        d = get(q + f"&offset={offset}")
        mols = d.get("molecules", [])
        out.extend(mols)
        print(f"  approved+filtered: {len(out)}", end="\r")
        if len(mols) < 1000 or offset > 6000:
            break
        offset += 1000
        time.sleep(0.2)
    print()
    keep = []
    for m in out:
        p = m.get("molecule_properties") or {}
        s = (m.get("molecule_structures") or {}).get("canonical_smiles")
        try:
            alogp = float(p.get("alogp"))
        except (TypeError, ValueError):
            continue
        if not s or "." in s:                      # skip salts / mixtures
            continue
        if not (f["alogp_min"] <= alogp <= f["alogp_max"]):
            continue
        # the ATP site reads a flat heteroaromatic that can hydrogen-bond to the
        # hinge, so require at least one aromatic ring and one acceptor
        if (p.get("aromatic_rings") or 0) < 1 or (p.get("hba") or 0) < 1:
            continue
        keep.append({"chembl_id": m["molecule_chembl_id"],
                     "name": (m.get("pref_name") or m["molecule_chembl_id"]).title(),
                     "smiles": s, "mw": p.get("full_mwt"), "alogp": alogp,
                     "hbd": p.get("hbd"), "hba": p.get("hba"),
                     "tpsa": p.get("psa"), "rotb": p.get("rtb"),
                     "block": "A_approved_drug"})
    print(f"  {len(keep)} approved drugs pass the CNS-like filter")
    return keep


def diversity_pick(cands, n, seed=SEED):
    from rdkit import Chem
    from rdkit.Chem import rdFingerprintGenerator
    from rdkit.SimDivFilters.rdSimDivPickers import MaxMinPicker
    gen = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)
    fps, ok = [], []
    for c in cands:
        m = Chem.MolFromSmiles(c["smiles"])
        if m is None:
            continue
        fps.append(gen.GetFingerprint(m)); ok.append(c)
    from rdkit import DataStructs

    def dist(i, j, fps=fps):
        return 1.0 - DataStructs.TanimotoSimilarity(fps[i], fps[j])

    picker = MaxMinPicker()
    idx = picker.LazyPick(dist, len(ok), n, seed=seed)
    return [ok[i] for i in idx]


def fetch_natural_products():
    rows = []
    for name in NATURAL_PRODUCTS:
        q = urllib.parse.quote(name)
        m = None
        for url in (f"{BASE}/molecule.json?pref_name__iexact={q}&limit=3",
                    f"{BASE}/molecule.json?molecule_synonyms__molecule_synonym__iexact={q}&limit=3"):
            try:
                d = get(url, tries=2)
            except Exception:
                continue
            for cand in d.get("molecules", []):
                s = (cand.get("molecule_structures") or {}).get("canonical_smiles")
                if s and "." not in s:
                    m = cand
                    break
            if m:
                break
            time.sleep(0.2)
        if not m:
            print(f"  !! {name} not resolved")
            continue
        p = m.get("molecule_properties") or {}
        rows.append({"chembl_id": m["molecule_chembl_id"], "name": name.title(),
                     "smiles": m["molecule_structures"]["canonical_smiles"],
                     "mw": p.get("full_mwt"), "alogp": p.get("alogp"),
                     "hbd": p.get("hbd"), "hba": p.get("hba"),
                     "tpsa": p.get("psa"), "rotb": p.get("rtb"),
                     "block": "B_natural_product"})
    print(f"  {len(rows)} natural products resolved")
    return rows


def main():
    os.makedirs(LIB, exist_ok=True)
    random.seed(SEED)
    print("Block A: approved drugs from ChEMBL ...")
    approved = fetch_approved()
    picked = diversity_pick(approved, N_DRUGS)
    print(f"  diversity-picked {len(picked)} of {len(approved)}")

    print("Block B: natural products ...")
    nps = fetch_natural_products()

    # positive controls, so we can drop any accidental duplicate
    ctrl_ids = set()
    ap = os.path.join(ACT, "known_actives.csv")
    if os.path.exists(ap):
        ctrl_ids = {r["chembl_id"] for r in csv.DictReader(open(ap, encoding="utf8"))}

    seen, rows = set(), []
    for r in picked + nps:
        if r["chembl_id"] in ctrl_ids or r["chembl_id"] in seen:
            continue
        seen.add(r["chembl_id"])
        rows.append(r)

    # stable, readable ligand ids
    for i, r in enumerate(rows, 1):
        r["ligand_id"] = f"LIG{i:03d}"

    cols = ["ligand_id", "name", "chembl_id", "block", "smiles", "mw", "alogp",
            "hbd", "hba", "tpsa", "rotb"]
    with open(os.path.join(LIB, "screening_library.csv"), "w", newline="", encoding="utf8") as fh:
        w = csv.DictWriter(fh, fieldnames=cols); w.writeheader()
        for r in rows:
            w.writerow({c: r.get(c, "") for c in cols})
    print(f"\nscreening library: {len(rows)} compounds -> "
          f"{os.path.join(LIB, 'screening_library.csv')}")
    from collections import Counter
    print(Counter(r["block"] for r in rows))


if __name__ == "__main__":
    main()

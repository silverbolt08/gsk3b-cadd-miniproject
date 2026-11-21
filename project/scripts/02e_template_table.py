"""
Stage 2 - summarise the template search.

Two searches were run against the PDB with the full-length P49841 sequence:
  * NCBI BLASTp (blastp, database pdb, E < 1e-5) - the classic route, run
    through the NCBI URL API; the raw report, when the NCBI queue returns it,
    is in 02_structure/alignment/blast_pdb_report.txt
  * the RCSB sequence search service (MMseqs2 profile search), which returns
    the same ranking immediately and is used here as the primary record:
    02_structure/alignment/rcsb_sequence_search.json

This script turns the hit list into the template table used in the report:
identity, coverage, method, resolution, bound ligand, and the modelled range,
and flags the three templates taken forward.

Output: 02_structure/alignment/template_table.csv
"""
import json, os, urllib.request, csv, time

HERE = os.path.dirname(os.path.abspath(__file__))
ALN = os.path.join(HERE, "..", "02_structure", "alignment")
CHOSEN = {"1Q5K", "1J1B", "4ACC"}
N_SHOW = 14


def get(url):
    req = urllib.request.Request(url, headers={"User-Agent": "cadd-miniproject/1.0"})
    with urllib.request.urlopen(req, timeout=60) as r:
        return json.load(r)


def main():
    hits = json.load(open(os.path.join(ALN, "rcsb_sequence_search.json")))
    # one row per PDB entry, best entity only
    seen, rows = set(), []
    for h in hits:
        pdb = h["entity"].split("_")[0]
        if pdb in seen:
            continue
        seen.add(pdb)
        rows.append(h)
        if len(rows) >= 60:
            break

    out = []
    for h in rows:
        pdb = h["entity"].split("_")[0]
        try:
            e = get(f"https://data.rcsb.org/rest/v1/core/entry/{pdb}")
        except Exception:
            continue
        method = (e.get("exptl") or [{}])[0].get("method", "?")
        res = (e.get("rcsb_entry_info") or {}).get("resolution_combined", [None])[0]
        title = (e.get("struct") or {}).get("title", "")[:70]
        ligs = (e.get("rcsb_entry_info") or {}).get("nonpolymer_bound_components") or []
        out.append({"pdb_id": pdb, "identity_pct": round(100 * h["identity"], 1),
                    "evalue": h["evalue"], "aligned_length": h["aln_len"],
                    "method": method, "resolution_A": res,
                    "ligands": ";".join(ligs[:4]), "title": title,
                    "used_as_template": "yes" if pdb in CHOSEN else ""})
        time.sleep(0.1)
        if len(out) >= N_SHOW and all(p in seen for p in CHOSEN):
            if all(any(r["pdb_id"] == p for r in out) for p in CHOSEN):
                break

    # make sure the three chosen templates are in the table
    for p in CHOSEN:
        if not any(r["pdb_id"] == p for r in out):
            e = get(f"https://data.rcsb.org/rest/v1/core/entry/{p}")
            h = next((x for x in hits if x["entity"].split("_")[0] == p), None)
            out.append({"pdb_id": p,
                        "identity_pct": round(100 * h["identity"], 1) if h else "",
                        "evalue": h["evalue"] if h else "", "aligned_length": h["aln_len"] if h else "",
                        "method": (e.get("exptl") or [{}])[0].get("method", "?"),
                        "resolution_A": (e.get("rcsb_entry_info") or {}).get("resolution_combined", [None])[0],
                        "ligands": ";".join(((e.get("rcsb_entry_info") or {}).get("nonpolymer_bound_components") or [])[:4]),
                        "title": (e.get("struct") or {}).get("title", "")[:70],
                        "used_as_template": "yes"})

    out.sort(key=lambda r: (-(r["identity_pct"] or 0), r["resolution_A"] or 99))
    path = os.path.join(ALN, "template_table.csv")
    with open(path, "w", newline="", encoding="utf8") as fh:
        w = csv.DictWriter(fh, fieldnames=list(out[0].keys())); w.writeheader(); w.writerows(out)
    print(f"wrote {path} ({len(out)} entries)")
    for r in out[:N_SHOW]:
        print(f"  {r['pdb_id']}  id {r['identity_pct']:5}%  {str(r['resolution_A']):>5} A  "
              f"{r['ligands'][:18]:18s} {r['used_as_template']:3s} {r['title'][:46]}")


if __name__ == "__main__":
    main()

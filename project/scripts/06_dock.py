"""
Stage 5 - run the virtual screen with AutoDock Vina 1.2.5.

Every ligand is docked into the rigid receptor with the box from
03_pocket/box.json. One Vina process per ligand, N processes in parallel,
each with --cpu 1 so the run is reproducible given the fixed --seed.

Usage:
    python 06_dock.py --receptor model     # production screen (homology model)
    python 06_dock.py --receptor crystal   # cross-check on 1Q5K chain A
    python 06_dock.py --receptor model --only ACT10   # single ligand

Outputs: 05_docking/<tag>/poses/<id>.pdbqt
         05_docking/<tag>/logs/<id>.log
         05_docking/<tag>/config_<id>.txt   (one per ligand, same box)
         05_docking/<tag>/scores.csv
"""
import argparse, csv, json, os, subprocess, sys, time
from concurrent.futures import ProcessPoolExecutor, as_completed

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, ".."))
VINA = r"C:\CADD\tools\vina\vina_1.2.5_win.exe"
LIGS = os.path.join(ROOT, "04_ligands", "pdbqt")
BOX = json.load(open(os.path.join(ROOT, "03_pocket", "box.json")))
RECEPTORS = {
    "model":   os.path.join(ROOT, "05_docking", "receptor", "gsk3b_model.pdbqt"),
    "crystal": os.path.join(ROOT, "05_docking", "receptor", "1q5kA.pdbqt"),
}
EXHAUSTIVENESS = 32
NUM_MODES = 12
ENERGY_RANGE = 4
SEED = 20241121


def write_config(path, receptor, ligand, out):
    c, s = BOX["center"], BOX["size"]
    txt = (f"receptor = {receptor}\n"
           f"ligand = {ligand}\n"
           f"out = {out}\n\n"
           f"center_x = {c[0]}\ncenter_y = {c[1]}\ncenter_z = {c[2]}\n"
           f"size_x = {s[0]}\nsize_y = {s[1]}\nsize_z = {s[2]}\n\n"
           f"exhaustiveness = {EXHAUSTIVENESS}\n"
           f"num_modes = {NUM_MODES}\n"
           f"energy_range = {ENERGY_RANGE}\n"
           f"seed = {SEED}\ncpu = 1\n")
    open(path, "w").write(txt)
    return path


def run_one(args):
    lig_id, receptor, outdir = args
    lig = os.path.join(LIGS, lig_id + ".pdbqt")
    pose = os.path.join(outdir, "poses", lig_id + ".pdbqt")
    logf = os.path.join(outdir, "logs", lig_id + ".log")
    cfg = write_config(os.path.join(outdir, "configs", f"config_{lig_id}.txt"),
                       receptor, lig, pose)
    t0 = time.time()
    try:
        p = subprocess.run([VINA, "--config", cfg], capture_output=True,
                           text=True, timeout=1800)
        open(logf, "w").write(p.stdout + "\n" + p.stderr)
        if p.returncode != 0:
            return lig_id, None, p.stderr.strip()[:200], time.time() - t0
    except subprocess.TimeoutExpired:
        return lig_id, None, "timeout", time.time() - t0
    best, modes = None, []
    for line in open(pose):
        if line.startswith("REMARK VINA RESULT"):
            f = line.split()
            modes.append(float(f[3]))
    if modes:
        best = modes[0]
    return lig_id, best, modes, time.time() - t0


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--receptor", choices=list(RECEPTORS), default="model")
    ap.add_argument("--only", default=None, help="comma-separated ligand ids")
    ap.add_argument("--workers", type=int, default=14)
    a = ap.parse_args()

    receptor = RECEPTORS[a.receptor]
    outdir = os.path.join(ROOT, "05_docking", a.receptor)
    for sub in ("poses", "logs", "configs"):
        os.makedirs(os.path.join(outdir, sub), exist_ok=True)

    man = list(csv.DictReader(open(os.path.join(ROOT, "04_ligands", "ligand_manifest.csv"),
                                   encoding="utf8")))
    ids = [r["id"] for r in man]
    if a.only:
        want = set(a.only.split(","))
        ids = [i for i in ids if i in want]
    meta = {r["id"]: r for r in man}

    print(f"receptor : {receptor}")
    print(f"box      : centre {BOX['center']} size {BOX['size']}")
    print(f"ligands  : {len(ids)}   workers {a.workers}   "
          f"exhaustiveness {EXHAUSTIVENESS}  seed {SEED}\n")

    rows, t0 = [], time.time()
    with ProcessPoolExecutor(max_workers=a.workers) as ex:
        futs = {ex.submit(run_one, (i, receptor, outdir)): i for i in ids}
        for n, fut in enumerate(as_completed(futs), 1):
            lig_id, best, modes, dt = fut.result()
            m = meta[lig_id]
            if best is None:
                print(f"[{n:3d}/{len(ids)}] {lig_id} FAILED: {modes}")
                continue
            rows.append({"id": lig_id, "name": m["name"], "set": m["set"],
                         "chembl_id": m["chembl_id"], "potency_nM": m["potency_nM"],
                         "best_score": best,
                         "mode_scores": ";".join(f"{x:.2f}" for x in modes),
                         "n_modes": len(modes), "mw": m["mw"], "clogp": m["clogp"],
                         "tpsa": m["tpsa"], "hbd": m["hbd"], "hba": m["hba"],
                         "rotb": m["rotb"], "heavy_atoms": m["heavy"],
                         "ligand_efficiency": round(-best / int(m["heavy"]), 3),
                         "seconds": round(dt, 1)})
            print(f"[{n:3d}/{len(ids)}] {lig_id:7s} {m['name'][:26]:26s} "
                  f"{best:7.2f} kcal/mol  ({dt:.0f}s)")

    rows.sort(key=lambda r: r["best_score"])
    out = os.path.join(outdir, "scores.csv")
    with open(out, "w", newline="", encoding="utf8") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys())); w.writeheader(); w.writerows(rows)
    print(f"\n{len(rows)} ligands docked in {(time.time()-t0)/60:.1f} min -> {out}")


if __name__ == "__main__":
    main()

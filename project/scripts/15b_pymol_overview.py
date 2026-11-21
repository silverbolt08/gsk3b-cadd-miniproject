"""
Render the two remaining structural figures in PyMOL:

  fig1_model_overview  the comparative model, its lobes, the rebuilt loops,
                       the P2Rank pocket-1 residues, the crystallographic
                       AR-A014418 position and the 26 A search box (drawn as
                       a CGO wireframe)
  fig9_redock          the redocking test: crystal pose (yellow) against
                       Vina's top pose (cyan) and the best of the 12 modes
                       (magenta)

Both are ray-traced headless. The .pml scripts are written to
06_analysis/figures/pymol/ so they can be reopened in the PyMOL GUI.
"""
import json, os, subprocess, sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, ".."))
FIG = os.path.join(ROOT, "06_analysis", "figures")
PML = os.path.join(FIG, "pymol")
MODEL = os.path.join(ROOT, "02_structure", "models", "gsk3b_model_aligned.pdb")
REF = os.path.join(ROOT, "05_docking", "ref_ligand_TMU.pdb")
POSE = os.path.join(ROOT, "05_docking", "crystal", "poses", "ACT10.pdbqt")
REC = os.path.join(ROOT, "05_docking", "receptor", "1q5kA_H.pdb")
BOX = json.load(open(os.path.join(ROOT, "03_pocket", "box.json")))
RED = json.load(open(os.path.join(ROOT, "05_docking", "redock_validation.json")))
OFFSET = 34

PRE = """reinitialize
bg_color white
set ray_opaque_background, 1
set antialias, 2
set ray_shadows, 0
set specular, 0.15
set ambient, 0.30
set direct, 0.55
set cartoon_fancy_helices, 1
set label_size, 16
set label_color, black
set label_font_id, 7
set dash_gap, 0.30
set dash_width, 2.2
"""


def p(path):
    return path.replace("\\", "/")


def cgo_box(centre, size, color=(0.89, 0.71, 0.0), name="searchbox"):
    """wireframe box as a CGO line set"""
    cx, cy, cz = centre
    hx, hy, hz = [s / 2.0 for s in size]
    v = [(cx - hx, cy - hy, cz - hz), (cx + hx, cy - hy, cz - hz),
         (cx + hx, cy + hy, cz - hz), (cx - hx, cy + hy, cz - hz),
         (cx - hx, cy - hy, cz + hz), (cx + hx, cy - hy, cz + hz),
         (cx + hx, cy + hy, cz + hz), (cx - hx, cy + hy, cz + hz)]
    edges = [(0, 1), (1, 2), (2, 3), (3, 0), (4, 5), (5, 6), (6, 7), (7, 4),
             (0, 4), (1, 5), (2, 6), (3, 7)]
    out = ["python",
           "from pymol.cgo import BEGIN, LINES, COLOR, VERTEX, END, LINEWIDTH",
           "from pymol import cmd",
           f"obj = [LINEWIDTH, 2.0, BEGIN, LINES, COLOR, {color[0]}, {color[1]}, {color[2]},"]
    for a, b in edges:
        out.append(f"         VERTEX, {v[a][0]:.3f}, {v[a][1]:.3f}, {v[a][2]:.3f},"
                   f" VERTEX, {v[b][0]:.3f}, {v[b][1]:.3f}, {v[b][2]:.3f},")
    out.append("         END]")
    out.append(f"cmd.load_cgo(obj, '{name}')")
    out.append("python end")
    return out


def overview():
    pocket = BOX["pocket_residues_model"]
    resi = "+".join(str(n) for n in pocket)
    loops = "86+87+253+254+255"          # rebuilt loops in model numbering
    lines = [PRE,
             f"load {p(MODEL)}, gsk3b",
             f"load {p(REF)}, tmu",
             "hide everything",
             "show cartoon, gsk3b",
             "color grey75, gsk3b",
             "color skyblue, gsk3b and resi 1-103",          # N-lobe (35-137)
             "color firebrick, gsk3b and resi " + loops,
             f"show sticks, gsk3b and resi {resi} and not (name C+N+O and not resn PRO)",
             "color palegreen, gsk3b and resi " + resi + " and elem C",
             "show sticks, tmu",
             "color cyan, tmu and elem C",
             "set stick_radius, 0.20, tmu",
             "set cartoon_transparency, 0.12, gsk3b",
             "util.cnc('tmu or (gsk3b and resi " + resi + ")')",
             ]
    lines += cgo_box(BOX["center"], BOX["size"])
    lines += ["orient gsk3b",
              "turn z, 180",
              "turn y, 10",
              "zoom gsk3b, 3",
              f"png {p(os.path.join(FIG, 'fig1_model_overview.png'))}, "
              "width=1900, height=1400, dpi=300, ray=1"]
    path = os.path.join(PML, "fig1_model_overview.pml")
    open(path, "w").write("\n".join(lines) + "\n")
    return path



def write_pose_pdb(conf_index, out_path, resname="LIG"):
    """one Vina mode from the PDBQT as a plain PDB, so PyMOL loads it cleanly"""
    from meeko import PDBQTMolecule, RDKitMolCreate
    from rdkit import Chem, RDLogger
    RDLogger.DisableLog("rdApp.*")
    pm = PDBQTMolecule.from_file(POSE, skip_typing=True)
    mol = Chem.RemoveHs(RDKitMolCreate.from_pdbqt_mol(pm, only_cluster_leads=False)[0])
    conf = mol.GetConformer(conf_index)
    with open(out_path, "w") as fh:
        for i, a in enumerate(mol.GetAtoms(), 1):
            q = conf.GetAtomPosition(a.GetIdx())
            el = a.GetSymbol()
            fh.write(f"HETATM{i:5d} {(el + str(i)):<4s} {resname:<3s} Z   1    "
                     f"{q.x:8.3f}{q.y:8.3f}{q.z:8.3f}  1.00  0.00          {el:>2s}\n")
        fh.write("END\n")
    return out_path


def redock():
    best = RED["receptors"]["crystal"]["best_rmsd_mode"]
    tmp = os.path.join(PML, "_pose_top.pdb")
    tmp2 = os.path.join(PML, "_pose_best.pdb")
    os.makedirs(PML, exist_ok=True)
    top_pdb = write_pose_pdb(0, tmp)
    best_pdb = write_pose_pdb(best - 1, tmp2)
    show = [62, 70, 83, 85, 132, 133, 134, 135, 138, 188, 199, 200]
    resi = "+".join(str(n) for n in show)
    lines = [PRE,
             f"load {p(REC)}, rec",
             f"load {p(REF)}, xtal",
             f"load {p(top_pdb)}, posetop",
             f"load {p(best_pdb)}, posebest",
             "hide everything",
             "show cartoon, rec",
             "color grey80, rec",
             "set cartoon_transparency, 0.55, rec",
             f"show sticks, rec and resi {resi} and not (name C+N+O and not resn PRO)",
             "color grey55, rec and elem C",
             "show sticks, posetop",
             "show sticks, posebest",
             "color cyan, posetop and elem C",
             "color magenta, posebest and elem C",
             "show sticks, xtal",
             "color yellow, xtal and elem C",
             "set stick_radius, 0.20, xtal",
             "set stick_radius, 0.15, posetop",
             "set stick_radius, 0.15, posebest",
             "hide everything, hydro",
             "util.cnc('xtal or posetop or posebest')",
             'label rec and resi 99 and name CB, "Asp133"',
             'label rec and resi 101 and name CB, "Val135"',
             "orient xtal",
             "zoom xtal, 6.5",
             f"png {p(os.path.join(FIG, 'fig9_redock.png'))}, "
             "width=1800, height=1200, dpi=300, ray=1"]
    path = os.path.join(PML, "fig9_redock.pml")
    open(path, "w").write("\n".join(lines) + "\n")
    return path


def run(pml):
    exe = os.path.join(os.path.expanduser("~"), "AppData", "Roaming", "Python",
                       f"Python{sys.version_info.major}{sys.version_info.minor}",
                       "Scripts", "pymol.exe")
    cmd = [exe, "-cq", pml] if os.path.exists(exe) else [sys.executable, "-m", "pymol", "-cq", pml]
    r = subprocess.run(cmd, capture_output=True, text=True, timeout=900)
    return r.returncode, (r.stdout + r.stderr)[-600:]


if __name__ == "__main__":
    os.makedirs(PML, exist_ok=True)
    for maker, fig in ((overview, "fig1_model_overview"), (redock, "fig9_redock")):
        pml = maker()
        rc, out = run(pml)
        png = os.path.join(FIG, fig + ".png")
        size = os.path.getsize(png) if os.path.exists(png) else 0
        print(f"  {fig}.png  rc={rc}  {size/1e3:.0f} kB")
        if not size:
            print("   ", out)

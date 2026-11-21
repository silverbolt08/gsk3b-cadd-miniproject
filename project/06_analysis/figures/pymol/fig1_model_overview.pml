reinitialize
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

load C:/CADD/project/02_structure/models/gsk3b_model_aligned.pdb, gsk3b
load C:/CADD/project/05_docking/ref_ligand_TMU.pdb, tmu
hide everything
show cartoon, gsk3b
color grey75, gsk3b
color skyblue, gsk3b and resi 1-103
color firebrick, gsk3b and resi 86+87+253+254+255
show sticks, gsk3b and resi 28+29+30+31+33+36+49+51+63+76+98+99+101+102+103+104+107+151+152+154+165+166 and not (name C+N+O and not resn PRO)
color palegreen, gsk3b and resi 28+29+30+31+33+36+49+51+63+76+98+99+101+102+103+104+107+151+152+154+165+166 and elem C
show sticks, tmu
color cyan, tmu and elem C
set stick_radius, 0.20, tmu
set cartoon_transparency, 0.12, gsk3b
util.cnc('tmu or (gsk3b and resi 28+29+30+31+33+36+49+51+63+76+98+99+101+102+103+104+107+151+152+154+165+166)')
python
from pymol.cgo import BEGIN, LINES, COLOR, VERTEX, END, LINEWIDTH
from pymol import cmd
obj = [LINEWIDTH, 2.0, BEGIN, LINES, COLOR, 0.89, 0.71, 0.0,
         VERTEX, 9.328, 11.372, -3.841, VERTEX, 35.328, 11.372, -3.841,
         VERTEX, 35.328, 11.372, -3.841, VERTEX, 35.328, 37.372, -3.841,
         VERTEX, 35.328, 37.372, -3.841, VERTEX, 9.328, 37.372, -3.841,
         VERTEX, 9.328, 37.372, -3.841, VERTEX, 9.328, 11.372, -3.841,
         VERTEX, 9.328, 11.372, 22.159, VERTEX, 35.328, 11.372, 22.159,
         VERTEX, 35.328, 11.372, 22.159, VERTEX, 35.328, 37.372, 22.159,
         VERTEX, 35.328, 37.372, 22.159, VERTEX, 9.328, 37.372, 22.159,
         VERTEX, 9.328, 37.372, 22.159, VERTEX, 9.328, 11.372, 22.159,
         VERTEX, 9.328, 11.372, -3.841, VERTEX, 9.328, 11.372, 22.159,
         VERTEX, 35.328, 11.372, -3.841, VERTEX, 35.328, 11.372, 22.159,
         VERTEX, 35.328, 37.372, -3.841, VERTEX, 35.328, 37.372, 22.159,
         VERTEX, 9.328, 37.372, -3.841, VERTEX, 9.328, 37.372, 22.159,
         END]
cmd.load_cgo(obj, 'searchbox')
python end
orient gsk3b
turn z, 180
turn y, 10
zoom gsk3b, 3
png C:/CADD/project/06_analysis/figures/fig1_model_overview.png, width=1900, height=1400, dpi=300, ray=1

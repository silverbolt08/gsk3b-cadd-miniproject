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

load C:/CADD/project/05_docking/receptor/1q5kA_H.pdb, rec
load C:/CADD/project/05_docking/ref_ligand_TMU.pdb, xtal
load C:/CADD/project/06_analysis/figures/pymol/_pose_top.pdb, posetop
load C:/CADD/project/06_analysis/figures/pymol/_pose_best.pdb, posebest
hide everything
show cartoon, rec
color grey80, rec
set cartoon_transparency, 0.55, rec
show sticks, rec and resi 62+70+83+85+132+133+134+135+138+188+199+200 and not (name C+N+O and not resn PRO)
color grey55, rec and elem C
show sticks, posetop
show sticks, posebest
color cyan, posetop and elem C
color magenta, posebest and elem C
show sticks, xtal
color yellow, xtal and elem C
set stick_radius, 0.20, xtal
set stick_radius, 0.15, posetop
set stick_radius, 0.15, posebest
hide everything, hydro
util.cnc('xtal or posetop or posebest')
label rec and resi 99 and name CB, "Asp133"
label rec and resi 101 and name CB, "Val135"
orient xtal
zoom xtal, 6.5
png C:/CADD/project/06_analysis/figures/fig9_redock.png, width=1800, height=1200, dpi=300, ray=1

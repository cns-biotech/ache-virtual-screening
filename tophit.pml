
load receptor_prep.pdb, rec
load poses/CHEMBL4214707.sdf, hit
hide everything
bg_color white
show surface, rec
color grey70, rec
set two_sided_lighting, on
set ray_interior_color, grey50
set transparency, 0.15
show sticks, hit and not hydro
color orange, hit and elem C
set ray_shadows, 0
set antialias, 2
orient hit
zoom hit, 6
png fig_tophit.png, width=1600, height=1200, dpi=300, ray=1


load receptor_prep.pdb, rec
load ligand_ref.pdb, cryst
load redock_poses.sdf, redock
set state, 1
hide everything
bg_color white
show cartoon, rec
color grey80, rec
set cartoon_transparency, 0.65
select pocket, byres (polymer within 5 of cryst)
show sticks, pocket and sidechain and not hydro
color palecyan, pocket and elem C
show sticks, cryst and not hydro
show sticks, redock and not hydro
color yellow, cryst and elem C
color magenta, redock and elem C
set ray_shadows, 0
set antialias, 2
orient cryst
zoom cryst, 4
png fig_redock.png, width=1600, height=1200, dpi=300, ray=1

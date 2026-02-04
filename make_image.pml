load halo_run/VN1R1_model.pdb, receptor
load halo_run/HaloV3_seed3.pdbqt, ligand
hide everything
show cartoon, receptor and resi 1-350
show sticks, ligand
color tv_blue, receptor
color orange, ligand
color red, ligand and name O
zoom resi 100-200+250-300
rotate y, 90
ray 1200, 800
png halo_run/winner_pose.png
quit

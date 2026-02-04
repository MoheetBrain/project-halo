load halo_run/VN1R1_model.pdb, receptor
load halo_run/HaloV3_seed3.pdbqt, ligand
hide everything
show cartoon, receptor and resi 1-350
show sticks, ligand
color spectrum, receptor
color orange, ligand
color red, ligand and name O
zoom ligand, 5.0
orient
ray 1200, 800
png halo_run/halo_v3_winner.png
quit

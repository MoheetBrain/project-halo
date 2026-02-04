load halo_run/VN1R1_model.pdb
load halo_run/HaloV3_out.pdbqt, ligand
hide everything
show cartoon, all and resi 1-350
show sticks, ligand
color tv_blue, all and not ligand
color orange, ligand
zoom ligand
ray 1200,800
png halo_run/final_winner.png
quit

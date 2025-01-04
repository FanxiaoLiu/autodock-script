#!/usr/bin/env python

from vina import Vina
import numpy as np

# Function to dock the ligand to the receptor
# Input: prepped_receptor - path to the prepped receptor pdbqt file
#        prepped_ligand - path to the prepped ligand pdbqt file
#        unprepped_pdb - path to the unprepped pdb file
#        box_size - size of each dimension of the docking box in Angstroms, an integer
#        exhaustiveness - exhaustiveness of the docking search, an integer
#        n_poses_docking - number of poses to output docking, an integer
#        n_poses_write - number of poses to write to the output file, an integer
#        output_file - path to the output file

def docking(prepped_receptor, prepped_ligand, unprepped_pdb, box_size, exhaustiveness, n_poses_docking, n_poses_write, output_file):
    
    a = np.genfromtxt(unprepped_pdb, invalid_raise=False, usecols=[6, 7, 8])
    a[np.isnan(a)] = 0
    center_array = a.mean(axis=0)
    print(center_array)

    v = Vina(sf_name='vina')

    v.set_receptor(prepped_receptor)
    v.set_ligand_from_file(prepped_ligand)

    v.compute_vina_maps(center=[-32.938, 33.253, -69.471], box_size=[box_size, box_size, box_size])

    # Dock the ligand
    v.dock(exhaustiveness, n_poses_docking)
    v.write_poses(output_file, n_poses_write, overwrite=True)

def testing():
    docking('protein-pdbqt/nmdareceptor.pdbqt', 'ligand-pdbqt/ketamine.pdbqt', 'proteins/4pe5.pdb', 50, 32, 20, 5, 'output-files/nmda-ketamine.pdbqt')

testing()
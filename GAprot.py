#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb 13 03:52:23 2023

@author: lucas
"""

import os

#from genetic_algorithm_rosetta import *
from genetic_algorithm_rosetta import genetic_algo
from apt_function import *
import apt_function

from pyrosetta import *
from rosetta.core.pack.task import TaskFactory
from rosetta.core.pack.task import operation

import numpy as np
from numpy.random import uniform
from random import sample

# Initialize PyRosetta
pyrosetta.init()

# Create a scoring function using the "ref2015_cart.wts" weight set
scorefxn = pyrosetta.create_score_function("ref2015_cart.wts")

# Creates pose object from input PDB
starting_pose = pose_from_pdb('1lzt.pdb')

# Relax the starting pose by packing and relaxing it iteratively for 3 times
starting_pose_relaxed = apt_function.pack_relax(starting_pose = starting_pose, scorefxn = scorefxn, times_to_relax = 3)
scorefxn(starting_pose)

# Define a list of single-letter amino acid codes
gene_values = ['A','C','D','E','F','G','H','I','K','L','M','N','P','Q','R','S','T','V','W','Y']
starting_pose_seq = [x for x in starting_pose.sequence()]



# List of residues to be fixed during optimization and population generation
fixed_residues_list = [1,59,60,61,62,63,64,67,90,91,92,93,94,132,133,134,135,136,137,138,139,156,157,158,159,160,161,177,180,181,182,183,184,187,188,211,212,213,215,247,263,267,268,269,270,271,272]
# Chain identifier for the fixed residues
chain = "A"

# Generate an initial population of protein structures for optimization
init_population, list_fixed_index = apt_function.Generate_random_population(starting_pose = starting_pose, 
                                                             pop_size = 5,
                                                             fixed_residues_list = fixed_residues_list,
                                                             chain = chain)


# Initiates GA object
GA = genetic_algo(pose=starting_pose, opt_direction='down',initial_population = init_population, gene_values=gene_values, gene_type='discrete',
             vector_size=len(starting_pose_seq), pop_size=len(init_population), mutation_rate=0.025, segment_fluctuation=0,
             apt_function=apt, selection_method='tournament', threads=False,
             convergence_threshold=0, n_cycles=4, tournament_cycles=int(np.round(len(init_population)/4)), tournament_size=2, benchmark=False, 
             lista_fixed=list_fixed_index, crossing_over_type='mask', file_name="teste_1.txt", dg_method="fold")

# Run the Genetic Algorithm
GA.execute()










#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb 13 03:52:23 2023

@author: lucas
"""

import os

from genetic_algorithm_rosetta import genetic_algo, genetic_algo_sequence

from apt_function import *
import apt_function

from pyrosetta import *
from rosetta.core.pack.task import TaskFactory
from rosetta.core.pack.task import operation

import numpy as np
from numpy.random import uniform
from random import sample
import random



############ Running GAPO - Structure

# Initialize PyRosetta
pyrosetta.init()

# Create a scoring function using the "ref2015_cart.wts" weight set
scorefxn = pyrosetta.create_score_function("ref2015_cart.wts")

# Creates pose object from input PDB
starting_pose = pose_from_pdb('ab_trimed_relax.pdb')
# Relax the starting pose by packing and relaxing it iteratively for 3 times
scorefxn(starting_pose)

# Define a list of single-letter amino acid codes
gene_values = ['A','C','D','E','F','G','H','I','K','L','M','N','P','Q','R','S','T','V','W','Y']
starting_pose_seq = [x for x in starting_pose.sequence()]
len(starting_pose_seq)


# Residues to mutate during optimization
CDRs = [412, 413, 414, 415, 416, 431, 432, 433, 434, 435, 436, 437, 438, 439, 440, 441, 442, 443, 444, 445, 446, 479, 480, 481, 482, 483, 484, 485, 486, 487, 488, 489, 194, 195, 196, 197, 198, 199, 200, 201, 202, 203, 204, 220, 221, 222, 223, 224, 225, 226, 259, 260, 261, 262, 263, 264, 265, 266, 267]

# Chain identifier for the fixed residues

# Generate an initial population of protein structures for optimization and return the list of residues to be locked during the evolution
init_population, list_fixed_index, listrosetta = apt_function.Generate_random_population(starting_pose = starting_pose, 
                                                              pop_size = 5,
                                                              fixed_residues_list = CDRs,
                                                              chains = ["C", "D"])
# Initiates GA object
GA = genetic_algo(pose=starting_pose, opt_direction='down',initial_population = init_population, gene_values=gene_values, gene_type='discrete',
              vector_size=len(starting_pose_seq), pop_size=len(init_population), mutation_rate=0.9, segment_fluctuation=0,
              apt_function=apt_rosetta, selection_method='tournament', threads=False,
              convergence_threshold=0, n_cycles=2, tournament_cycles=int(np.round(len(init_population)/4)), tournament_size=4, benchmark=False, 
              lista_fixed=list_fixed_index, crossing_over_type='mask', file_name="teste_esm_mut.txt", cpus  = 5, mutation_type = "esm")

# Run the Genetic Algorithm
GA.execute()


############ Running GAPO - Sequence

def generate_protein_sequence(length):
    if length < 1:
        raise ValueError("Length must be at least 1 to include the starting Methionine.")
    
    # List of amino acid residues excluding Methionine, which will be added at the start
    amino_acids = 'ACDEFGHIKLMNPQRSTVWY'
    
    # Initialize the sequence with Methionine
    sequence = 'M'
    
    # Generate the remaining sequence
    for _ in range(length - 1):
        sequence += random.choice(amino_acids)
    
    return sequence

pop_size = 50
init_population = [generate_protein_sequence(50) for i in range(pop_size)]
gene_values = ['A','C','D','E','F','G','H','I','K','L','M','N','P','Q','R','S','T','V','W','Y']
cdrs = [27, 28, 29, 30, 31, 32, 33, 34, 35, 50, 51, 52, 53, 54, 55, 56 ,57, 99, 100, 101, 102, 103, 104, 105,106,107,108,109,110,111,112,113]


# Initiates GA object
GA = genetic_algo_sequence(opt_direction='up',initial_population = init_population, gene_values=gene_values, gene_type='discrete',
             vector_size=len(init_population[0]), pop_size=len(init_population), mutation_rate=0.9, segment_fluctuation=0,
             apt_function=apt_esm, selection_method='tournament', threads=False,
             convergence_threshold=0, n_cycles=5, lista_fixed = cdrs, tournament_cycles=int(np.round(len(init_population)/4)), tournament_size=4, benchmark=False,  
             crossing_over_type='mask', file_name="teste_new_algo.txt", mutation_type = "esm")

# Run the Genetic Algorithm
GA.execute()







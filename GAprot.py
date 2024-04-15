#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb 13 03:52:23 2023

@author: lucas
"""

import os

from genetic_algorithm_rosetta import *
from apt_function import *
import apt_function

from pyrosetta import *
from rosetta.core.pack.task import TaskFactory
from rosetta.core.pack.task import operation

import numpy as np
from numpy.random import uniform
from random import sample

pyrosetta.init()


scorefxn = pyrosetta.create_score_function("ref2015_cart.wts")
starting_pose = pose_from_pdb('1lzt.pdb')
scorefxn(starting_pose)
starting_pose_seq = [x for x in starting_pose.sequence()]
gene_values = ['A','C','D','E','F','G','H','I','K','L','M','N','P','Q','R','S','T','V','W','Y']
scorefxn(starting_pose)
pop_size = 100
vector_size = len(starting_pose_seq)
init_pop = [list(starting_pose_seq)]

positions_starting_pose = list(range(0,vector_size))

#### Residues to be fixed during optimization and new population generation
fixed_residues_list = []


lista_fixed_rosetta = apt_function.PDB_to_Pose(starting_pose, fixed_residues_list, "A")

#### Residues to be fixed during optimization and new population generation

#lista_fixed_rosetta = [x - 1 for x in lista_fixed_rosetta]
new_indiv = [[gene_values[int(np.round(uniform(low=0, high=len(gene_values)-1)))] for i in range(vector_size)] for indiv in range(pop_size)]

#### Fixing the fixed positions to the new inds
for individual in new_indiv:
    for i in lista_fixed_rosetta:
        individual[i-1] = starting_pose_seq[i-1]

init_popu = list(new_indiv)


### If using threads, use apt_threads, if not use APT


GA = genetic_algo(pose=starting_pose, opt_direction='down',initial_population = init_popu, gene_values=gene_values, gene_type='discrete',
             vector_size=len(starting_pose_seq), pop_size=len(init_pop), mutation_rate=0.025, segment_fluctuation=0,
             apt_function=apt, selection_method='tournament', threads=False,
             convergence_threshold=0, n_cycles=150, tournament_cycles=int(np.round(pop_size/4)), tournament_size=4, benchmark=False, lista_fixed=lista_fixed_rosetta,crossing_over_type='mask', file_name="cond_1.txt", dg_method="fold")


GA.execute()

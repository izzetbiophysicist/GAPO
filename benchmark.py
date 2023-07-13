#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr 11 13:38:06 2023

@author: lucas
"""

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb 13 03:52:23 2023

@author: lucas
"""

import os


os.chdir('/home/lucas/genetic_algo/') 

from genetic_algorithm_rosetta import *
from apt_function import *

from pyrosetta import *
from rosetta.core.pack.task import TaskFactory
from rosetta.core.pack.task import operation

import numpy as np
from numpy.random import uniform
from random import sample

pyrosetta.init()


gene_values = ['A','C','D','E','F','G','H','I','K','L','M','N','P','Q','R','S','T','V','W','Y']

GA = genetic_algo(pose=[],opt_direction='down', gene_values=gene_values, gene_type='discrete', 
             vector_size=20, pop_size=20, mutation_rate=0.6, segment_fluctuation=0, 
             apt_function=apt_benchmark, selection_method='tournament', threads=False,
             convergence_threshold=0, n_cycles=100, tournament_cycles=4, tournament_size=4,benchmark=True, crossing_over_type='mask')

GA.execute()
print(GA.score_history)

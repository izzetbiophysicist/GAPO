# GAPO - Genetic algorithm for Protein Optimization

<div align="center">
  <img src="./GAPO_logo.png" alt="Slinky Logo" width="400">
</div>

GAPO is a _in silico_ genetic algorithm used to optimize proteins for a desired function, such as stability and affinity. The algorithm mimetizes the evolutionary process, by recombining and adding mutations to the best sequences in order to generate a new population with higher diversity and optimized for the given objetive function. 
<br/>

<p align="center">
  <em>In silico protein optimization through genetic algorithms.</em>
  <br />
  <br />
  <a href="#about-the-project">About The Project</a> •
  <a href="#getting-started">Getting Started</a> •
  <a href="#usage">Usage</a> •
  <a href="#contributing">Contributing</a>
</p>

---

### Table of Contents

- [About The Project](#about-the-project)
- [Getting Started](#getting-started)
  - [Prerequisites](#prerequisites)
  - [Installation](#installation)
- [Usage](#usage)
  - [Example 1: Structure-Based Optimization](#example-1-structure-based-optimization)
  - [Example 2: Sequence-Based Optimization](#example-2-sequence-based-optimization)
  - [Algorithm Parameters](#algorithm-parameters)
- [Contributing](#contributing)
- [License](#license)
- [Contact](#contact)

# Installing dependencies

## Getting started

## Download & Install
#### Step 1 - Install Conda (if not already installed)

If Conda is not already installed on your system, follow the instructions below to download and install it:

1. Visit the [Anaconda installation page](https://www.anaconda.com/download/success).
2. Download the Anaconda installer for your operating system.
3. Follow the installation instructions on the page to complete the setup.

After installing Anaconda, you should be able to use `conda` commands in your terminal.

#### Step 2 - Clone the repository

```
git clone https://github.com/izzetbiophysicist/prot_eng_GA.git
```
#### Step 3 - Create gapo_env in Conda
```
cd prot_eng_GA
conda env create -f environment.yml
conda activate gapo_env
```
As a prerequisite, you must have PyTorch installed to use this repository. In order to install the right Pytorch version to your machine, please visit [Pytorch installation page](https://pytorch.org/get-started/locally/).
#### Step 4 - Install PyRosetta in gapo_env
```
pip install pyrosetta-installer
python -c 'import pyrosetta_installer; pyrosetta_installer.install_pyrosetta()'
```


## Running the code

genetic_algorithm.py = contains main function
<br/>
apt_functions.py = contains a collection of objective functions for the genetic algorithm
<br/>
GAprot.py imports genetic_algorithmpy and apt_functions.py, configures initial population and carries out the optimization
<br/>
### GA parameters and descriptions
  
  | Parameter | Description  | 
  | :---:   | :---: |
  |opt_direction | Selects if the objective function goes up or down
  |apt_function | Select between "apt_thread" for pbee and "apt_rosetta" for rosetta.
  |gene_values | Values that vectors can assume
  |gene_type | Selects between discrete or continuous genes
  |vector_size | Size of the vectors
  |selection_methods | Elitist or tournament selection
  |n_survivor | Number of survivors in each selection
  |tournament_size | Number of individuals selected in each tournament
  |crossing_over_type | Type of crossover during optimization
  |lista_fixed | List of residues to maintain during evolution
  |initial_population | An initial population can be given. Otherwise, a random population is created
  |file_name | Output log file, Sequences - dG - population
  |cpus | Numbers of CPU usage to paralelize

<br/>

### GA function
<br/>
Follow below a usage example of the function used to optimize the CDRs from a scFv complexed to CD19 (script GAprot.py)
<br/>

```python

####### Loading libraries
import os
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
starting_pose = pose_from_pdb('CD19_scFv_relax.pdb')
# Relax the starting pose by packing and relaxing it iteratively for 3 times

starting_pose = apt_function.pack_relax(starting_pose = starting_pose, scorefxn = scorefxn, times_to_relax = 3)
scorefxn(starting_pose)

# Define a list of single-letter amino acid codes
gene_values = ['A','C','D','E','F','G','H','I','K','L','M','N','P','Q','R','S','T','V','W','Y']
starting_pose_seq = [x for x in starting_pose.sequence()]



# Residues to mutate during optimization
CDRs = [62,63,64,65,66,67,68,69,70,71,72,88,89,90,91,92,93,94,127,128,129,130,131,132,133,134,135,186,187,188,189,190,191,192,212,213,214,215,216,257,258,259,260,261,262,263,264,265,266,267,268,269]# Chain identifier for the fixed residues


# Generate an initial population of protein structures for optimization and return the list of residues to maintain locked during the evolution
# Pop size receives the size of the population to be generated
# fixed_residues_list receives the list of residues to mutate during the evolution
# chains receive a list of the chains from the respective residues to mutate

init_population, list_fixed_index = apt_function.Generate_random_population(starting_pose = starting_pose, 
                                                             pop_size = 50,
                                                             fixed_residues_list = CDRs,
                                                             chains = ["C", "D"])

GA = genetic_algo(pose=starting_pose, opt_direction='down',initial_population = init_population, gene_values=gene_values, gene_type='discrete',
             vector_size=len(starting_pose_seq), pop_size=len(init_population), mutation_rate=0.025, segment_fluctuation=0,
             apt_function=apt, selection_method='tournament', threads=False,
             convergence_threshold=0, n_cycles=4, tournament_cycles=int(np.round(len(init_population)/4)), tournament_size=4, benchmark=False, 
             lista_fixed=list_fixed_index, crossing_over_type='mask', file_name="teste_1.txt", dg_method="bind",  cpus  = 5)
GA.execute()
```



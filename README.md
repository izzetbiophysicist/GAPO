# GAPEng - Genetic algorithm for protein engineering
In silico protein engineering using Genetic Algorithms and rosetta energy calculations

genetic_algorithm.py = contains main function
apt_functions.py = contains a collection of objective functions for the genetic algorithm

GAprot.py imports genetic_algorithmpy and apt_functions.py, configures initial population and carries out the optimization

-Threads not supported yet

### GA function

GA = genetic_algo(pose=starting_pose, opt_direction='down',initial_population = init_population, gene_values=gene_values, gene_type='discrete',
             vector_size=len(starting_pose_seq), pop_size=len(init_population), mutation_rate=0.025, segment_fluctuation=0,
             apt_function=apt, selection_method='tournament', threads=False,
             convergence_threshold=0, n_cycles=4, tournament_cycles=int(np.round(len(init_population)/4)), tournament_size=2, benchmark=False, 
             lista_fixed=list_fixed_index, crossing_over_type='mask', file_name="teste_1.txt", dg_method="fold")
             
  pose = starting PDB structure
  
  opt_direction = selects if the objective function goes up or down
  
  gene_values = Values that vectors can assume
  
  gene_type = Selects between discrete or continuous genes
  
  vector_size = size of the vectors
  
  selection methods = elitist or tournament selection
  
  threahds = parallel processing (not supported)
  
  n_cycles = number of optimization cycles
  
  n_survivor = number of survivors in each selection
  
  tournament_size =number of individuals selected in each tournament
  
  initial_population = An initial population can be given. Otherwise, a random population is created

  file_name= output log file, Sequences - dG - population
  
  dg_method = select between "fold" and "bind".

# Installing dependencies

## Getting started

## Clone the reposity
```
git clone https://github.com/izzetbiophysicist/prot_eng_GA.git
```
## Downloading pyRosetta
First of all, you must download PyRosetta. To download , you need to get a license.
<br />
License and Downloads links:
<br />
[License](https://www.rosettacommons.org/software/license-and-download)
<br />
[PyRosetta4 Download](https://graylab.jhu.edu/download/PyRosetta4/archive/release/PyRosetta4.Release.python310.linux/PyRosetta4.Release.python310.linux.release-370.tar.bz2)


## Creating a conda env to run GA
```
conda create --name GA_env --file requirements.txt
```
## Installing PyRosetta in Conda ENV
### After downloading, unzip PyRosetta's and enter the setup directory to install it
```
tar -xvf PyRosetta4.Release.python310.linux.release-370.tar.bz2
cd PyRosetta4.Release.python310.linux.release-370/setup
conda activate GA_env
python setup.py install
```

Additional help for downloading and installing and PyRosetta (source:Sari Sabban youtube channel )

[![IMAGE ALT TEXT HERE](https://img.youtube.com/vi/UEaFmUMEL9c/0.jpg)](https://www.youtube.com/watch?v=UEaFmUMEL9c)


# Genetic algorithm for protein optimization
In silico protein engineering using Genetic Algorithms and rosetta energy calculations

genetic_algorithm.py = main function
apt_functions.py = contains the objective functions

GAprot.py imports both in order to run the optimization

-Threads not supported yet

###########
### GA function

genetic_algo(pose=starting_pose, opt_direction='down', gene_values=gene_values, gene_type='discrete', 
             vector_size=len(starting_pose_seq), pop_size=len(init_pop), mutation_rate=0.2, segment_fluctuation=0, 
             apt_function=apt, selection_method='tournament', threads=False,
             convergence_threshold=0, n_cycles=20, n_survivors=1, tournament_size=4,
             initial_population=init_pop)
             
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
  

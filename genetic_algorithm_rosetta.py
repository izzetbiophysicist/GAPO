import numpy as np
from numpy.random import uniform
from random import sample

from threading import Thread
from time import sleep
from pyrosetta import *
from rosetta.core.pack.task import TaskFactory
from rosetta.core.pack.task import operation
from datetime import datetime
import pandas as pd
from apt_function import correct_multi_input, batchs_to_run

#### ESM stuff
import torch
import esm
import random
import math

# Load ESM-2 model
model, alphabet = esm.pretrained.esm2_t33_650M_UR50D()
batch_converter = alphabet.get_batch_converter()
model.eval()  # disables dropout for deterministic results
def insert_mask(sequence, position, mask="<mask>"):
    """
    Replaces a character in a given position of a sequence with a mask.

    Parameters:
    - sequence (str or list): The sequence to replace the character in.
    - position (int): The position in the sequence where the character should be replaced.
    - mask (str): The mask to insert (default is "<mask>").

    Returns:
    - str or list: The sequence with the mask replacing the character at the specified position.
    """
    
    if not (0 <= position < len(sequence)):
        raise ValueError("Position is out of bounds.")
    
    if isinstance(sequence, str):
        return sequence[:position] + mask + sequence[position + 1:]
    elif isinstance(sequence, list):
        return sequence[:position] + [mask] + sequence[position + 1:]
    else:
        raise TypeError("Sequence must be a string or list.")


def complete_mask(input_sequence, posi, temperature=1.0):
    
    
    exclude = [alphabet.get_idx('X'), alphabet.get_idx('B'),
               alphabet.get_idx('U'), alphabet.get_idx('Z'),
               alphabet.get_idx('O')]
    
    data = [
        ("protein1", insert_mask(input_sequence, posi, mask="<mask>"))
    ]

    batch_labels, batch_strs, batch_tokens = batch_converter(data)
    batch_lens = (batch_tokens != alphabet.padding_idx).sum(1)

    # Predict masked tokens
    with torch.no_grad():
        token_probs = model(batch_tokens, repr_layers=[33])["logits"]

    # Apply temperature
    token_probs /= temperature
    
    softmax = torch.nn.Softmax(dim=-1)
    probabilities = softmax(token_probs)

    # Get the index of the <mask> token
    mask_idx = (batch_tokens == alphabet.mask_idx).nonzero(as_tuple=True)
    
    # Zero out probabilities for excluded tokens
    for token in exclude:
        probabilities[:, :, token] = 0.0

    # Sample from the probability distribution
    predicted_tokens = torch.multinomial(probabilities[mask_idx], num_samples=1).squeeze(-1)

    # Replace the <mask> token with the predicted token
    batch_tokens[mask_idx] = predicted_tokens


    predicted_residues = [alphabet.get_tok(pred.item()) for pred in batch_tokens[0]]
    
    seq_predicted = ''.join(predicted_residues[1:-1])
    
    if input_sequence != seq_predicted:
        print("Mutation added!! 😉")
    
    return seq_predicted
class GeneticAlgoBase:
    def __init__(self, opt_direction, gene_values, gene_type, vector_size, threads, pop_size, mutation_rate, segment_fluctuation, apt_function, selection_method, convergence_threshold, n_cycles, benchmark, crossing_over_type, tournament_cycles, file_name, mutation_type, esm_tmp=1.0, initial_population=[], lista_fixed=[], tournament_size=2):
        self.initial_population = initial_population
        self.population = initial_population
        self.gene_values = gene_values  ### if gene_type = 'continuous' then gene_values should contain the upper and lower bounds
        self.pop_size = pop_size
        self.mutation_rate = mutation_rate
        self.segment_fluctuation = segment_fluctuation
        self.apt_function = apt_function
        self.selection_method = selection_method
        self.tournament_size = tournament_size
        self.tournament_cycles = tournament_cycles
        self.convergence_threshold = convergence_threshold
        self.n_cycles = n_cycles
        self.gene_type = gene_type
        self.opt_direction = opt_direction
        self.vector_size = vector_size
        self.threads = threads
        self.benchmark = benchmark
        self.crossing_over_type = crossing_over_type
        self.file_name = file_name
        self.lista_fixed = lista_fixed
        self.t = "init"
        self.mutation_type = mutation_type
        self.esm_tmp = esm_tmp

        if len(self.population) != 0:
            self.pop_size = len(self.population)

    def initialize_population(self):
        if self.gene_type == 'discrete':
            self.population = [[self.gene_values[int(np.round(uniform(low=0, high=len(self.gene_values)-1)))] for i in range(self.vector_size)] for indiv in range(self.pop_size)]
            self.scores = ['None' for _ in range(self.pop_size)]
            self.first_population = self.population

        elif self.gene_type == 'continuous':
            self.population = [[uniform(low=self.gene_values[0], high=self.gene_values[1])] for _ in range(self.vector_size) for _ in range(self.pop_size)]
            self.scores = ['None' for _ in range(self.pop_size)]
            self.first_population = self.population

    def calculate_scores(self, population, pre_calc=[]):
        if len(pre_calc) == 0:
            if self.benchmark:
                print('CALCULATING SCORES!')
                scores = [self.apt_function(population[x]) for x in range(len(population))]
            else:
                if not self.threads:
                    print('CALCULATING SCORES!')
                    population = correct_multi_input(population)
                    scores = batchs_to_run(self.pose, self.apt_function, population, self.cpus, self.t)
                else:
                    print('CALCULATING SCORES!')
                    t1 = thread_rosetta(population, self.pose, self.scorefxn, self.apt_function, self.dg_method)
                    t1.run()
                    scores = [t1.return_results[x][0] for x in range(len(t1.return_results))]

        else:
            scores = list(pre_calc)
            scores_to_append = [self.apt_function(population[x]) for x in range(len(population)) if x >= len(pre_calc)]
            scores = scores + scores_to_append

        return scores

    def crossing_over(self, ind1, ind2, crossing_over_type):
        if crossing_over_type == 'punctual':
            if self.segment_fluctuation + (len(ind1) / 2) > len(ind1):
                print('segment fluctuation is too long')

            fluct = int(np.round(uniform(low=0, high=self.segment_fluctuation)))
            newind1 = ind1[0:int((len(ind1) / 2) + fluct)] + ind2[int((len(ind2) / 2) + fluct):]
            newind2 = ind2[0:int((len(ind2) / 2) + fluct)] + ind1[int((len(ind1) / 2) + fluct):]

        elif crossing_over_type == 'mask':
            mask = [int(np.round(uniform(low=0, high=1))) for _ in range(self.vector_size)]
            newind1 = [ind1[x] if mask[x] == 1 else ind2[x] for x in range(len(ind2))]
            newind2 = [ind1[x] if mask[x] != 1 else ind2[x] for x in range(len(ind2))]

        elif crossing_over_type == 'multi_segment':
            num_segments = np.random.randint(1, self.max_segments + 1)
            segment_lengths = np.diff(np.sort(np.random.choice(range(1, len(ind1)), num_segments - 1, replace=False)))
            segment_lengths = np.insert(segment_lengths, 0, np.random.randint(1, len(ind1) // num_segments))
            segment_lengths = np.append(segment_lengths, len(ind1) - np.sum(segment_lengths))

            newind1, newind2 = [], []
            start = 0
            for i, seg_len in enumerate(segment_lengths):
                end = start + seg_len
                if i % 2 == 0:
                    newind1.extend(ind1[start:end])
                    newind2.extend(ind2[start:end])
                else:
                    newind1.extend(ind2[start:end])
                    newind2.extend(ind1[start:end])
                start = end

        return newind1, newind2

    def mutate(self, ind, lista_fixed, mutation_type):
        lista_len_seq = list(range(1, len(ind)))
        inds_to_mut = [i for i in lista_len_seq if i not in lista_fixed]
        if mutation_type == "esm":
            position = int(np.random.choice(inds_to_mut, 1))
            ind = list(complete_mask(''.join(ind), position, self.esm_tmp))
        if mutation_type == "random":
            position = int(np.random.choice(inds_to_mut, 1))
            gene = self.gene_values[int(np.round(uniform(low=0, high=(len(self.gene_values)-1))))]
            ind[position] = gene
        return ind

    def breed(self, population, scores):
        if self.selection_method == 'tournament':
            pop_scores = {'indv': population, 'score': scores}
            init_pop = pd.DataFrame(pop_scores)

            if self.opt_direction == 'up':
                init_pop = init_pop.sort_values('score', ascending=False)
            elif self.opt_direction == 'down':
                init_pop = init_pop.sort_values('score', ascending=True)

            offspring = []
            for _ in range(self.tournament_cycles):
                to_select = init_pop.iloc[sample(range(len(population)), self.tournament_size)]
                to_select = to_select.sort_values('score', ascending=(self.opt_direction == 'down'))

                co = self.crossing_over(to_select.iloc[0][0], to_select.iloc[1][0], self.crossing_over_type)
                for newindv in co:
                    offspring.append(newindv)

            for trymut in range(len(offspring)):
                if uniform(low=0, high=1) < self.mutation_rate:
                    offspring[trymut] = self.mutate(offspring[trymut], self.lista_fixed, self.mutation_type)

            offspring_scores = self.calculate_scores(offspring)

            whole_pop = {'indv': population + offspring, 'score': scores + offspring_scores}
            whole_pop = pd.DataFrame(whole_pop)
            whole_pop = whole_pop.sort_values('score', ascending=(self.opt_direction == 'down'))

            return whole_pop['indv'].to_list(), whole_pop['score'].to_list()

    def write_out(self, t):
        with open(self.file_name, 'a') as file:
            for indiv_index in range(len(self.population)):
                h_individual = ''.join(self.population[indiv_index])
                h_score = str(self.scores[indiv_index])
                h_population = t
                final_text = ','.join([h_individual, h_score, str(h_population)])
                file.write(final_text + '\n')
    def opt_cycle(self):

        ### Initialize population if none was given
        if len(self.initial_population) == 0:
            self.initialize_population()
        self.scores = self.calculate_scores(self.population)
        self.initial_scores = self.scores

        self.score_history = [self.initial_scores]
        self.pop_history = []
        self.best = []
        self.best_ind = []
        self.start_time = datetime.now()
        print('start')
        self.write_out(0)

        for t in range(self.n_cycles):
            print('Running round '+str(t))

            new_pop = self.breed(self.population, self.scores)

            self.population = new_pop[0]
            self.scores = new_pop[1]


            #### Tracking variables
            self.pop_history.append(self.population)
            self.t=t
            print(t)
            if self.opt_direction == "up":
                self.best_ind.append(self.population[np.argmax(self.scores)])
                self.best.append(self.scores[np.argmax(self.scores)])
            else:
                self.best_ind.append(self.population[np.argmin(self.scores)])
                self.best.append(self.scores[np.argmin(self.scores)])
            self.pop_history.append(self.population)
            self.score_history.append(self.scores)

            ### write population
            self.write_out(t+1)

        self.finish_time = datetime.now()
        self.exec_time = self.finish_time - self.start_time


    def execute(self):
        self.opt_cycle()

class genetic_algo(GeneticAlgoBase):
    def __init__(self, pose, opt_direction, gene_values, mutation_type, gene_type, vector_size, threads, pop_size, mutation_rate, segment_fluctuation, apt_function, selection_method, convergence_threshold, n_cycles, benchmark, crossing_over_type, tournament_cycles, file_name, lista_fixed, cpus, tournament_size=2, esm_tmp=1.0, initial_population=[]):
        super().__init__(opt_direction, gene_values, gene_type, vector_size, threads, pop_size, mutation_rate, segment_fluctuation, apt_function, selection_method, convergence_threshold, n_cycles, benchmark, crossing_over_type, tournament_cycles, file_name, mutation_type, esm_tmp, initial_population, lista_fixed, tournament_size)
        self.pose = pose
        self.scorefxn = pyrosetta.create_score_function("ref2015_cart.wts")
        self.cpus = cpus

    def calculate_scores(self, population, pre_calc=[]):

        if len(pre_calc) == 0:

            if self.benchmark==True:
                print('CALCULATING SCORES!')

            #### Iterate over population and fill self.scores
                scores = [self.apt_function(population[x]) for x in range(len(population))]

            else:
                if self.threads==False:
                    print('CALCULATING SCORES!')

                #### Iterate over population and fill self.scores
                    #scores = [self.apt_function(population[x], self.pose, self.scorefxn, x, self.t) for x in range(len(population))]
                    population = correct_multi_input(population)
                    scores = batchs_to_run(self.pose, self.apt_function, population, self.cpus, self.t)
                    #scores = self.apt_function(self.pose, population, self.dg_method, self.cpus, self.t)


        if len(pre_calc) != 0:

            if self.benchmark == True:
                scores = list(pre_calc)
                scores_to_append = [self.apt_function(population[x]) for x in range(len(population)) if x >= len(pre_calc)]
                scores = scores + scores_to_append

            else:

                if self.threads == False:
                    scores = list(pre_calc)
                    scores_to_append = [self.apt_function(population[x], self.pose, self.scorefxn) for x in range(len(population)) if x >= len(pre_calc)]
                    scores = scores + scores_to_append

                if self.threads == True:
                    scores = list(pre_calc)
                    t1 = thread_rosetta([self.population[x] for x in range(len(population)) if x >= len(pre_calc)], self.pose, self.scorefxn, self.apt_function)
                    t1.run()

                    scores_to_append = [t1.return_results[x][0] for x in range(len(t1.return_results))]


                scores = scores + scores_to_append


        return scores
    
class genetic_algo_sequence(GeneticAlgoBase):
    def __init__(self, opt_direction, gene_values, mutation_type, gene_type, vector_size, threads, pop_size, mutation_rate, segment_fluctuation, apt_function, selection_method, convergence_threshold, n_cycles, benchmark, crossing_over_type, tournament_cycles, file_name, lista_fixed, tournament_size=2, esm_tmp=1.0, initial_population=[]):
        super().__init__(opt_direction, gene_values, gene_type, vector_size, threads, pop_size, mutation_rate, segment_fluctuation, apt_function, selection_method, convergence_threshold, n_cycles, benchmark, crossing_over_type, tournament_cycles, file_name, mutation_type, esm_tmp, initial_population, lista_fixed, tournament_size)
        

    def calculate_scores(self, population, pre_calc=[]):

        if len(pre_calc) == 0:

            if self.benchmark==True:
                print('CALCULATING SCORES!')

            #### Iterate over population and fill self.scores
                scores = [self.apt_function(population[x]) for x in range(len(population))]

            else:
                if self.threads==False:
                    print('CALCULATING SCORES!')

                #### Iterate over population and fill self.scores
                    
                    population = correct_multi_input(population)
                    scores = [self.apt_function(population[x]) for x in range(len(population))]
                    
                    #scores = batchs_to_run(self.pose, population, self.dg_method, self.cpus, self.t)
                    
                ## TEST APT-FUNCTION WHEN NOT USING THREADS

                if self.threads==True:
                    print('CALCULATING SCORES!')
                    t1 = thread_rosetta(population, self.pose, self.scorefxn, self.apt_function, self.dg_method)
                    t1.run()

                    scores = [t1.return_results[x][0] for x in range(len(t1.return_results))]


        if len(pre_calc) != 0:

            if self.benchmark == True:
                scores = list(pre_calc)
                scores_to_append = [self.apt_function(population[x]) for x in range(len(population)) if x >= len(pre_calc)]
                scores = scores + scores_to_append

            else:

                if self.threads == False:
                    scores = list(pre_calc)
                    scores_to_append = [self.apt_function(population[x], self.pose, self.scorefxn) for x in range(len(population)) if x >= len(pre_calc)]
                    scores = scores + scores_to_append

                if self.threads == True:
                    scores = list(pre_calc)
                    t1 = thread_rosetta([self.population[x] for x in range(len(population)) if x >= len(pre_calc)], self.pose, self.scorefxn, self.apt_function)
                    t1.run()

                    scores_to_append = [t1.return_results[x][0] for x in range(len(t1.return_results))]


                scores = scores + scores_to_append


        return scores

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Feb 10 18:34:46 2023

@author: lucas
"""

import numpy as np
from numpy.random import uniform
from random import sample

class genetic_algo:
    def __init__(self, opt_direction, gene_values, gene_type, vector_size, pop_size, mutation_rate, segment_fluctuation, apt_function, selection_method, convergence_threshold, n_cycles, n_survivors, tournament_size=0, initial_population=[]):
        self.population = initial_population
        self.gene_values = gene_values  ### if gene_type = 'continuous' then gene_values should contain the upper and lower bounds
        self.pop_size = pop_size
        self.mutation_rate = mutation_rate
        self.segment_fluctuation = segment_fluctuation
        self.apt_function = apt_function
        self.selection_method = selection_method
        self.tournament_size = tournament_size
        self.convergence_threshold = convergence_threshold
        self.n_cycles = n_cycles
        self.gene_type = gene_type
        self.n_survivors = n_survivors
        self.opt_direction = opt_direction
        self.vector_size = vector_size

    
        ### recalculate pop_size if the initial population is not randomly generated
        if(len(self.population) != 0):
            self.pop_size = len(self.population)
        
    def initialize_population(self):
        
        #### select gene type
        if(self.gene_type == 'discrete'):
            
            self.population = [[self.gene_values[int(np.round(uniform(low=0, high=len(self.gene_values)-1)))] for i in range(self.vector_size)] for indiv in range(self.pop_size)]
            self.scores = ['None' for i in range(self.pop_size)]
            self.first_population = self.population
            
        if(self.gene_type == 'continuous'):
            self.population = [[uniform(low=self.gene_values[0], high=self.gene_values[1])] for i in range(self.vector_size) for indiv in range(self.pop_size)]
            self.scores = ['None' for i in range(self.pop_size)]
            self.first_population = self.population
            
            
            
    def calculate_scores(self, population, pre_calc=[]):
        if len(pre_calc) == 0:
        #### Iterate over population and fill self.scores
            scores = [self.apt_function(population[x]) for x in range(len(population))]
        else:
            scores = list(pre_calc) 
            scores_to_append = [self.apt_function(population[x]) for x in range(len(population)) if x >= len(pre_calc)]
            scores = scores + scores_to_append
        return scores
    
   
    def select(self, population, scores):
        
        ##### Tournament rounds take N individuals and select them. When N > the remaining population (not subjected to selection yet)
        ##### the remaining individuals follow to the next step without selection
        
            
        
        
        selected_ind = []
        selected_scores = []
        
        dead_ind = []
        dead_scores = []
        
        if self.selection_method == 'tournament':
            
            tmp_pop = list(population)
            tmp_scores = list(scores)
            
            while len(tmp_pop) > self.tournament_size:
                to_select = sample(range(len(tmp_pop)), self.tournament_size)
                
                tournament_scores = [tmp_scores[x] for x in to_select]
                tournament_pop = [tmp_pop[x] for x in to_select]
                
                for rounds in range(self.n_survivors):
                    
                    if self.opt_direction == 'up':
                        selected_ind = selected_ind + [tournament_pop[np.argmax(tournament_scores)]]
                        selected_scores.append(tournament_scores[np.argmax(tournament_scores)]) ## append here
                    if self.opt_direction == 'down':
                        selected_ind = selected_ind + [tournament_pop[np.argmin(tournament_scores)]]
                        selected_scores.append(tournament_scores[np.argmin(tournament_scores)]) ## append here
                        
                    if self.opt_direction == 'up':
                        tournament_pop.pop(np.argmax(tournament_scores))                    
                        tournament_scores.pop(np.argmax(tournament_scores))
                    
                    if self.opt_direction == 'down':
                        tournament_pop.pop(np.argmin(tournament_scores))                    
                        tournament_scores.pop(np.argmin(tournament_scores))
                    
                    
                dead_ind = dead_ind + tournament_pop
                dead_scores = dead_scores + tournament_scores
                
                
                tmp_scores = [tmp_scores[x] for x in range(len(tmp_scores)) if x not in to_select]
                tmp_pop = [tmp_pop[x] for x in range(len(tmp_pop)) if x not in to_select]
            
            return selected_ind, selected_scores, dead_ind, dead_scores
            
        if self.selection_method == 'elitist':
            selected_ind = []
            selected_scores = []
            
            dead_ind = []
            dead_scores = []
            
            tmp_pop = list(population)
            tmp_scores = list(scores)
            
            while len(tmp_pop) != self.n_survivors:
                
                
                if self.opt_direction == 'up':

                    selected_ind = selected_ind + [tmp_pop[np.argmax(tmp_scores)]]
                    selected_scores.append(tmp_scores[np.argmax(tmp_scores)]) ## append here
                    
                    tmp_pop.pop(np.argmax(tmp_scores))                    
                    tmp_scores.pop(np.argmax(tmp_scores))
                
                if self.opt_direction == 'down':

                    selected_ind = selected_ind + [tmp_pop[np.argmin(tmp_scores)]]
                    selected_scores.append(tmp_scores[np.argmin(tmp_scores)]) ## append here
                    
                    tmp_pop.pop(np.argmin(tmp_scores))                    
                    tmp_scores.pop(np.argmin(tmp_scores))
                    
                dead_ind = list(tmp_pop)
                dead_scores = list(tmp_scores)

            
            return selected_ind, selected_scores, dead_ind, dead_scores
        
    
    def crossing_over(self, ind1, ind2):
        if self.segment_fluctuation+(len(ind1)/2) > len(ind1):
            print('segment fluctuation is too long') 
            
        #sample(segment_fluctuation, 1)
        fluct=int(np.round(uniform(low=0, high=self.segment_fluctuation)))
        newind1 = ind1[0:int((len(ind1)/2)+fluct)] + ind2[(int((len(ind2)/2)+fluct)):len(ind2)+1]
        newind2 = ind2[0:int((len(ind2)/2)+fluct)] + ind1[(int((len(ind1)/2)+fluct)):len(ind1)+1]        
        
        return newind1, newind2
    
    def mutate(self, ind):
        #### sample gene and value and mutate
        #### use mutation_rate as probability
        if self.gene_type == 'discrete':
            position=int(np.round(uniform(low=0, high=len(ind)-1)))
            gene = self.gene_values[int(np.round(uniform(low=0, high=(len(self.gene_values)-1))))]
            ind[position] = gene
            
        if self.gene_type == 'continuous':
            position=int(np.round(uniform(low=0, high=len(ind)-1)))
            gene = [uniform(low=self.gene_values[0], high=(self.gene_values[1]-1))]
            ind[position] = gene
            
        return ind
        
    
    
    def repopulate(self, selected, selected_scores, dead, dead_scores):
        
        
        new_population = selected
        
        
        while len(new_population) != self.pop_size:
            
            ind1=selected[int(np.round(uniform(low=[0], high=(len(selected)-1))))]
            ind2=dead[int(np.round(uniform(low=[0], high=(len(dead)-1))))]
            
            offspring = self.crossing_over(ind1, ind2)
            
            i=0
            for i in range(2):
                if len(new_population) != self.pop_size:
                    new_population.append(list(offspring[i]))
                    if uniform(low=0, high=1) < self.mutation_rate:
                        new_population[len(new_population)-1] = self.mutate(new_population[len(new_population)-1])
                
                i=+1
            
            
        return new_population
        
            
                
    def opt_cycle(self):
        
        self.initialize_population()
        self.scores = self.calculate_scores(self.population)
        self.initial_scores = self.scores
        
        
        self.pop_history = []
        self.best = []
        self.best_ind = []
        print('start')

        
        for t in range(self.n_cycles):
            select_round = self.select(self.population, self.scores)

            self.population = self.repopulate(selected=select_round[0], selected_scores=select_round[1], dead=select_round[2], dead_scores=select_round[3])
            self.pop_history.append(self.population)

            self.scores = self.calculate_scores(self.population, pre_calc=select_round[1])

            
            self.t=t
            print(t)
            self.best.append(self.scores[np.argmax(self.scores)])
            self.best_ind.append(self.population[np.argmax(self.scores)])
            self.pop_history.append(self.population)
    
    def execute(self):
        self.opt_cycle()
        
        
 

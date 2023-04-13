#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Feb 10 18:34:46 2023

@author: lucas
"""

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

pyrosetta.init()

class thread_rosetta:
    def __init__(self, pop, starting_pose, scorefxn, apt_function):
        self.threads = list()
        self.return_results = list()
        
        for i in range(len(pop)):
            print(i)
            self.return_results.append([None]*1)
            self.threads.append(Thread(target=apt_function, args=(pop[i], starting_pose, scorefxn, self.return_results[i])))

    def run(self):
        for i in range(len(self.threads)):
            print(i)
            sleep(1)
            self.threads[i].start()
            sleep(1)



class genetic_algo:
    def __init__(self, pose, opt_direction, gene_values, gene_type, vector_size, threads, pop_size, mutation_rate, segment_fluctuation, apt_function, selection_method, convergence_threshold, n_cycles, benchmark, crossing_over_type, tournament_size=0, initial_population=[]):
        self.initial_population=initial_population 
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
        self.opt_direction = opt_direction
        self.vector_size = vector_size
        self.pose = pose
        self.scorefxn = pyrosetta.create_score_function("ref2015_cart.wts")
        self.threads = threads
        self.benchmark = benchmark
        self.crossing_over_type = crossing_over_type
        
        
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
        
            if self.benchmark==True:
                print('CALCULATING SCORES!')
    
            #### Iterate over population and fill self.scores
                scores = [self.apt_function(population[x]) for x in range(len(population))]
            
            else:
                if self.threads==False:
                    print('CALCULATING SCORES!')
        
                #### Iterate over population and fill self.scores
                    scores = [self.apt_function(population[x], self.pose, self.scorefxn) for x in range(len(population))]
                
                ## TEST APT-FUNCTION WHEN NOT USING THREADS
                
                if self.threads==True:
                    print('CALCULATING SCORES!')
                    t1 = thread_rosetta(population, self.pose, self.scorefxn, self.apt_function)
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
    
    
    def crossing_over(self, ind1, ind2, crossing_over_type):
        
        if crossing_over_type == 'punctual':
            if self.segment_fluctuation+(len(ind1)/2) > len(ind1):
                print('segment fluctuation is too long') 
                
            #sample(segment_fluctuation, 1)
            fluct=int(np.round(uniform(low=0, high=self.segment_fluctuation)))
            newind1 = ind1[0:int((len(ind1)/2)+fluct)] + ind2[(int((len(ind2)/2)+fluct)):len(ind2)+1]
            newind2 = ind2[0:int((len(ind2)/2)+fluct)] + ind1[(int((len(ind1)/2)+fluct)):len(ind1)+1]        
        
        if crossing_over_type == 'mask':
            
            newind1 = ind1
            newind2 = ind2
            
            mask = [int(np.round(uniform(low=0, high=1))) for x in range(self.vector_size)]
            #positive = [x for x in range(len(mask)) if mask[x] == 1]
            #negative = [x for x in range(len(mask)) if mask[x] != 1]
            
            newind1 = ind1
            newind1 = [ind1[x] if mask[x] == 1 else ind2[x] for x in range(len(ind2))]
            
            newind2 = ind2
            newind1 = [ind1[x] if mask[x] != 1 else ind2[x] for x in range(len(ind2))]
            
            
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
        
    def breed(self, population, scores):
        if self.selection_method == 'tournament':
            
            # Create pandas dict to facilitate sorting
            pop_scores = {'indv': population , 'score': scores}
            init_pop = pd.DataFrame(pop_scores)
            
            if self.opt_direction == 'up':
                init_pop = init_pop.sort_values('score', ascending=False)
            
            if self.opt_direction == 'down':
                init_pop = init_pop.sort_values('score', ascending=True)
        
            
            ### Initiate offspring dataframe
            #offspring = {'indv': [] , 'score': []}
            #offspring = pd.DataFrame(offspring)
            offspring = []
            
            #### number of tournament rounds are population/2 
            for rounds in range(int(len(population)/2)):
                to_select = init_pop.iloc[sample(range(len(population)), 2)]               
                if self.opt_direction == 'up':
                    to_select = to_select.sort_values('score', ascending=False)
                
                if self.opt_direction == 'down':
                    to_select = to_select.sort_values('score', ascending=True)
    
                
                co = self.crossing_over(to_select.iloc[0][0], to_select.iloc[1][0], self.crossing_over_type)
                for newindv in co:
                    offspring.append(newindv)
                    
            
            for trymut in range(len(offspring)):
                if uniform(low=0, high=1) < self.mutation_rate:
                    offspring[trymut] = self.mutate(offspring[trymut])
                        
            offspring_scores = self.calculate_scores(offspring)
            
            
            
            
            whole_pop = {'indv': population+offspring , 'score': scores+offspring_scores}
            whole_pop = pd.DataFrame(whole_pop)
            
            if self.opt_direction == 'up':
                whole_pop = whole_pop.sort_values('score', ascending=False)
            
            if self.opt_direction == 'down':
                whole_pop = whole_pop.sort_values('score', ascending=True)
                
                
            whole_pop = whole_pop[0:self.pop_size]

            return whole_pop['indv'].to_list(), whole_pop['score'].to_list()
            
                
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

        
        for t in range(self.n_cycles):
            print('Running round '+str(t))
            
            new_pop = self.breed(self.population, self.scores)
            
            self.population = new_pop[0]
            self.scores = new_pop[1]


            #### Tracking variables            
            self.pop_history.append(self.population)
            self.t=t
            print(t)
            self.best.append(self.scores[np.argmax(self.scores)])
            self.best_ind.append(self.population[np.argmax(self.scores)])
            self.pop_history.append(self.population)
            self.score_history.append(self.scores)
            
            self.finish_time = datetime.now()


    
    def execute(self):
        self.opt_cycle()
        
        

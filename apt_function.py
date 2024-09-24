#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Feb 12 19:52:36 2023

@author: lucas
"""

import pandas as pd
import subprocess
#Core Includes
import csv

import time
import os
from numpy.random import uniform
from random import sample
import numpy as np
import multiprocessing
from tqdm import tqdm

#### ESM stuff
import torch
import esm
import random
import math

from collections import Counter


# Load ESM-2 model
model, alphabet = esm.pretrained.esm2_t33_650M_UR50D()
batch_converter = alphabet.get_batch_converter()
model.eval()  # disables dropout for deterministic results

### END of ESM stuff


def Generate_random_population(starting_pose, pop_size, fixed_residues_list, chains):
    """
    Generate a random population of protein structures for optimization.

    Parameters:
    - starting_pose: PyRosetta Pose object representing the initial protein structure
    - pop_size: Number of individuals in the population
    - fixed_residues_list: List of residue indices to be fixed during optimization
    - chain: Chain identifier for fixed residues

    Returns:
    A tuple containing the initial population and the list of fixed residue indices.
    """    
    starting_pose_seq = [x for x in starting_pose.sequence()]
    gene_values = ['A','C','D','E','F','G','H','I','K','L','M','N','P','Q','R','S','T','V','W','Y']
    scorefxn(starting_pose)
    vector_size = len(starting_pose_seq)

    positions_starting_pose = list(range(0,vector_size))
    #lista_fixed_rosetta = PDB_to_Pose(starting_pose, fixed_residues_list, chains)
    lista_fixed_rosetta = fixed_residues_list
    #### Residues to be fixed during optimization and new population generation
    result_list = [position for position in positions_starting_pose if position not in lista_fixed_rosetta]

    #### Residues to be fixed during optimization and new population generation
    new_indiv = [[gene_values[int(np.round(uniform(low=0, high=len(gene_values)-1)))] for i in range(vector_size)] for indiv in range(pop_size)]

    #### Fixing the fixed positions to the new inds
    for individual in new_indiv:
        for i in result_list:
            individual[i-1] = starting_pose_seq[i-1]

    init_popu = list(new_indiv)
    
    return init_popu, result_list, lista_fixed_rosetta
    
    
def apt_benchmark(seq):
    """
    Benchmark a sequence for a specific trait.

    Parameters:
    - seq: Sequence to be benchmarked

    Returns:
    A benchmark score indicating the strength of the trait in the sequence.
    """ 
    score=0
    for i in range(len(seq)):
        if seq[i] == 'A':
            score = score+1
    
    ## add pack_relax?
    return score

def apt_rosetta(seq, starting_pose, scorefxn, dg_method, index_ind, index_cycle):
    """
    Perform threading optimization for a sequence.

    Parameters:
    - seq: Target sequence for threading
    - starting_pose: PyRosetta Pose object representing the initial protein structure
    - scorefxn: Score function to evaluate the energy of the structure
    - returning_val: List to store the result

    Returns:
    None (Result is stored in the returning_val list).
    """ 
    ###define starting pose outside of the function
    print("APT THREAD CHAMADA COM SUCESSO/nAPT THREAD CHAMADA COM SUCESSO/nAPT THREAD CHAMADA COM SUCESSO/nAPT THREAD CHAMADA COM SUCESSO/nAPT THREAD CHAMADA COM SUCESSO/n")
    scorefxn = pyrosetta.create_score_function("ref2015_cart.wts")
    
    resids, index = Get_residues_from_pose(pose = starting_pose)

    indexes = index
    resids_from_seq = resids
    
    to_mutate = Compare_sequences(resids_from_seq, seq, indexes)
    
    new_pose = starting_pose.clone()  
    for index in to_mutate:
        new_pose = mutate_repack(starting_pose = new_pose, posi = index, amino = to_mutate[index], scorefxn = scorefxn)
    #new_pose = pack_relax(starting_pose = new_pose, scorefxn = scorefxn, times_to_relax = 1)
    new_pose.dump_pdb(f"PDBs/{index_ind}_{index_cycle}.pdb")
    if dg_method == "bind":
        score = Dg_bind(new_pose, "A_D", scorefxn)
        data = pd.DataFrame({'Sequence': [new_pose.sequence()],'dG': [score]})
        data.to_csv(f'temp_{index_ind}.csv')
        return score
    if dg_method == "fold":
        data = pd.DataFrame({'Sequence': [new_pose.sequence()],'dG': [scorefxn(new_pose)]})
        data.to_csv(f'temp_{index_ind}.csv')
        return scorefxn(new_pose)

def apt(seq, starting_pose, scorefxn, dg_method, index_ind, index_cycle):
    """
    Perform an optimization process on a protein structure using an adaptive protein threading approach.

    Parameters:
    - seq: Target sequence to optimize towards
    - starting_pose: PyRosetta Pose object representing the initial protein structure
    - scorefxn: Score function to evaluate the energy of the structure
    - dg_method: Method to calculate the energy difference (dG)
    - index_ind: Index identifier
    - index_cycle: Index cycle

    Returns:
    The optimized score or dG value for the protein structure.
    """
    ###define starting pose outside of the function
    scorefxn = pyrosetta.create_score_function("ref2015_cart.wts")
    
    resids, index = Get_residues_from_pose(pose = starting_pose)

    indexes = index
    resids_from_seq = resids
    
    to_mutate = Compare_sequences(resids_from_seq, seq, indexes)
    
    new_pose = starting_pose.clone()  
    for index in to_mutate:
        new_pose = mutate_repack(starting_pose = new_pose, posi = index, amino = to_mutate[index], scorefxn = scorefxn)
    pack_relax(starting_pose = new_pose, scorefxn = scorefxn, times_to_relax = 1)
    new_pose.dump_pdb(f"PDBs/{index_ind}_{index_cycle}.pdb")
    if dg_method == "bind":
        
        score = Dg_bind(new_pose, "A_D", scorefxn)
        return score
    if dg_method == "fold":
        return scorefxn(new_pose)
    
    ## add pack_relax?
    #return score

def apt_thread(seq, starting_pose, scorefxn, dg_method, index_ind, index_cycle):
    """
    Perform threading optimization for a sequence.

    Parameters:
    - seq: Target sequence for threading
    - starting_pose: PyRosetta Pose object representing the initial protein structure
    - scorefxn: Score function to evaluate the energy of the structure
    - returning_val: List to store the result

    Returns:
    None (Result is stored in the returning_val list).
    """ 
    ###define starting pose outside of the function
    scorefxn = pyrosetta.create_score_function("ref2015_cart.wts")
    
    resids, index = Get_residues_from_pose(pose = starting_pose)

    indexes = index
    resids_from_seq = resids
    
    to_mutate = Compare_sequences(resids_from_seq, seq, indexes)
    
    new_pose = starting_pose.clone()  
    for index in to_mutate:
        new_pose = mutate_repack(starting_pose = new_pose, posi = index, amino = to_mutate[index], scorefxn = scorefxn)
    #new_pose = pack_relax(starting_pose = new_pose, scorefxn = scorefxn, times_to_relax = 1)
    new_pose.dump_pdb(f"PDBs/{index_ind}_{index_cycle}.pdb")
    if dg_method == "bind":
        #### Trocar apenas essa linha pelo calculo de dG usando o pbee
        command = f"python3 pbee.py --ipdb PDBs/{index_ind}_{index_cycle}.pdb --partner1 CD --partner2 A --odir . --force_mode"
        print(command)
        subprocess.run(command, stdout=subprocess.PIPE, shell=True)
        temp_df = pd.read_csv(f"outputs_pbee/{index_ind}_{index_cycle}/dG_pred.csv")
        score = temp_df["dG_pred"][0]
        data = pd.DataFrame({'Sequence': [new_pose.sequence()],'dG': [score]})
        data.to_csv(f'temp_{index_ind}.csv')
        return score
    if dg_method == "fold":
        data = pd.DataFrame({'Sequence': [new_pose.sequence()],'dG': [scorefxn(new_pose)]})
        data.to_csv(f'temp_{index_ind}.csv')
        return scorefxn(new_pose)

def unbind(pose, partners, scorefxn):
    """
    Simulate unbinding by applying a rigid body translation to a dummy pose and performing FastRelax.

    Parameters:
    - pose: PyRosetta Pose object representing the complex
    - partners: Vector1 specifying the interacting partners
    - scorefxn: Score function to evaluate the energy of the structure

    Returns:
    A tuple containing two Pose objects, one for the unbound state and one for the bound state.
    """
    #### Generates dummy pose to maintain original pose
    pose_dummy = pose.clone()
    pose_binded = pose.clone()
    STEP_SIZE = 100
    JUMP = 1
    docking.setup_foldtree(pose_dummy, partners, Vector1([-1,-1,-1]))
    trans_mover = rigid.RigidBodyTransMover(pose_dummy,JUMP)
    trans_mover.step_size(STEP_SIZE)
    trans_mover.apply(pose_dummy)
    pack_relax(pose_dummy, scorefxn, times_to_relax = 1)
    #### Return a tuple containing:
    #### Pose binded = [0] | Pose separated = [1]
    return pose_binded , pose_dummy

def dG_v2_0(pose_Sep, pose_bind, scorefxn):
    """
    Calculate the binding energy difference between the unbound and bound states.

    Parameters:
    - pose_Sep: Pose object representing the unbound state
    - pose_bind: Pose object representing the bound state
    - scorefxn: Score function to evaluate the energy of the structure

    Returns:
    The binding energy difference (dG).
    """
    bound_score = mc(pose_bind, scorefxn, 1000, 30)
    unbound_score = mc(pose_Sep, scorefxn, 1000, 30)
    dG = bound_score - unbound_score
    return dG

def Dg_bind(pose, partners, scorefxn):
    """
    Calculate the binding energy difference for a given complex.

    Parameters:
    - pose: PyRosetta Pose object representing the complex
    - partners: Vector1 specifying the interacting partners
    - scorefxn: Score function to evaluate the energy of the structure

    Returns:
    The binding energy difference (dG) for the complex.
    """
    pose_dummy = pose.clone()
    unbinded_dummy = unbind(pose_dummy, partners, scorefxn)

    return (dG_v2_0(unbinded_dummy[1], unbinded_dummy[0], scorefxn))

def Compare_sequences(before_seq, after_seq, indexes):
    """
    Compare two sequences and identify mutations. Print new mutations.

    Parameters:
    - before_seq: Original sequence
    - after_seq: Mutated sequence
    - indexes: Residue indexes

    Returns:
    A dictionary containing the mutated residues.
    """
    wt = before_seq
    mut = after_seq
    
    mutation = dict()
    for index, (res1, res2) in enumerate(zip(wt, mut)):
        if res1 != res2:
            mutation[indexes[index]] = res2
            print(f"New mutation \U0001f600 : {res1}{indexes[index]}{res2}")
    return mutation
def Get_residues_from_pose(pose):
    """
    Get the sequence and residue numbers for a specific chain in a given pose.

    Parameters:
    - pose: PyRosetta Pose object
    - chain: Chain identifier (e.g., 'A')

    Returns:
    A tuple containing the chain sequence and a list of residue numbers.
    """
    
    residue_numbers = [residue for residue in range(1, pose.size() + 1)]
    sequence = ''.join([pose.residue(residue).name1() for residue in residue_numbers])

    # residue_numbers = [residue for residue in range(1, pose.size() + 1) if pose.pdb_info().chain(residue) == chain]
    # chain_sequence = ''.join([pose.residue(residue).name1() for residue in residue_numbers])
    
    return sequence, residue_numbers

def correct_multi_input(population):
    correct = []
    for inds in population:
        correct.append(''.join(inds))
        
    return correct
        
def Run_all_batches(pose, list_seq,dg_method, batch_indexes, index_cycle):
    pose_init = pose.clone()
    scorefxn = pyrosetta.create_score_function("ref2015_cart.wts")
    
    processes = []

    for x in tqdm(range(len(batch_indexes)), desc="Processing list of sequences"):
        p = multiprocessing.Process(target=apt_thread,args=[list_seq[x], pose_init, scorefxn, dg_method, batch_indexes[x], index_cycle])

        p.start()
        processes.append(p)
        
    for p in processes:
        p.join()
        
def remove_files(file_list):
    for file_name in file_list:
        try:
            os.remove(file_name)
            print(f"File '{file_name}' removed successfully.")
        except OSError as e:
            print(f"Error removing file '{file_name}': {e}")


def check_files_in_directory(directory, file_list):
    """
    Check if all files in the given list exist in the specified directory.

    Args:
    directory (str): The directory to check for the files.
    file_list (list): A list of filenames to check for.

    Returns:
    bool: True if all files exist in the directory, False otherwise.
    """
    while True:
        # Check if the directory exists
        if not os.path.isdir(directory):
            print(f"Error: Directory '{directory}' does not exist.")
            time.sleep(10)  # Wait for 10 seconds before retrying
            continue

        # Check if all files exist in the directory
        missing_files = []
        for filename in file_list:
            file_path = os.path.join(directory, filename)
            if not os.path.isfile(file_path):
                missing_files.append(filename)
        
        if not missing_files:
            # All files exist in the directory
            to_save = retrieve_data(directory, file_list)
            remove_files(file_list)
            return to_save

        else:
            print(f"Not all files exist in directory '{directory}'. Missing files: {missing_files}")
            time.sleep(10)  # Wait for 10 seconds before retrying

def retrieve_data(directory, file_list):
    sequence = []
    dG = []
    for file in file_list:
        print(f"{file}")
        df_temp = pd.read_csv(f"{file}")
        sequence.append(df_temp.iloc[0,1])
        dG.append(df_temp.iloc[0,2])
    return sequence, dG



def batchs_to_run(pose, sequences, dg_method, batch_size, index_cycle):
    list_seq_final = []
    list_dg_final = []

    remaining_sequences = sequences[:]# Make a copy of the sequences
    batch_index = 0  # Initialize batch index
    while remaining_sequences:
        lista_seq = []
        lista_dg = []
        # Take a batch of sequences from the remaining ones
        batch = remaining_sequences[:batch_size]
        #list_files_temp = [f"temp_{i}.csv" for i in batch]
        list_files_temp = [f"temp_{batch_index * batch_size + i}.csv" for i in range(len(batch))]  # Adjusted line
        batch_indexes = [batch_index * batch_size + i for i in range(len(batch))] 
        remaining_sequences = remaining_sequences[batch_size:]  # Update remaining sequences
        # Process the current batch (here, just print the sequences)
        print("Processing batch:")
        Run_all_batches(pose, batch,dg_method, batch_indexes, index_cycle)
        ### Check it all files in list temp are in dir, if true, get sequences and dG from all 
        seqs, dgs = check_files_in_directory(".", list_files_temp)
        lista_seq.append(seqs)
        lista_dg.append(dgs)
        

        # Append the current batch size to the list
        list_seq_final += lista_seq
        list_dg_final += lista_dg
        time.sleep(5)
        batch_index += 1  # Increment batch index
    seqs_finais_ficou_de_verdade = [i for lista in list_seq_final for i in lista]
    dgs_finais_ficou_de_verdade = [i for lista in list_dg_final for i in lista]
    return dgs_finais_ficou_de_verdade

def mc(pose_ref, scorefxn, en_size, en_runs):
    
    
    logging.basicConfig(filename='mc.log', level=logging.DEBUG, format='%(asctime)s - %(levelname)s - %(message)s')

    logging.info("Starting Monte Carlo simulation")
    en_size=en_size ## ensemble_size
    kT=1.0
    n_moves=2 ### number of shear mover movements
    anglemax_ss=1 ### anglemax of shearmove in Alpha helix and beta-sheets
    anglemax=1 ### anglemax of shearmove in coil regions

    structure = pose_ref.clone()
    score = []
    score_final = []
    movemap = MoveMap()
    
    to_move = list(range(1, len(pose_ref.sequence())))
    
    move_map_list=[]
    for posi in to_move:
        mut_posi = pyrosetta.rosetta.core.select.residue_selector.ResidueIndexSelector()
        mut_posi.set_index(posi)
        
        # Select Neighbor Position
        nbr_selector = pyrosetta.rosetta.core.select.residue_selector.NeighborhoodResidueSelector()
        nbr_selector.set_focus_selector(mut_posi)
        nbr_selector.set_include_focus_in_subset(True)
        bla=pyrosetta.rosetta.core.select.get_residues_from_subset(nbr_selector.apply(pose_ref))
        for xres in bla:
            move_map_list.append(xres)

    
    movemap = MoveMap()
    for x in move_map_list:    
        movemap.set_bb(x,True)
        
    
    minmover = MinMover()
    minmover.movemap(movemap)
    minmover.score_function(scorefxn)
    
    
    ###################
    ## Version 4 changes
    ##################
    to_pack = standard_packer_task(structure)
    to_pack.restrict_to_repacking()    # prevents design, packing only
    to_pack.or_include_current(True)    # considers the original sidechains

    
    to_pack.temporarily_fix_everything()
    
    for xres in bla: 
        to_pack.temporarily_set_pack_residue(xres, True)
    
    packmover = PackRotamersMover(scorefxn, to_pack)
    
    
    shearmover = ShearMover(movemap, kT, n_moves) ## move phi e psi do mesmo residuo
    
    
    shearmover.angle_max('H', anglemax_ss) ##helices
    shearmover.angle_max('E', anglemax_ss) ## strands
    shearmover.angle_max('L', anglemax) ## loops
    
    
    
    combined_mover = SequenceMover()
    combined_mover.add_mover(shearmover)
    combined_mover.add_mover(packmover)



    before_pose = pose_ref.clone()
    for ensembles in range(en_runs):
        logging.info(f"Starting Monte Carlo simulation: run {ensembles}")
        for en in range(en_size):
            
            after_pose = before_pose.clone()
            
            #### APPLY
            combined_mover.apply(after_pose)
            #structure.dump_pdb("./ensembles/rbd_shear_"+str(en)+"_"+name_debug+".pdb")
    
           
            # The task factory accepts all the task operations
            
            #after_pose.dump_pdb("Dumps_ensemble/controle_mc_"+str(en)+".pdb")
            score.append(scorefxn(after_pose))
                
            before_pose = decision(before_pose, after_pose, scorefxn)
            
            
            #structure.dump_pdb("./ensembles/rbd_decision_"+str(en)+"_"+name_debug+".pdb")
    
            
            ###dump for debug    
           #before_pose.dump_pdb("Dumps_ensemble/teste_"+str(en)+".pdb")
        
        
    
            score_final.append(scorefxn(before_pose))
        
                
    
        
    return np.mean(score_final)

def decision(before_pose, after_pose, scorefxn):
    ### BEGIN SOLUTION
    E = scorefxn(after_pose) - scorefxn(before_pose)
    if E < 0:
        return after_pose
    elif random.uniform(0, 1) >= math.exp(-E/1):
        return before_pose
    else:
        return after_pose
    ### END SOLUTION

def penalize_contiguous_repetitions(sequence):
    """
    Penalizes contiguous repetitions in a sequence, normalized between 0 and 1.

    Args:
    sequence (list or str): The sequence to evaluate (could be a list of elements or a string).

    Returns:
    float: The normalized penalty for contiguous repetitions, between 0 and 1.
    """
    penalty = 0
    max_penalty = len(sequence) - 1  # Maximum penalty is if the entire sequence is a repetition

    # Iterate through the sequence and count contiguous repetitions
    for i in range(1, len(sequence)):
        if sequence[i] == sequence[i - 1]:
            penalty += 1

    # Normalize the penalty to be between 0 and 1
    normalized_penalty = penalty / max_penalty if max_penalty > 0 else 0

    return normalized_penalty

def apt_esm(seq, temperature=1):
    data = [
        ("protein1", seq)
    ]

    batch_labels, batch_strs, batch_tokens = batch_converter(data)
    batch_lens = (batch_tokens != alphabet.padding_idx).sum(1)

    # Predict token probabilities
    with torch.no_grad():
        token_probs = model(batch_tokens, repr_layers=[33])["logits"]

    # Apply temperature
    token_probs /= temperature

    softmax = torch.nn.Softmax(dim=-1)
    probabilities = softmax(token_probs)

    # Calculate the average probability of the tokens at each position, excluding start and end tokens
    average_probs = []
    for i in range(1, batch_lens[0] - 1):  # Exclude start and end tokens
        token_idx = batch_tokens[0, i].item()
        token_prob = probabilities[0, i, token_idx].item()
        average_probs.append(token_prob)

    avg_probability = sum(average_probs) / len(average_probs)

    # Apply repetition penalty if rep_penalty is True


    return avg_probability
                    


def apt_esm_penalty(seq, temperature=1):
    data = [
        ("protein1", seq)
    ]

    batch_labels, batch_strs, batch_tokens = batch_converter(data)
    batch_lens = (batch_tokens != alphabet.padding_idx).sum(1)

    # Predict token probabilities
    with torch.no_grad():
        token_probs = model(batch_tokens, repr_layers=[33])["logits"]

    # Apply temperature
    token_probs /= temperature

    softmax = torch.nn.Softmax(dim=-1)
    probabilities = softmax(token_probs)

    # Calculate the average probability of the tokens at each position, excluding start and end tokens
    average_probs = []
    for i in range(1, batch_lens[0] - 1):  # Exclude start and end tokens
        token_idx = batch_tokens[0, i].item()
        token_prob = probabilities[0, i, token_idx].item()
        average_probs.append(token_prob)

    avg_probability = sum(average_probs) / len(average_probs)
    
    # Apply repetition penalty if rep_penalty is True
    repetition_penalty = penalize_contiguous_repetitions(seq)
    avg_probability *= (1 - repetition_penalty)  # Apply penalty to average probability


    return avg_probability


def shannon_entropy(sequence):
    """Calculates the Shannon entropy of a sequence and normalizes it between 0 and 1."""
    counts = Counter(sequence)
    total = len(sequence)
    entropy = -sum((count / total) * np.log2(count / total) for count in counts.values())

    # Normalize entropy by the maximum possible entropy
    num_unique_tokens = len(set(sequence))
    max_entropy = np.log2(num_unique_tokens)

    # Avoid division by zero when there's only one unique token
    if max_entropy > 0:
        normalized_entropy = entropy / max_entropy
    else:
        normalized_entropy = 0

    return normalized_entropy

def apt_esm_shannon_penalty(seq, temperature=1):
    data = [
        ("protein1", seq)
    ]

    batch_labels, batch_strs, batch_tokens = batch_converter(data)
    batch_lens = (batch_tokens != alphabet.padding_idx).sum(1)

    # Predict token probabilities
    with torch.no_grad():
        token_probs = model(batch_tokens, repr_layers=[33])["logits"]

    # Apply temperature
    token_probs /= temperature

    softmax = torch.nn.Softmax(dim=-1)
    probabilities = softmax(token_probs)

    # Calculate the average probability of the tokens at each position, excluding start and end tokens
    average_probs = []
    for i in range(1, batch_lens[0] - 1):  # Exclude start and end tokens
        token_idx = batch_tokens[0, i].item()
        token_prob = probabilities[0, i, token_idx].item()
        average_probs.append(token_prob)

    avg_probability = sum(average_probs) / len(average_probs)

    # Calculate normalized Shannon entropy of the sequence
    normalized_entropy = shannon_entropy(seq)

    # Use normalized entropy directly to scale the probability
    adjusted_probability = avg_probability * normalized_entropy

    return adjusted_probability



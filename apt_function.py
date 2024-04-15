#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Feb 12 19:52:36 2023

@author: lucas
"""

from pyrosetta import *
from rosetta.core.pack.task import TaskFactory
from rosetta.core.pack.task import operation
from rosetta.core.kinematics import MoveMap
from rosetta.core.kinematics import FoldTree
from rosetta.core.pack.task import TaskFactory
from rosetta.core.pack.task import operation
from rosetta.core.simple_metrics import metrics
from rosetta.core.select import residue_selector as selections
from rosetta.core import select
from rosetta.core.select.movemap import *
from rosetta.protocols import minimization_packing as pack_min
from rosetta.protocols import relax as rel
from rosetta.protocols.antibody.residue_selector import CDRResidueSelector
from rosetta.protocols.antibody import *
from rosetta.protocols.loops import *
from rosetta.protocols.relax import FastRelax
from pyrosetta.rosetta.protocols.docking import setup_foldtree
from pyrosetta.rosetta.protocols import *
from rosetta.core.scoring.methods import EnergyMethodOptions
import pandas as pd
from Bio import SeqIO


pyrosetta.init()


scorefxn = pyrosetta.create_score_function("ref2015_cart.wts")

######### take sequence, mutate PDB, return score
def pack_relax(starting_pose, scorefxn):

    pose = starting_pose.clone()

    tf = TaskFactory()
    tf.push_back(operation.InitializeFromCommandline())
    tf.push_back(operation.RestrictToRepacking())

    # Set up a MoveMapFactory
    mmf = pyrosetta.rosetta.core.select.movemap.MoveMapFactory()
    mmf.all_bb(setting=True)
    mmf.all_bondangles(setting=True)
    mmf.all_bondlengths(setting=True)
    mmf.all_chi(setting=True)
    mmf.all_jumps(setting=True)
    mmf.set_cartesian(setting=True)


    fr = pyrosetta.rosetta.protocols.relax.FastRelax(scorefxn_in=scorefxn, standard_repeats=1)
    fr.cartesian(True)
    fr.set_task_factory(tf)
    fr.set_movemap_factory(mmf)
    fr.min_type("lbfgs_armijo_nonmonotone")
    fr.apply(pose)
    return pose

def mutate_repack(starting_pose, posi, amino, scorefxn):

    pose = starting_pose.clone()

     #Select position to mutate
    mut_posi = pyrosetta.rosetta.core.select.residue_selector.ResidueIndexSelector()
    mut_posi.set_index(posi)

    #Select neighbor positions
    nbr_selector = pyrosetta.rosetta.core.select.residue_selector.NeighborhoodResidueSelector()
    nbr_selector.set_focus_selector(mut_posi)
    nbr_selector.set_include_focus_in_subset(True)

    not_design = pyrosetta.rosetta.core.select.residue_selector.NotResidueSelector(mut_posi)

    tf = pyrosetta.rosetta.core.pack.task.TaskFactory()

    tf.push_back(pyrosetta.rosetta.core.pack.task.operation.InitializeFromCommandline())
    tf.push_back(pyrosetta.rosetta.core.pack.task.operation.IncludeCurrent())
    tf.push_back(pyrosetta.rosetta.core.pack.task.operation.NoRepackDisulfides())

    # Disable Packing
    prevent_repacking_rlt = pyrosetta.rosetta.core.pack.task.operation.PreventRepackingRLT()
    prevent_subset_repacking = pyrosetta.rosetta.core.pack.task.operation.OperateOnResidueSubset(prevent_repacking_rlt, nbr_selector, True )
    tf.push_back(prevent_subset_repacking)

    # Disable design
    tf.push_back(pyrosetta.rosetta.core.pack.task.operation.OperateOnResidueSubset(pyrosetta.rosetta.core.pack.task.operation.RestrictToRepackingRLT(),not_design))

    # Enable design
    aa_to_design = pyrosetta.rosetta.core.pack.task.operation.RestrictAbsentCanonicalAASRLT()
    aa_to_design.aas_to_keep(amino)
    tf.push_back(pyrosetta.rosetta.core.pack.task.operation.OperateOnResidueSubset(aa_to_design, mut_posi))

    # Create Packer
    packer = pyrosetta.rosetta.protocols.minimization_packing.PackRotamersMover(scorefxn)
    packer.task_factory(tf)
    packer.apply(pose)

    return pose

def Read_sequence(fasta):
    # Specify the path to your FASTA file
    fasta_file = fasta

    # Use SeqIO.parse to read the FASTA file
    sequences = SeqIO.parse(fasta_file, "fasta")

    seq_records = [record for record in sequences]

    seqs = [str(seq.seq) for seq in seq_records]


    sequences_listed = [[x for x in sequence] for sequence in seqs]

    return sequences_listed

def PDB_pose_dictionairy(pose_to):
    """
    Create a Pandas DataFrame representing the mapping of residue positions
    between PDB numbering and Pose numbering for a given PyRosetta pose.

    Parameters:
    - pose_to (PyRosetta pose): The PyRosetta pose object for which to create
                                the PDB to Pose numbering dictionary.

    Returns:
    - pandas.DataFrame: A DataFrame containing three columns - 'Chain',
                        'IndexPDB', and 'IndexPose' representing the mapping
                        of residue positions between PDB and Pose numbering.

    Example:
    >>> pose = ...  # Initialize your PyRosetta pose object
    >>> pdb_pose_dict = PDB_pose_dictionairy(pose)
    >>> print(pdb_pose_dict)
         Chain  IndexPDB  IndexPose
    0        A         1          1
    1        A         2          2
    ...      ...       ...        ...
    n        B       100         99
    """
    # List to store data during iteration
    vectorIndexChain = []
    vectorIndexPDB = []
    vectorIndexPose = []

    # Creating dictionary for residue position - PDB and Pose numbering
    for i in range(pose_to.total_residue()):
        vectorIndexChain.append(pose_to.pdb_info().chain(i + 1))
        vectorIndexPDB.append(pose_to.pdb_info().number(i + 1))
        vectorIndexPose.append(pose_to.pdb_info().pdb2pose(vectorIndexChain[i], vectorIndexPDB[i]))

    # Inserting values into a pandas dataframe
    df_pdbtopose_dictionary = {"Chain": vectorIndexChain,
                                "IndexPDB": vectorIndexPDB,
                                "IndexPose": vectorIndexPose}
    df_dictionary = pd.DataFrame(df_pdbtopose_dictionary)

    return df_dictionary

def PDB_to_Pose(pose, index_list, chain):

    df_indexes = PDB_pose_dictionairy(pose)
    lists = list(df_indexes[(df_indexes["IndexPDB"].isin(index_list)) & (df_indexes["Chain"] == chain)]["IndexPose"])
    return lists

def apt_benchmark(seq):
    score=0
    for i in range(len(seq)):
        if seq[i] == 'A':
            score = score+1

    ## add pack_relax?
    return score


def apt(seq, starting_pose, scorefxn, dg_method, index_ind, index_cycle):
        ###define starting pose outside of the function
        scorefxn = pyrosetta.create_score_function("ref2015_cart.wts")

        resids, index = Get_residues_from_pose(pose = starting_pose)

        indexes = index
        resids_from_seq = resids

        to_mutate = Compare_sequences(resids_from_seq, seq, indexes)

        new_pose = starting_pose.clone()
        for index in to_mutate:
            new_pose = mutate_repack(starting_pose = new_pose, posi = index, amino = to_mutate[index], scorefxn = scorefxn)
        pack_relax(starting_pose = new_pose, scorefxn = scorefxn)
        new_pose.dump_pdb(f"PDBs/{index_ind}_{index_cycle}.pdb")
        if dg_method == "bind":

            score = Dg_bind(new_pose, "A_D", scorefxn)
            return score
        if dg_method == "fold":
            return scorefxn(new_pose)

        ## add pack_relax?
        #return score

def apt_thread(seq, starting_pose, scorefxn, returning_val):

    ###define starting pose outside of the function

    for i in range(len(seq)):
        starting_pose = mutate_repack(starting_pose, posi=i+1, amino=seq[i])
        starting_pose = pack_relax(starting_pose, scorefxn)
        returning_val [0] = scorefxn(starting_pose)
    ## add pack_relax?
    #return score

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
    pack_relax(pose_dummy, scorefxn)
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
    bound_score = scorefxn(pose_bind)
    unbound_score = scorefxn(pose_Sep)
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
    #unbinded_dummy[1].dump_pdb("./Dumps/unbinded_teste.pdb")
    #unbinded_dummy[0].dump_pdb("./Dumps/binded_teste.pdb")
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

    return sequence, residue_numbers

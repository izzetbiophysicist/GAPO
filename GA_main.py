# GAprot.py - Versão Refatorada (Exemplo)

import os
import argparse

from genetic_algorithm_rosetta import genetic_algo, genetic_algo_sequence

from apt_function import *
import apt_function

from pyrosetta import *
from rosetta.core.pack.task import TaskFactory
from rosetta.core.pack.task import operation

import numpy as np
from numpy.random import uniform
from random import sample
import random
# ... outras importações

def main():
    parser = argparse.ArgumentParser(description="GAPO: Genetic Algorithm for Protein Optimization")
    subparsers = parser.add_subparsers(dest="command", required=True, help="Sub-command help")

    # --- Parser para modo "structure" ---
    parser_struct = subparsers.add_parser("structure", help="Structure-based optimization")
    parser_struct.add_argument("--pdb", type=str, required=True, help="Input PDB file")
    parser_struct.add_argument("--pop_size", type=int, default=50, help="Population size")
    parser_struct.add_argument("--cycles", type=int, default=50, help="Number of GA cycles")
    parser_struct.add_argument("--mutation_type", type=str, default= "esm", help="Mutation type during optimization")
    parser_struct.add_argument("--mutation_rate", type=float, default=0.9, help="Mutation rate")
    parser_struct.add_argument("--direction", type=str,default='down', help="Direction to optimize: 'up' or 'down'")
    parser_struct.add_argument("--apt_function", type=str, default='rosetta', choices=['rosetta', 'esm', 'esm_penalty', 'esm_shannon_penalty'], help="Aptitude function to use")
    parser_struct.add_argument("--residues_to_mut", type=int, nargs='+', required=True, help="List of residue indices (PDB numbering) to  mutate")
    parser_struct.add_argument("--temp", type=float, default= 1.5, help="ESM2 Temperature to control the randomness of the mutations")
    parser_struct.add_argument("--output_file", type=str, default="gapo_results", help="Output file name")
    parser_struct.add_argument("--cpus", type=int, default=1, help="Number of CPUs for parallel processing")

    # --- Parser para modo "sequence" ---
    parser_seq = subparsers.add_parser("sequence", help="Sequence-based optimization")
    parser_seq.add_argument("--seq", type=str, required=True, help="Input PDB file")
    parser_seq.add_argument("--cycles", type=int, default=50, help="Number of GA cycles")
    parser_seq.add_argument("--mutation_type", type=str, default= "esm", help="Mutation type during optimization")
    parser_seq.add_argument("--pop_size", type=int, default=50, help="Population size")
    parser_seq.add_argument("--mutation_rate", type=float, default=0.9, help="Mutation rate")
    parser_seq.add_argument("--direction", type=str, default='up', help="Direction to optimize: 'up' or 'down'")
    parser_seq.add_argument("--apt_function", type=str, default='esm', choices=['rosetta', 'esm', 'esm_penalty', 'esm_shannon_penalty'], help="Aptitude function to use")
    parser_seq.add_argument("--residues_to_mut", type=int, nargs='+', required=True, help="List of residue indices to  mutate")
    parser_seq.add_argument("--temp", type=float, default= 1.5, help="ESM2 Temperature to control the randomness of the mutations")
    parser_seq.add_argument("--output_file", type=str, default="gapo_results", help="Output file name")
    parser_seq.add_argument("--cpus", type=int, default=1, help="Number of CPUs for parallel processing")

    args = parser.parse_args()

    # --- Mapeamento de funções de aptidão ---
    apt_functions = {
        'rosetta': apt_rosetta,
        'esm': apt_esm,
        'esm_penalty': apt_esm_penalty,
        'esm_shannon_penalty': apt_esm_shannon_penalty
    }
    selected_apt_function = apt_functions[args.apt_function]

    gene_values = ['A','C','D','E','F','G','H','I','K','L','M','N','P','Q','R','S','T','V','W','Y']

    if args.command == "structure":
        pyrosetta.init(extra_options="\
    -mute core \
    -mute basic \
    ")
        
        starting_pose = pose_from_pdb(args.pdb)
        inds_sequence = list(range(1, len(starting_pose.sequence())+1))
        init_pop = apt_function.generate_population_esm(starting_pose.sequence(), args.residues_to_mut, population_size=args.pop_size, temperature=args.temp)
        fixed_residues = [i for i in inds_sequence if i not in args.residues_to_mut]

        GA = genetic_algo(
            pdb=args.pdb, 
            opt_direction=args.direction,
            initial_population=init_pop,
            gene_values=gene_values,
            gene_type='discrete',
            vector_size=len(starting_pose.sequence()),
            pop_size=len(init_pop),
            mutation_rate=args.mutation_rate,
            segment_fluctuation=0,
            selection_method='tournament',
            threads=False,
            convergence_threshold=0,
            tournament_cycles=int(np.round(len(init_pop)/4)),
            tournament_size=4,
            benchmark=False,
            lista_fixed=fixed_residues,
            crossing_over_type='mask',
            apt_function=selected_apt_function,
            n_cycles=args.cycles,
            file_name=args.output_file,
            cpus=args.cpus,
            mutation_type = args.mutation_type,
            # ...
        )
        GA.execute()

    elif args.command == "sequence":
        # Lógica para gerar população inicial de sequências
        starting_sequence = args.seq
        init_pop = apt_function.generate_population_esm(starting_sequence, args.residues_to_mut, population_size=args.pop_size, temperature=args.temp)
        inds_sequence = list(range(1, len(starting_sequence)+1))
        fixed_residues = [i for i in inds_sequence if i not in args.residues_to_mut]

        GA = genetic_algo_sequence(
            opt_direction=args.direction,
            initial_population=init_pop,
            gene_values=gene_values,
            gene_type='discrete',
            vector_size=len(starting_sequence),
            pop_size=len(init_pop),
            mutation_rate=args.mutation_rate,
            segment_fluctuation=0,
            selection_method='tournament',
            threads=False,
            convergence_threshold=0,
            tournament_cycles=int(np.round(len(init_pop)/4)),
            tournament_size=4,
            benchmark=False,
            crossing_over_type='mask',
            apt_function=selected_apt_function,
            n_cycles=args.cycles,
            lista_fixed=fixed_residues,
            file_name=args.output_file,
            mutation_type = args.mutation_type,

            # ...
        )
        GA.execute()

if __name__ == "__main__":
    main()

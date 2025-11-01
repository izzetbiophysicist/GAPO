# GAPO - Genetic Algorithm for Protein Optimization 🧬

<div align="center">
  <img src="./GAPO_logo.png" alt="GAPO Logo" width="400">
</div>

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
  - [Algorithm Parameters](#algorithm-parameters)
  - [Example 1: Structure-Based Optimization](#example-1-structure-based-optimization)
  - [Example 2: Sequence-Based Optimization](#example-2-sequence-based-optimization)
 -[Preprint DOI](https://www.biorxiv.org/content/10.1101/2025.10.29.685431v1) 
- [Contributing](#contributing)

---

## About The Project

**GAPO** is an _in silico_ genetic algorithm used to optimize proteins for a desired function, such as stability and affinity. The algorithm mimics the evolutionary process by recombining and adding mutations to the best sequences in order to generate a new population with higher diversity and optimized for the given objective function.

---

## Getting Started

Follow these steps to set up and run the project locally.

### Prerequisites

Before you begin, ensure you have **Conda** installed on your system.

- If you don't have Conda, follow the installation instructions on the official website: [**Anaconda Installation**](https://www.anaconda.com/download).

### Installation

1.  **Clone the Repository**
    ```sh
    git clone [https://github.com/izzetbiophysicist/prot_eng_GA.git](https://github.com/izzetbiophysicist/prot_eng_GA.git)
    cd prot_eng_GA
    ```

2.  **Create the Base Conda Environment**
    This command uses the `environment.yml` file to create a new environment named `gapo_env` with all the base dependencies.
    ```sh
    conda env create -f environment.yml
    ```

3.  **Activate the New Environment**
    ```sh
    conda activate gapo_env
    ```

4.  **Install PyTorch (⚠️ Crucial Step)**
    The `environment.yml` file **does not install PyTorch** to ensure you choose the correct version for your hardware. You must install it manually.

    * **🚀 For NVIDIA GPU Users (Highly Recommended):**
        Visit the **[Official PyTorch Website](https://pytorch.org/get-started/locally/)**. Select the settings that match your system (e.g., Conda, Python, your CUDA version) and run the generated command. It will look something like this:
        ```sh
        # This is an EXAMPLE command, get the correct one from the PyTorch website!
        conda install pytorch torchvision torchaudio pytorch-cuda=11.8 -c pytorch -c nvidia
        ```

    * **💻 For CPU-Only Users:**
        If you do not have a compatible GPU, install the CPU-only version of PyTorch with this command:
        ```sh
        conda install pytorch torchvision torchaudio cpuonly -c pytorch
        ```

5.  **Install PyRosetta**
    Finally, install PyRosetta using its dedicated installer.
    ```sh
    pip install pyrosetta-installer
    python -c 'import pyrosetta_installer; pyrosetta_installer.install_pyrosetta()'
    ```
    
---

## Usage

GAPO is run from the command line, specifying the optimization mode (`structure` or `sequence`) and the desired parameters.

### Algorithm Parameters 

Below are the command-line parameters for each mode of operation in GAPO.

#### Structure Mode Parameters

##### **Required Parameters**
| Parameter            | Description                                                 |
| :------------------- | :---------------------------------------------------------- |
| `--pdb`              | The input PDB file for the optimization.                    |
| `--residues_to_mute` | List of residue indices (PDB numbering) to be mutated.      |

##### **Optional Parameters**
| Parameter         | Description                                                    | Default Value |
| :---------------- | :------------------------------------------------------------- | :------------ |
| `--pop_size`      | Size of the population in each generation.                     | `50`          |
| `--cycles`        | Number of cycles (generations) for the genetic algorithm.      | `50`          |
| `--mutation_type` | Type of mutation to be used during optimization.               | `esm`         |
| `--mutation_rate` | The mutation rate applied to the population.                   | `0.9`         |
| `--direction`     | Optimization direction: `up` (maximize) or `down` (minimize).  | `down`        |
| `--apt_function`  | Aptitude function to be used.                                  | `rosetta`     |
| `--temp`          | ESM2 temperature to control the randomness of mutations.       | `1.5`         |
| `--output_file`   | Base name for the output file.                                 | `gapo_results`|
| `--cpus`          | Number of CPUs to use for parallel processing.                 | `1`           |

---

#### Sequence Mode Parameters

##### **Required Parameters**
| Parameter            | Description                                                      |
| :------------------- | :------------------------------------------------------------- |
| `--seq`              | The initial amino acid sequence for optimization.              |
| `--residues_to_mute` | List of indices for the residues in the sequence to be mutated.|

##### **Optional Parameters**
| Parameter         | Description                                                    | Default Value |
| :---------------- | :------------------------------------------------------------- | :------------ |
| `--pop_size`      | Size of the population in each generation.                     | `50`          |
| `--cycles`        | Number of cycles (generations) for the genetic algorithm.      | `50`          |
| `--mutation_type` | Type of mutation to be used during optimization.               | `esm`         |
| `--mutation_rate` | The mutation rate applied to the population.                   | `0.9`         |
| `--direction`     | Optimization direction: `up` (maximize) or `down` (minimize).  | `up`          |
| `--apt_function`  | Aptitude function to be used.                                  | `esm`         |
| `--temp`          | ESM2 temperature to control the randomness of mutations.       | `1.5`         |
| `--output_file`   | Base name for the output file.                                 | `gapo_results`|
 
### Example 1: Structure-Based Optimization

This example optimizes the CDRs of an scFv based on its PDB structure, using the Rosetta score as the objective function.

```bash
python GA_main.py structure \
    --pdb ab_trimed_relax.pdb \
    --residues_to_mute 62 63 64 65 66 67 68 69 70 71 72 88 89 90 91 92 93 94 127 128 129 130 131 132 133 134 135 186 187 188 189 190 191 192 212 213 214 215 216 257 258 259 260 261 262 263 264 265 266 267 268 269 \
    --apt_function rosetta \
    --pop_size 50 \
    --cycles 10 \
    --opt_direction down \
    --output_file rosetta_run_01
```
### Example 2: Sequence-Based Optimization

This example takes an initial peptide sequence and evolves it to maximize its likelihood according to the ESM-2 model, mutating only the core region.

```bash
python GA_main.py sequence \
    --seq "RKVCNGIGIGEFKDSLSINATNIKHFKNCTSISGDLHILPVAFRGDSFTHTPPLDPQELDILKTVKEITGFLLIQAW" \
    --residues_to_mute 3 4 5 6 7 8 9 \
    --apt_function esm \
    --pop_size 100 \
    --cycles 20 \
    --direction up \
    --output_file esm_run_peptide
```
## Contributing

Contributions are what make the open-source community such an amazing place to learn, inspire, and create. Any contributions you make are **greatly appreciated**.

1.  Fork the Project
2.  Create your Feature Branch (`git checkout -b feature/AmazingFeature`)
3.  Commit your Changes (`git commit -m 'Add some AmazingFeature'`)
4.  Push to the Branch (`git push origin feature/AmazingFeature`)
5.  Open a Pull Request

---


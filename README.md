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
  - [Example 1: Structure-Based Optimization](#example-1-structure-based-optimization)
  - [Example 2: Sequence-Based Optimization](#example-2-sequence-based-optimization)
  - [Algorithm Parameters](#algorithm-parameters)
- [Contributing](#contributing)
- [License](#license)
- [Contact](#contact)

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

1.  **Clone the repository:**
    ```sh
    git clone [https://github.com/izzetbiophysicist/prot_eng_GA.git](https://github.com/izzetbiophysicist/prot_eng_GA.git)
    cd prot_eng_GA
    ```

2.  **Create the Conda environment:**
    This command will create an environment named `gapo_env` with all the necessary dependencies listed in the `environment.yml` file.
    ```sh
    conda env create -f environment.yml
    ```

3.  **Activate the new environment:**
    ```sh
    conda activate gapo_env
    ```

4.  **Install PyRosetta:**
    PyRosetta requires a license and is installed separately. The `pyrosetta-installer` simplifies this process.
    ```sh
    pip install pyrosetta-installer
    python -c 'import pyrosetta_installer; pyrosetta_installer.install_pyrosetta()'
    ```
    
    > **Note on PyTorch:** The `environment.yml` file installs a **CPU-only** version of PyTorch. If you have an NVIDIA GPU, it is highly recommended to visit the [**PyTorch installation page**](https://pytorch.org/get-started/locally/) to get the optimized installation command for your CUDA version and run it after activating the environment.

---

## Usage

GAPO is run from the command line, specifying the optimization mode (`structure` or `sequence`) and the desired parameters.

### Example 1: Structure-Based Optimization

This example optimizes the CDRs of an scFv based on its PDB structure, using the Rosetta score as the objective function.

```bash
python GAprot.py structure \
    --pdb inputs/CD19_scFv_relax.pdb \
    --fixed_residues 62 63 64 65 66 67 68 69 70 71 72 88 89 90 91 92 93 94 127 128 129 130 131 132 133 134 135 186 187 188 189 190 191 192 212 213 214 215 216 257 258 259 260 261 262 263 264 265 266 267 268 269 \
    --chains C D \
    --apt_function rosetta \
    --pop_size 50 \
    --cycles 10 \
    --opt_direction down \
    --output_file results/rosetta_run_01.csv

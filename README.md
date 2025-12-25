# Alzheimer's Protein Docking Simulation

This project is a high-throughput screening and docking simulation model designed to investigate the binding interactions between the Alzheimer's Amyloid Beta (Aβ-42) protein and various potential therapeutic molecules.

## Project Goal
To build a predictive model, easy for visualization, and grounded in real databases to filter and select known molecules that have the highest probability to bind with B-amyloid protein and prevent accumulation.

## Methodology

### 1. Docking Simulation
We use **AutoDock Vina**, a state-of-the-art docking software, to calculate the binding affinity (kcal/mol) between the protein and ligands.
*   **Target Protein**: Amyloid Beta (Aβ-42) (PDB ID: 1IYT).
*   **Target Site**: 'KLVFF' region (residues 16-20), known for its role in aggregation.

### 2. Molecule Categories
We compare three groups of molecules:
*   **FDA Approved Drugs**: e.g., Donepezil, Galantamine.
*   **Flavonoids**: Natural compounds found in fruits/vegetables (e.g., Curcumin, Quercetin).
*   **Known Ligands**: Compounds known to bind Aβ (e.g., Congo Red, Thioflavin T).

### 3. Evolutionary Optimization (AI)
We implemented a Genetic Algorithm (Evolutionary "Back Propagation") to optimize ligand structures. The system:
1.  Takes a parent molecule.
2.  Mutates its structure (adds atoms, substitutes elements).
3.  Simulates docking for the new structures.
4.  Selects the best binders for the next generation.

## Installation & Usage

### Prerequisites
*   Python 3.8+
*   AutoDock Vina (Binary included or installed in path)

### Install Dependencies
```bash
pip install streamlit rdkit meeko py3Dmol stmol pandas numpy scipy
```

### Run the Simulation
```bash
streamlit run app.py
```

## Structure
*   `app.py`: Main Streamlit web application.
*   `docking.py`: Core docking engine handling Vina execution.
*   `evolution.py`: Genetic algorithm for ligand optimization.
*   `prep_protein.py`: Utility to prepare the PDB structure.
*   `data/`: Contains molecule database and protein structures.

## Sources
*   **Protein**: [RCSB PDB - 1IYT](https://www.rcsb.org/structure/1IYT)
*   **Molecules**: Structures from [PubChem](https://pubchem.ncbi.nlm.nih.gov/).
*   **Docking**: [AutoDock Vina](https://vina.scripps.edu/).
*   **Ligand Prep**: [Meeko](https://github.com/forlilab/meeko).

## Credits
Developed for Science Fair Project.

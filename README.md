# Tau Protein Docking Simulation

This project is a high-throughput screening and docking simulation model designed to investigate the binding interactions between the **Tau Protein Fibril** (PDB: 5O3L) and various potential therapeutic molecules.

## Project Goal
To build a predictive model, easy for visualization, and grounded in real databases to filter and select known molecules that have the highest probability to bind with Tau protein aggregates and prevent neurodegeneration.

## Methodology

### 1. Docking Simulation
We use **AutoDock Vina**, a state-of-the-art docking software, to calculate the binding affinity (kcal/mol) between the protein and ligands.
*   **Target Protein**: Tau Fibril (Alzheimer's Brain) (PDB ID: 5O3L).
*   **Target Site**: **VQIVYK (PHF6)** motif (Residues 306-311), a critical region for filament aggregation.

### 2. Molecule Categories
We compare three groups of molecules:
*   **FDA Approved Drugs**: e.g., Methylene Blue (TRx0237), Epothilone D.
*   **Flavonoids**: Natural compounds (e.g., Curcumin, EGCG).
*   **Known Ligands**: PET tracers and inhibitors (e.g., PBB3, T807/AV-1451).

### 3. Evolutionary Optimization (AI)
We implemented a Genetic Algorithm (Evolutionary "Back Propagation") to optimize ligand structures. The system:
1.  Takes a parent molecule.
2.  Mutates its structure (adds atoms, substitutes elements).
3.  Simulates docking for the new structures.
4.  Selects the best binders for the next generation.

## Installation & Usage

### Prerequisites
*   Python 3.8+

### 1. Install AutoDock Vina
Since `pip install vina` can be difficult on MacOS, we recommend downloading the executable directly.

1.  **Download Vina**: Go to the [AutoDock Vina Releases Page](https://github.com/ccsb-scripps/AutoDock-Vina/releases).
2.  **Select Version**: Download the version for your OS (e.g., `vina_1.2.5_mac_x86_64` for Intel Mac, or look for an arm64 version if available, though x86_64 works on M1/M2/M3 via Rosetta).
3.  **Place in Folder**:
    *   Rename the downloaded file to just `vina` (or `vina_mac`).
    *   Move this file into the root folder of this project (where `app.py` is).
    *   **Permission (Mac)**: You might need to run `chmod +x vina` in the terminal to make it executable. If MacOS blocks it, right-click and "Open" to allow it once.

### 2. Install Python Dependencies
```bash
pip install streamlit rdkit meeko py3Dmol stmol pandas numpy scipy
```

> **Note**: During installation of `meeko`, you might see a `SyntaxError` related to `compute_water_map.py`. This is due to a legacy file included in the package and **can be safely ignored**.

### 3. Run the Simulation
```bash
streamlit run app.py
```

## Structure
*   `app.py`: Main Streamlit web application.
*   `docking.py`: Core docking engine handling Vina execution.
*   `evolution.py`: Genetic algorithm for ligand optimization.
*   `prep_protein.py`: Utility to prepare the PDB structure.
*   `data/`: Contains molecule database and protein structures (5O3L).

## Sources
*   **Protein**: [RCSB PDB - 5O3L](https://www.rcsb.org/structure/5O3L)
*   **Molecules**: Structures from [PubChem](https://pubchem.ncbi.nlm.nih.gov/).
*   **Docking**: [AutoDock Vina](https://vina.scripps.edu/).
*   **Ligand Prep**: [Meeko](https://github.com/forlilab/meeko).

## Credits
Developed for Science Fair Project.

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

### 1. Install AutoDock Vina (The Engine)
Since `pip install vina` fails on many Macs, **do not use pip for Vina**.

1.  **Download**: Go to [AutoDock Vina Releases](https://github.com/ccsb-scripps/AutoDock-Vina/releases).
2.  **Select**: Download `vina_1.2.5_mac_x86_64` (or the version for your OS).
3.  **Setup**:
    *   Rename the downloaded file to `vina`.
    *   Move it into this project folder (same folder as `app.py`).
    *   **Mac Permission**: Open terminal in this folder and run `chmod +x vina`. (If it's blocked by security, right-click the file in Finder and select "Open" to allow it).

### 2. Install Python Dependencies
Run these commands one by one to ensure everything installs correctly:

1.  **Install the Web Interface:**
    ```bash
    pip install streamlit
    ```
2.  **Install Scientific Libraries & Utilities:**
    ```bash
    pip install rdkit meeko py3Dmol stmol pandas numpy scipy ipywidgets ipython_genutils
    ```

> **Troubleshooting**:
> *   If you see "No module named 'ipython_genutils'", simply run:
>     ```bash
>     pip install ipython_genutils
>     ```
> *   If `meeko` fails or you need `gemmi`, ensure your virtual environment is active and install it manually:
>     ```bash
>     source .venv/bin/activate
>     pip install gemmi
>     ```
> *   If you see a `SyntaxError` about `compute_water_map.py` during `meeko` install, **ignore it**. It is a harmless warning.

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

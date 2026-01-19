# Computational Analysis of Tau Protein Aggregation

## Abstract
This project performs an *in silico* screening of ligands targeting the Tau protein fibril, specifically the VQIVYK (PHF6) aggregation motif. Using AutoDock Vina for molecular docking and an evolutionary algorithm for generative design, we identify potential inhibitors that could prevent the formation of neurofibrillary tangles associated with Alzheimer's disease. The study compares FDA-approved drugs, flavonoids, and known inhibitors, and further optimizes lead compounds to maximize binding affinity across multiple protein conformations.

## Installation

### Prerequisites
*   Python 3.8 or higher
*   Git (optional, for cloning)

### Quick Start (Recommended)

It is highly recommended to run this project in a virtual environment to avoid dependency conflicts.

**1. Clone or Download the Repository**
```bash
git clone <repository_url>
cd <repository_folder>
```

**2. Create and Activate a Virtual Environment**

*   **Windows:**
    ```bash
    python -m venv venv
    .\venv\Scripts\activate
    ```
*   **macOS / Linux:**
    ```bash
    python3 -m venv venv
    source venv/bin/activate
    ```

**3. Install Dependencies**
Once the virtual environment is active, install the required packages:
```bash
pip install -r requirements.txt
```

> **Note:** This installation includes `rdkit`, `meeko`, `scipy`, `streamlit`, and other necessary scientific libraries. `gemmi` and `ipython_genutils` are also included to ensure compatibility with structure preparation and the interface.

## Usage

To launch the interactive dashboard, run one of the following commands in your terminal:

```bash
streamlit run app.py
```

**Alternative Command:**
If the command above does not work (e.g., "command not found"), try running via the Python module:

```bash
python3 -m streamlit run app.py
```

> **Important:** Do **not** run the application using `python app.py` or `python3 app.py`. This will not start the web server correctly.

This will open the application in your default web browser (usually at `http://localhost:8501`).

### Dashboard Features:
*   **Research Findings:** Explore the data comparing single-structure vs. ensemble docking strategies, including binding energy distributions and ligand stability metrics.
*   **Interactive Validation:**
    *   **3D Visualization:** Interactively view 3D docked poses of ligands and the receptor with customizable styles.
    *   **AI Optimization:** Real-time demonstration of the Evolutionary Algorithm optimizing a molecule structure.

## Reproducing Science Fair Results

To generate the raw data used for the "Shapeshifter" analysis—which compares standard single-structure docking against an ensemble-based evolutionary approach—run the following script:

```bash
python experiment_phase3_shapeshifter.py
```

### What this script does:
1.  **Phase 1 (Data Prep):** Downloads the NMR ensemble (PDB: 2MZ7), splits it into individual models, and calculates geometric centers for the binding pocket.
2.  **Phase 2 (Optimization):**
    *   **Experiment A:** Evolves a ligand to bind to a *single* protein conformation (Model 1).
    *   **Experiment B:** Evolves a ligand to bind to an *ensemble* of conformations (Models 1, 5, 10), rewarding robustness.
3.  **Phase 3 (Validation):** Cross-docks both evolved molecules against the ensemble to determine which strategy yields a more reliable inhibitor.

## License
This project is for educational and research purposes.

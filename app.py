import streamlit as st
import pandas as pd
import os
import py3Dmol
from stmol import showmol
import plotly.express as px
import docking
import evolution
import scoring
from typing import Tuple

st.set_page_config(page_title="Tau Protein Docking Sim", layout="wide")

# Load Config for Visualization (Legacy single receptor view)
@st.cache_resource
def load_config() -> Tuple[float, float, float]:
    if not os.path.exists("data/config.txt"):
        return (0.0, 0.0, 0.0)
    with open("data/config.txt") as f:
        c = f.read().strip().split(',')
        return (float(c[0]), float(c[1]), float(c[2]))

CENTER = load_config()
RECEPTOR_PATH = "data/receptor.pdbqt" # Use this for visualizer default

# Title
st.title("In Silico Screening of Ligands for Tau Protein Aggregation")

# Shapeshifter Text
SHAPESHIFTER_TEXT = """
**The "Shapeshifter" Hypothesis:**

The Tau protein is an **Intrinsically Disordered Protein (IDP)**, meaning it does not have a single fixed structure but fluctuates between many conformations. Traditional drug discovery targets a single rigid structure, often leading to 'brittle' drug candidates that fail when the protein shifts shape.

This project utilizes an **Ensemble Docking** approach. By optimizing ligands against a diverse set of Tau conformations simultaneously (represented by NMR models from PDB: 2MZ7), I hypothesize that we can identify **Robust Capping Inhibitors**. These molecules may have a slightly lower peak affinity for any single shape, but they maintain a consistent, stable binding interaction across the protein's entire structural landscape, effectively 'capping' the aggregation-prone PHF6 motif regardless of thermal fluctuations.
"""

tabs = st.tabs(["Research Findings", "Interactive Validation"])

# --- TAB 1: RESEARCH FINDINGS ---
with tabs[0]:
    st.header("Research Findings: Single vs. Ensemble Optimization")

    st.markdown("### Abstract")
    st.markdown("""
    This project performs an *in silico* screening of ligands targeting the Tau protein fibril, specifically the VQIVYK (PHF6) aggregation motif. Using AutoDock Vina for molecular docking and an evolutionary algorithm for generative design, we identify potential inhibitors that could prevent the formation of neurofibrillary tangles associated with Alzheimer's disease. The study compares FDA-approved drugs, flavonoids, and known inhibitors, and further optimizes lead compounds to maximize binding affinity across multiple protein conformations.
    """)

    st.markdown("---")
    st.markdown(SHAPESHIFTER_TEXT)
    st.markdown("---")

    # Load Results
    results_path = "results_ensemble_data.csv"
    if os.path.exists(results_path):
        df_results = pd.read_csv(results_path)

        st.subheader("Data Analysis")

        col1, col2 = st.columns(2)

        with col1:
            st.markdown("#### Binding Energy Distribution")
            st.caption("Lower (more negative) is better. Compares the stability of the drug across different protein shapes.")
            fig_box = px.box(
                df_results,
                x="Experiment_Mode",
                y="Binding_Energy",
                color="Experiment_Mode",
                points="all",
                title="Binding Energy Consistency (Single vs. Ensemble)",
                labels={"Binding_Energy": "Binding Energy (kcal/mol)", "Experiment_Mode": "Strategy"}
            )
            st.plotly_chart(fig_box, use_container_width=True)

        with col2:
            st.markdown("#### Ligand Stability (Distance to Zipper)")
            st.caption("Distance from the ligand center to the target VQIVYK motif center. Consistent low distance implies robust targeting.")
            fig_scatter = px.scatter(
                df_results,
                x="Model_ID",
                y="Distance_to_Zipper",
                color="Experiment_Mode",
                symbol="Experiment_Mode",
                size=[10]*len(df_results),
                title="Ligand Deviation across Protein Conformations",
                labels={"Distance_to_Zipper": "Distance to VQIVYK Center (Å)", "Model_ID": "Protein Model (Conformation)"}
            )
            # Add line to show trend if desired, or just scatter
            fig_scatter.update_traces(mode='markers+lines')
            st.plotly_chart(fig_scatter, use_container_width=True)

        st.dataframe(df_results)

    else:
        st.warning("No experimental data found. Please run 'experiment_phase3_shapeshifter.py' to generate results.")

# --- TAB 2: INTERACTIVE VALIDATION ---
with tabs[1]:
    st.header("Interactive Validation Tools")
    st.markdown("Experiment with the models and the evolutionary algorithm in real-time.")

    col_tools = st.columns([1, 1])

    with col_tools[0]:
        st.subheader("3D Visualization")
        st.markdown("Visualize the molecular interactions. Select a molecule to see its docked pose.")

        # Simple Molecule Selector for Visualization
        # We can let them type a SMILES or pick a default
        smi_input = st.text_input("Enter SMILES for Visualization", "c1ccccc1")
        mol_name = st.text_input("Molecule Name", "Benzene")

        receptor_style = st.selectbox("Receptor Style", ["Cartoon", "Stick", "Line", "Sphere"], index=0)
        ligand_style = st.selectbox("Ligand Style", ["Stick", "Line", "Sphere"], index=0)

        if st.button("Dock & Visualize"):
            with st.spinner("Calculating pose..."):
                pdbqt = docking.prep_ligand(smi_input, mol_name)
                score = None
                docked_pose = None
                if pdbqt:
                    score, docked_pose = docking.run_docking(pdbqt, RECEPTOR_PATH, CENTER)
                else:
                    st.error(f"Failed to prepare ligand {mol_name}.")

            if docked_pose:
                st.metric("Binding Score", f"{score} kcal/mol")

                view = py3Dmol.view(width=600, height=500)

                # Load Receptor
                with open("data/5O3L.pdb") as f: # Use the nice PDB for visuals
                    pdb_content = f.read()
                view.addModel(pdb_content, "pdb")

                 # Apply Receptor Style
                rec_style_map = {
                    "Cartoon": {"cartoon": {'color': 'spectrum'}},
                    "Stick": {"stick": {}},
                    "Line": {"line": {}},
                    "Sphere": {"sphere": {}}
                }
                view.setStyle({'model': -1}, rec_style_map.get(receptor_style, {"cartoon": {'color': 'spectrum'}}))

                # Add Ligand
                view.addModel(docked_pose, "pdbqt")

                # Apply Ligand Style
                lig_style_map = {
                    "Stick": {"stick": {'colorscheme': 'greenCarbon'}},
                    "Line": {"line": {'colorscheme': 'greenCarbon'}},
                    "Sphere": {"sphere": {'colorscheme': 'greenCarbon'}}
                }
                view.setStyle({'model': -1}, lig_style_map.get(ligand_style, {"stick": {'colorscheme': 'greenCarbon'}}))

                view.zoomTo()
                showmol(view, height=500, width=600)

    with col_tools[1]:
        st.subheader("AI Evolutionary Optimization")
        st.markdown("Run the Genetic Algorithm live to optimize a starting molecule.")

        start_smi = st.text_input("Starting SMILES", "c1ccccc1", key="evo_smi")
        generations = st.slider("Generations", 1, 10, 3)

        if st.button("Start Evolution"):
            st.write(f"Evolving from: {start_smi}")

            with st.spinner("Evolving..."):
                # Use single receptor for speed in live demo
                best_score, best_smi, best_pose, history = evolution.evolve_ligand(
                    start_smi, RECEPTOR_PATH, CENTER, generations=generations, population_size=4
                )

            st.success("Evolution Complete!")

            col_res1, col_res2 = st.columns(2)
            with col_res1:
                st.metric("Best Score", f"{best_score} kcal/mol")
            with col_res2:
                st.code(best_smi, language="text")

            hist_df = pd.DataFrame(history)
            fig_opt = px.line(
                hist_df,
                x="generation",
                y="best_score",
                title="Optimization Trajectory",
                markers=True
            )
            st.plotly_chart(fig_opt, use_container_width=True)

st.markdown("---")
st.caption("Developed for Science Fair Project | Sources: PDB, PubChem, PubMed.")

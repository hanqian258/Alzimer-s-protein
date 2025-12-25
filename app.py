import streamlit as st
import pandas as pd
import os
import py3Dmol
from stmol import showmol
import docking
import evolution
import time

st.set_page_config(page_title="Tau Protein Docking Sim", layout="wide")

# Load Config
@st.cache_resource
def load_config():
    if not os.path.exists("data/config.txt"):
        return (0.0, 0.0, 0.0)
    with open("data/config.txt") as f:
        c = f.read().strip().split(',')
        return (float(c[0]), float(c[1]), float(c[2]))

CENTER = load_config()
RECEPTOR_PATH = "data/receptor.pdbqt"

st.title("In Silico Screening of Ligands for Tau Protein Aggregation")
st.markdown("""
**Research Aim:** This simulation compares the binding affinity of various ligands to the **Tau Protein Fibril** (PDB: 5O3L), specifically targeting the **VQIVYK (PHF6)** aggregation motif.
We analyze three groups of molecules: **FDA Approved Drugs** (Repurposed), **Flavonoids**, and **Known Ligands** (PET Tracers/Inhibitors).
Additionally, we employ an **Evolutionary Algorithm** to virtually evolve and optimize ligand structures for better binding.
""")

tabs = st.tabs(["Dashboard & Molecules", "Docking Simulation", "Visualization", "AI Optimization"])

# Load Data
@st.cache_data
def load_data():
    return pd.read_csv("data/molecules.csv")

df = load_data()

with tabs[0]:
    st.header("Molecule Library")
    st.markdown("The following molecules have been selected based on literature reviews for their potential to bind Tau aggregates.")

    col1, col2, col3 = st.columns(3)
    with col1:
        st.subheader("FDA/Clinical Drugs")
        st.dataframe(df[df['category'] == 'FDA Drug'][['name', 'citation']])
    with col2:
        st.subheader("Flavonoids")
        st.dataframe(df[df['category'] == 'Flavonoid'][['name', 'citation']])
    with col3:
        st.subheader("Known Ligands")
        st.dataframe(df[df['category'] == 'Known Ligand'][['name', 'citation']])

    st.info("Navigate to the **Docking Simulation** tab to run the binding affinity calculations.")

with tabs[1]:
    st.header("High-Throughput Docking Simulation")
    st.markdown(f"**Target:** Tau Fibril (5O3L) | **Grid Center:** {CENTER}")

    if st.button("Run Docking Simulation"):
        st.write("Initializing AutoDock Vina...")
        progress_bar = st.progress(0)
        results = []

        total = len(df)
        for i, row in df.iterrows():
            st.text(f"Docking {row['name']} ({row['category']})...")

            # Prepare Ligand
            pdbqt = docking.prep_ligand(row['smiles'], row['name'])

            score = None
            if pdbqt:
                score, _ = docking.run_docking(pdbqt, RECEPTOR_PATH, CENTER)

            results.append({
                "Name": row['name'],
                "Category": row['category'],
                "Binding Score (kcal/mol)": score if score else 0
            })

            progress_bar.progress((i + 1) / total)

        results_df = pd.DataFrame(results)
        st.success("Simulation Complete!")

        # Save results for visualization
        results_df.to_csv("data/docking_results.csv", index=False)
        st.session_state['docking_results'] = results_df

    if 'docking_results' in st.session_state:
        res = st.session_state['docking_results']

        st.subheader("Results Analysis")
        st.dataframe(res.sort_values("Binding Score (kcal/mol)"))

        st.subheader("Comparative Binding Affinity")
        st.bar_chart(res, x="Name", y="Binding Score (kcal/mol)", color="Category")

        st.write("**Note:** Lower (more negative) scores indicate stronger binding to the Tau fibril.")

with tabs[2]:
    st.header("3D Visualization")
    st.markdown("Visualize the molecular interactions. Select a molecule to see its docked pose with the Tau Fibril.")

    # Needs docking to be run first for specific poses, or we re-run on demand for visualization
    selected_mol_name = st.selectbox("Select Molecule", df['name'].unique())

    if st.button("Visualize Interaction"):
        row = df[df['name'] == selected_mol_name].iloc[0]

        with st.spinner("Calculating pose..."):
            pdbqt = docking.prep_ligand(row['smiles'], row['name'])
            score, docked_pose = docking.run_docking(pdbqt, RECEPTOR_PATH, CENTER)

        if docked_pose:
            st.metric("Binding Score", f"{score} kcal/mol")

            # View
            view = py3Dmol.view(width=800, height=600)

            # Load Receptor (cleaned PDB or PDBQT)
            # We display the original PDB for better visuals (cartoon)
            with open("data/5O3L.pdb") as f:
                pdb_content = f.read()
            view.addModel(pdb_content, "pdb")
            view.setStyle({'model': -1}, {"cartoon": {'color': 'spectrum'}})

            # Add Ligand
            view.addModel(docked_pose, "pdbqt")
            view.setStyle({'model': -1}, {"stick": {'colorscheme': 'greenCarbon'}})

            view.zoomTo()
            showmol(view, height=600, width=800)
        else:
            st.error("Docking failed for visualization.")

with tabs[3]:
    st.header("AI-Driven Optimization (Evolutionary Algorithm)")
    st.markdown("""
    This module uses a Genetic Algorithm to "evolve" a molecule structure to maximize its binding affinity to the **Tau VQIVYK motif**.
    It starts with a parent molecule, applies random mutations (Back Propagation concept), and selects the best binder.
    """)

    parent_name = st.selectbox("Select Parent Molecule", df['name'].unique(), index=0)
    generations = st.slider("Generations", 1, 10, 5)

    if st.button("Start Evolution"):
        row = df[df['name'] == parent_name].iloc[0]
        initial_smiles = row['smiles']

        st.write(f"Starting evolution from: **{parent_name}**")
        st.code(initial_smiles)

        placeholder = st.empty()
        chart_placeholder = st.empty()

        with st.spinner("Evolving..."):
             best_score, best_smi, best_pose, history = evolution.evolve_ligand(
                 initial_smiles, RECEPTOR_PATH, CENTER, generations=generations, population_size=4
             )

        st.success("Evolution Complete!")

        col1, col2 = st.columns(2)
        with col1:
            st.metric("Original Score", f"{history[0]['best_score'] if history else 0} kcal/mol") # Approx
        with col2:
            st.metric("Optimized Score", f"{best_score} kcal/mol")

        st.subheader("Optimization Trajectory")
        hist_df = pd.DataFrame(history)
        st.line_chart(hist_df.set_index("generation")['best_score'])

        st.subheader("Optimized Structure")
        st.code(best_smi)

        # Visualize Optimized
        view = py3Dmol.view(width=800, height=600)
        with open("data/5O3L.pdb") as f:
            pdb_content = f.read()
        view.addModel(pdb_content, "pdb")
        view.setStyle({'model': -1}, {"cartoon": {'color': 'spectrum'}})
        view.addModel(best_pose, "pdbqt")
        view.setStyle({'model': -1}, {"stick": {'colorscheme': 'cyanCarbon'}})
        view.zoomTo()
        showmol(view, height=600, width=800)

st.markdown("---")
st.caption("Developed for Science Fair Project | Sources: PDB, PubChem, PubMed.")

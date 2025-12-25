import streamlit as st
import pandas as pd
import os
import py3Dmol
from stmol import showmol
import plotly.express as px
import docking
import evolution
import scoring
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
        df_fda = df[df['category'] == 'FDA Drug'][['name', 'citation']]
        st.dataframe(df_fda)
        st.download_button(
            label="Download CSV",
            data=df_fda.to_csv(index=False).encode('utf-8'),
            file_name='fda_drugs.csv',
            mime='text/csv',
            key='download-fda'
        )
    with col2:
        st.subheader("Flavonoids")
        df_flav = df[df['category'] == 'Flavonoid'][['name', 'citation']]
        st.dataframe(df_flav)
        st.download_button(
            label="Download CSV",
            data=df_flav.to_csv(index=False).encode('utf-8'),
            file_name='flavonoids.csv',
            mime='text/csv',
            key='download-flav'
        )
    with col3:
        st.subheader("Known Ligands")
        df_known = df[df['category'] == 'Known Ligand'][['name', 'citation']]
        st.dataframe(df_known)
        st.download_button(
            label="Download CSV",
            data=df_known.to_csv(index=False).encode('utf-8'),
            file_name='known_ligands.csv',
            mime='text/csv',
            key='download-known'
        )

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

            # Calculate Metrics
            kd_val = scoring.calculate_kd(score) if score else 0
            affinity_str = scoring.interpret_affinity(score) if score else "N/A"

            results.append({
                "Name": row['name'],
                "Category": row['category'],
                "Binding Score (kcal/mol)": score if score else 0,
                "Kd (uM)": f"{kd_val:.4f}" if score else "N/A",
                "Affinity Level": affinity_str
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
        st.download_button(
            label="Download Results CSV",
            data=res.to_csv(index=False).encode('utf-8'),
            file_name='docking_results.csv',
            mime='text/csv',
            key='download-results'
        )

        st.subheader("Comparative Binding Affinity")
        fig = px.bar(
            res,
            x="Name",
            y="Binding Score (kcal/mol)",
            color="Category",
            title="Binding Affinity by Molecule"
        )
        st.plotly_chart(fig, use_container_width=True)

        st.write("**Note:** Lower (more negative) scores indicate stronger binding to the Tau fibril.")

with tabs[2]:
    st.header("3D Visualization")
    st.markdown("Visualize the molecular interactions. Select a molecule to see its docked pose with the Tau Fibril.")

    # Needs docking to be run first for specific poses, or we re-run on demand for visualization
    col_sel1, col_sel2, col_sel3 = st.columns(3)
    with col_sel1:
        selected_mol_name = st.selectbox("Select Molecule", df['name'].unique())
    with col_sel2:
        receptor_style = st.selectbox("Receptor Style", ["Cartoon", "Stick", "Line", "Sphere"], index=0)
    with col_sel3:
        ligand_style = st.selectbox("Ligand Style", ["Stick", "Line", "Sphere"], index=0)

    if st.button("Visualize Interaction"):
        row = df[df['name'] == selected_mol_name].iloc[0]

        with st.spinner("Calculating pose..."):
            pdbqt = docking.prep_ligand(row['smiles'], row['name'])
            score = None
            docked_pose = None
            if pdbqt:
                score, docked_pose = docking.run_docking(pdbqt, RECEPTOR_PATH, CENTER)
            else:
                st.error(f"Failed to prepare ligand {row['name']}.")

        if docked_pose:
            col_res1, col_res2, col_res3 = st.columns(3)
            with col_res1:
                st.metric("Binding Score", f"{score} kcal/mol")
            with col_res2:
                kd_val = scoring.calculate_kd(score)
                st.metric("Dissociation Constant (Kd)", f"{kd_val:.2f} uM")
            with col_res3:
                interp = scoring.interpret_affinity(score)
                color = scoring.get_score_color(score)
                st.markdown(f"### Affinity: <span style='color:{color}'>{interp}</span>", unsafe_allow_html=True)

            # View
            view = py3Dmol.view(width=800, height=600)

            # Load Receptor (cleaned PDB or PDBQT)
            # We display the original PDB for better visuals (cartoon)
            with open("data/5O3L.pdb") as f:
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
            showmol(view, height=600, width=800)

            st.markdown("### Download Structure Data")
            col_dl1, col_dl2 = st.columns(2)
            with col_dl1:
                st.download_button(
                    label="Download Docked Ligand (PDBQT)",
                    data=docked_pose,
                    file_name=f"{selected_mol_name}_docked.pdbqt",
                    mime="text/plain",
                    key='dl-docked'
                )
            with col_dl2:
                st.download_button(
                    label="Download Receptor (PDB)",
                    data=pdb_content,
                    file_name="5O3L_receptor.pdb",
                    mime="text/plain",
                    key='dl-receptor'
                )

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

        initial_score = history[0]['best_score'] if history else 0
        initial_pose = history[0]['best_pose'] if history and 'best_pose' in history[0] else None

        improvement_pct = 0.0
        if initial_score != 0:
            improvement_pct = ((best_score - initial_score) / initial_score) * 100

        col1, col2, col3, col4 = st.columns(4)
        with col1:
            st.metric("Original Score", f"{initial_score} kcal/mol")
        with col2:
            st.metric("Optimized Score", f"{best_score} kcal/mol")
        with col3:
             st.metric("Improvement", f"{improvement_pct:.2f}%")
        with col4:
             kd_val = scoring.calculate_kd(best_score)
             st.metric("Optimized Kd", f"{kd_val:.4f} uM")

        st.subheader("Optimization Trajectory")
        hist_df = pd.DataFrame(history)
        fig_opt = px.line(
            hist_df,
            x="generation",
            y="best_score",
            title="Binding Score Optimization over Generations",
            markers=True
        )
        st.plotly_chart(fig_opt, use_container_width=True)

        st.subheader("Optimized Structure")
        st.code(best_smi)

        if best_pose:
            st.download_button(
                label="Download Optimized Ligand (PDBQT)",
                data=best_pose,
                file_name=f"optimized_{parent_name}.pdbqt",
                mime="text/plain",
                key='dl-optimized'
            )

        # Visual Comparison
        st.subheader("Visual Comparison")
        col_v1, col_v2 = st.columns(2)

        with open("data/5O3L.pdb") as f:
            pdb_content = f.read()

        with col_v1:
            st.markdown("### Original Molecule")
            if initial_pose:
                view1 = py3Dmol.view(width=400, height=400)
                view1.addModel(pdb_content, "pdb")
                view1.setStyle({'model': -1}, {"cartoon": {'color': 'spectrum'}})
                view1.addModel(initial_pose, "pdbqt")
                view1.setStyle({'model': -1}, {"stick": {'colorscheme': 'grayCarbon'}})
                view1.zoomTo()
                showmol(view1, height=400, width=400)
            else:
                st.write("No pose available for original molecule.")

        with col_v2:
            st.markdown("### Optimized Molecule")
            if best_pose:
                view2 = py3Dmol.view(width=400, height=400)
                view2.addModel(pdb_content, "pdb")
                view2.setStyle({'model': -1}, {"cartoon": {'color': 'spectrum'}})
                view2.addModel(best_pose, "pdbqt")
                view2.setStyle({'model': -1}, {"stick": {'colorscheme': 'cyanCarbon'}})
                view2.zoomTo()
                showmol(view2, height=400, width=400)
            else:
                st.write("No pose available for optimized molecule.")

st.markdown("---")
st.caption("Developed for Science Fair Project | Sources: PDB, PubChem, PubMed.")

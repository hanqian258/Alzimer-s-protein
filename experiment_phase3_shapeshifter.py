import json
import os
import prep_shapes
import evolution
import docking

def main():
    print("=== Phase 1: Data Preparation ===")
    prep_shapes.main()

    # Load centers
    centers_path = "data/2MZ7_centers.json"
    if not os.path.exists(centers_path):
        print("Error: Centers file not found.")
        return

    with open(centers_path, 'r') as f:
        centers = json.load(f)

    print(f"Loaded centers for {len(centers)} models.")

    # Define Configurations

    # Exp A: Model 1 only
    config_a = [{
        'path': "data/2MZ7_model_1.pdbqt",
        'center': centers.get("1")
    }]

    # Exp B: Models 1, 5, 10
    config_b = [
        {
            'path': "data/2MZ7_model_1.pdbqt",
            'center': centers.get("1")
        },
        {
            'path': "data/2MZ7_model_5.pdbqt",
            'center': centers.get("5")
        },
        {
            'path': "data/2MZ7_model_10.pdbqt",
            'center': centers.get("10")
        }
    ]

    # Validation Config (Models 1, 5, 10) - same as B for scoring
    config_val = config_b

    initial_smiles = "c1ccccc1" # Benzene start

    print("\n=== Phase 2: Experiment A (Single Structure) ===")
    print("Running GA on Model 1...")
    score_a, smi_a, pose_a, hist_a = evolution.evolve_ligand(
        initial_smiles,
        config_a,
        generations=2,
        population_size=4
    )
    print(f"Best Molecule A: {smi_a} (Score on Model 1: {score_a})")

    print("\n=== Phase 2: Experiment B (Shapeshifter Ensemble) ===")
    print("Running GA on Models 1, 5, 10...")
    score_b, smi_b, pose_b, hist_b = evolution.evolve_ligand(
        initial_smiles,
        config_b,
        generations=2,
        population_size=4
    )
    print(f"Best Molecule B: {smi_b} (Avg Score on Ensemble: {score_b})")

    print("\n=== Phase 3: The Gauntlet (Validation) ===")

    # Dock A against Ensemble
    pdbqt_a = docking.prep_ligand(smi_a, "Mol_A")
    val_score_a, _ = docking.run_docking_on_list(pdbqt_a, config_val)

    # Dock B against Ensemble
    pdbqt_b = docking.prep_ligand(smi_b, "Mol_B")
    val_score_b, _ = docking.run_docking_on_list(pdbqt_b, config_val)

    print("-" * 60)
    print(f"{'Metric':<30} | {'Drug A (Single)':<12} | {'Drug B (Ensemble)':<12}")
    print("-" * 60)
    print(f"{'SMILES':<30} | {smi_a:<12} | {smi_b:<12}")
    print(f"{'Training Score':<30} | {score_a:<12.2f} | {score_b:<12.2f}")
    print(f"{'Validation Score (Avg 1,5,10)':<30} | {val_score_a:<12.2f} | {val_score_b:<12.2f}")
    print("-" * 60)

    if val_score_b < val_score_a:
        print("\nConclusion: The Ensemble-trained drug (B) is more robust across conformations.")
    else:
        print("\nConclusion: The Single-structure drug (A) performed better or equal (Result might vary with short run).")

if __name__ == "__main__":
    main()

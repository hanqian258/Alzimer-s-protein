import json
import os
import pandas as pd
import prep_shapes
import evolution
import docking
import scoring

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

    # Validation Config (Models 1, 5, 10)
    # We will manually iterate over these for detailed logging
    validation_models = ["1", "5", "10"]

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

    results = []

    # Function to run validation
    def run_validation(smiles, experiment_mode, label):
        pdbqt = docking.prep_ligand(smiles, label)
        if not pdbqt:
            print(f"Failed to prep {label}")
            return

        for model_id in validation_models:
            model_path = f"data/2MZ7_model_{model_id}.pdbqt"
            center = centers.get(model_id)

            if not os.path.exists(model_path) or not center:
                print(f"Skipping Model {model_id}, missing file or center.")
                continue

            # Run Docking
            score, pose = docking.run_docking(pdbqt, model_path, center)

            dist = 0.0
            if pose:
                # Calculate Centroid of Ligand
                lig_centroid = scoring.get_ligand_centroid(pose)
                if lig_centroid:
                    # Calculate Distance to Receptor Zipper Center
                    dist = scoring.calculate_distance(lig_centroid, center)

            # Log Data
            results.append({
                "Experiment_Mode": experiment_mode,
                "Model_ID": f"Model_{model_id}",
                "Binding_Energy": score if score is not None else 0.0,
                "Distance_to_Zipper": dist,
                "Ligand_SMILES": smiles
            })
            print(f"  > {label} vs Model {model_id}: Score={score}, Dist={dist:.2f}")

    print("Validating Single-Structure Drug (A)...")
    run_validation(smi_a, "Single (A)", "Mol_A")

    print("Validating Ensemble Drug (B)...")
    run_validation(smi_b, "Ensemble (B)", "Mol_B")

    # Save to CSV
    csv_path = "results_ensemble_data.csv"
    df = pd.DataFrame(results)
    df.to_csv(csv_path, index=False)
    print(f"\nResults saved to {csv_path}")
    print(df)

if __name__ == "__main__":
    main()

import random
from typing import List, Tuple, Any, Optional, Dict
from rdkit import Chem
import docking


def mutate_molecule(smiles: str, mutation_rate: float = 0.1) -> str:
    """
    Applies random mutations to a molecule string.

    This function attempts to mimic evolutionary mutations by modifying the molecular structure.
    Supported mutations:
    - Add a carbon atom (methyl group) to a random atom.
    - Remove a random terminal atom.
    - Change an atom type (substitute C, N, O, F).

    Args:
        smiles (str): The SMILES string of the molecule to mutate.
        mutation_rate (float): The probability of mutation (currently unused but reserved for future logic).

    Returns:
        str: The SMILES string of the mutated molecule, or the original SMILES if mutation fails.
    """
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return smiles

        # Clone molecule
        rw_mol = Chem.RWMol(mol)
        atoms = [a for a in rw_mol.GetAtoms()]
        num_atoms = len(atoms)

        if num_atoms == 0:
            return smiles

        mutation_type = random.choice(['add', 'remove', 'substitute'])

        if mutation_type == 'add':
            # Add a methyl group to a random atom that has valency
            idx = random.randint(0, num_atoms - 1)
            atom = rw_mol.GetAtomWithIdx(idx)
            if atom.GetImplicitValence() > 0:
                new_idx = rw_mol.AddAtom(Chem.Atom(6))  # Add Carbon
                rw_mol.AddBond(idx, new_idx, Chem.BondType.SINGLE)

        elif mutation_type == 'remove':
            # Remove a terminal atom (degree 1)
            terminal_atoms = [a.GetIdx() for a in atoms if a.GetDegree() == 1]
            if terminal_atoms:
                idx_to_remove = random.choice(terminal_atoms)
                rw_mol.RemoveAtom(idx_to_remove)

        elif mutation_type == 'substitute':
            # Change random atom to C, N, O, F
            idx = random.randint(0, num_atoms - 1)
            atom = rw_mol.GetAtomWithIdx(idx)
            current_atomic_num = atom.GetAtomicNum()
            choices = [6, 7, 8, 9]  # C, N, O, F
            if current_atomic_num in choices:
                choices.remove(current_atomic_num)
            new_atomic_num = random.choice(choices)
            atom.SetAtomicNum(new_atomic_num)

        # Sanitize to check validity
        try:
            Chem.SanitizeMol(rw_mol)
            new_smiles = Chem.MolToSmiles(rw_mol)
            return new_smiles
        except Exception:
            return smiles  # Return original if mutation failed

    except Exception:
        # print(f"Mutation failed: {e}")
        return smiles

def evolve_ligand(initial_smiles, receptor_configs, center=None, generations=5, population_size=5):
    """
    Genetic Algorithm to evolve the ligand.

    receptor_configs: Can be:
      - A single path string (legacy support).
      - A list of dicts [{'path': str, 'center': tuple}, ...] for multi-target docking.

    center: Only used if receptor_configs is a string path.
    """

    # Normalize receptor_configs
    if isinstance(receptor_configs, str):
        if center is None:
            raise ValueError("Center must be provided if receptor_configs is a path string.")
        configs = [{'path': receptor_configs, 'center': center}]
    else:
        configs = receptor_configs

    current_population = [initial_smiles]
    best_overall_score = 0.0
    best_overall_smiles = initial_smiles
    best_overall_pdbqt = None

    history = []

    # Initialize Score
    print("Evaluating initial structure...")
    pdbqt = docking.prep_ligand(initial_smiles, "Initial")
    if pdbqt:
        score, docked_pose = docking.run_docking_on_list(pdbqt, configs)
        if score is not None:
            best_overall_score = score
            best_overall_pdbqt = docked_pose
            print(f"Initial Score: {score}")
            history.append({
                "generation": 0,
                "best_score": best_overall_score,
                "best_smiles": initial_smiles,
                "best_pose": best_overall_pdbqt
            })
        else:
            best_overall_score = 0.0  # Failed

    for gen in range(generations):
        print(f"--- Generation {gen+1} ---")
        # Mutate to create offspring
        offspring = []
        for _ in range(population_size):
            parent = random.choice(current_population)
            child_smiles = mutate_molecule(parent)
            offspring.append(child_smiles)

        # Evaluate offspring
        results = []
        for i, smi in enumerate(offspring):
            pdbqt = docking.prep_ligand(smi, f"Gen{gen}_Mol{i}")
            if pdbqt:
                score, docked_pose = docking.run_docking_on_list(pdbqt, configs)
                if score is not None:
                    results.append((score, smi, docked_pose))

        # Sort by score (lower is better for affinity)
        results.sort(key=lambda x: x[0])

        if results:
            best_gen_score = results[0][0]
            print(f"Best Score in Gen {gen+1}: {best_gen_score}")

            # Update overall best
            if best_gen_score < best_overall_score:
                best_overall_score = best_gen_score
                best_overall_smiles = results[0][1]
                best_overall_pdbqt = results[0][2]

            # Selection: Keep top 50% for next generation
            # Ensure diversity?
            current_population = [r[1] for r in results[:max(1, population_size // 2)]]

            history.append({
                "generation": gen + 1,
                "best_score": best_gen_score,
                "best_smiles": results[0][1],
                "best_pose": results[0][2]
            })
        else:
            print("No valid offspring in this generation.")

    return best_overall_score, best_overall_smiles, best_overall_pdbqt, history


if __name__ == "__main__":
    # Test Evolution
    try:
        with open("data/config.txt") as f:
            c = f.read().strip().split(',')
            center_coords = (float(c[0]), float(c[1]), float(c[2]))

        smi = "c1ccccc1"  # Benzene
        score, best_smi, _, hist = evolve_ligand(
            smi, "data/receptor.pdbqt", center_coords, generations=2, population_size=3
        )
        print(f"Final Best Score: {score}")
        print(f"Final Best SMILES: {best_smi}")
    except FileNotFoundError:
        print("Config file not found. Skipping main test.")

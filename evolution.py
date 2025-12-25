import random
from rdkit import Chem
from rdkit.Chem import AllChem
import docking
import pandas as pd

def mutate_molecule(smiles, mutation_rate=0.1):
    """
    Applies random mutations to a molecule string.
    Supported mutations:
    - Add a carbon atom (methyl group) to a random atom.
    - Remove a random terminal atom.
    - Change an atom type (C -> N, N -> O, etc).

    This is a simplified mutation strategy for demonstration.
    """
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None: return smiles

        # Clone molecule
        rw_mol = Chem.RWMol(mol)
        atoms = [a for a in rw_mol.GetAtoms()]
        num_atoms = len(atoms)

        if num_atoms == 0: return smiles

        mutation_type = random.choice(['add', 'remove', 'substitute'])

        if mutation_type == 'add':
            # Add a methyl group to a random atom that has valency
            idx = random.randint(0, num_atoms - 1)
            atom = rw_mol.GetAtomWithIdx(idx)
            if atom.GetImplicitValence() > 0:
                new_idx = rw_mol.AddAtom(Chem.Atom(6)) # Add Carbon
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
            choices = [6, 7, 8, 9] # C, N, O, F
            if current_atomic_num in choices:
                choices.remove(current_atomic_num)
            new_atomic_num = random.choice(choices)
            atom.SetAtomicNum(new_atomic_num)

        # Sanitize to check validity
        try:
            Chem.SanitizeMol(rw_mol)
            new_smiles = Chem.MolToSmiles(rw_mol)
            return new_smiles
        except:
            return smiles # Return original if mutation failed

    except Exception as e:
        # print(f"Mutation failed: {e}")
        return smiles

def evolve_ligand(initial_smiles, receptor_pdbqt, center, generations=5, population_size=5):
    """
    Genetic Algorithm to evolve the ligand.
    """
    current_population = [initial_smiles]
    best_overall_score = 0
    best_overall_smiles = initial_smiles
    best_overall_pdbqt = None

    # Initialize Score
    print("Evaluating initial structure...")
    pdbqt = docking.prep_ligand(initial_smiles, "Initial")
    if pdbqt:
        score, docked_pose = docking.run_docking(pdbqt, receptor_pdbqt, center)
        if score:
            best_overall_score = score
            best_overall_pdbqt = docked_pose
            print(f"Initial Score: {score}")
        else:
            best_overall_score = 0 # Failed

    history = []

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
                score, docked_pose = docking.run_docking(pdbqt, receptor_pdbqt, center)
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
            current_population = [r[1] for r in results[:max(1, population_size//2)]]

            history.append({
                "generation": gen + 1,
                "best_score": best_gen_score,
                "best_smiles": results[0][1]
            })
        else:
            print("No valid offspring in this generation.")

    return best_overall_score, best_overall_smiles, best_overall_pdbqt, history

if __name__ == "__main__":
    # Test Evolution
    with open("data/config.txt") as f:
        c = f.read().strip().split(',')
        center = (float(c[0]), float(c[1]), float(c[2]))

    smi = "c1ccccc1" # Benzene
    score, best_smi, _, hist = evolve_ligand(smi, "data/receptor.pdbqt", center, generations=2, population_size=3)
    print(f"Final Best Score: {score}")
    print(f"Final Best SMILES: {best_smi}")

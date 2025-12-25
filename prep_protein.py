import os
from rdkit import Chem
from rdkit.Chem import AllChem

def prep_protein(input_pdb, output_pdbqt):
    """
    Extracts the first model from the PDB file and converts it to PDBQT format.
    Calculates the center of the VQIVYK region (Residues 306-311 in full Tau numbering).

    For 5O3L (AD Tau Fibril), the chains are stacked.
    We will calculate the center of the VQIVYK motif (Val-Gln-Ile-Val-Tyr-Lys) across all chains
    to target the core of the filament, or pick a specific chain interface.

    In 5O3L PDB file, the residues are renumbered.
    We need to identify the VQIVYK sequence (Val-Gln-Ile-Val-Tyr-Lys).
    """

    # Load PDB
    mol = Chem.MolFromPDBFile(input_pdb, removeHs=False)
    if mol is None:
        raise ValueError(f"Could not load {input_pdb}")

    coords = []
    conf = mol.GetConformer()

    # VQIVYK corresponds to residues 306-311 in full Tau.
    # In PDB 5O3L, residue numbers might be preserved (306-378).
    # Let's check for residues 306 to 311.

    target_res_nums = set([306, 307, 308, 309, 310, 311])

    found_atoms = 0
    for atom in mol.GetAtoms():
        res_info = atom.GetPDBResidueInfo()
        if res_info:
            res_num = res_info.GetResidueNumber()
            if res_num in target_res_nums:
                pos = conf.GetAtomPosition(atom.GetIdx())
                coords.append((pos.x, pos.y, pos.z))
                found_atoms += 1

    if not coords:
        print("Warning: Residues 306-311 (VQIVYK) not found by ID. Attempting to target Center of Mass.")
        coords = [conf.GetAtomPosition(i) for i in range(mol.GetNumAtoms())]
    else:
        print(f"Found {found_atoms} atoms in the VQIVYK region.")

    xs = [c[0] for c in coords]
    ys = [c[1] for c in coords]
    zs = [c[2] for c in coords]

    center = (
        sum(xs) / len(xs),
        sum(ys) / len(ys),
        sum(zs) / len(zs)
    )

    print(f"Grid Center (X, Y, Z): {center}")

    # Writing PDBQT (Simplified heuristic)
    pdb_block = Chem.MolToPDBBlock(mol)
    pdbqt_lines = []

    atom_types = {
        'C': 'C', 'N': 'N', 'O': 'O', 'S': 'S', 'H': 'H', 'F': 'F',
        'Cl': 'Cl', 'Br': 'Br', 'I': 'I', 'P': 'P'
    }

    for line in pdb_block.split('\n'):
        if line.startswith("ATOM") or line.startswith("HETATM"):
            element = line[76:78].strip()
            if not element:
                 atom_name = line[12:16].strip()
                 element = atom_name[0]

            ad_type = atom_types.get(element, 'A')

            base = line[:60]
            base = f"{base:<70}"
            base += " 0.000"
            base += f" {ad_type:<2}"

            pdbqt_lines.append(base)
        elif line.startswith("TER") or line.startswith("END"):
            pdbqt_lines.append(line)

    with open(output_pdbqt, 'w') as f:
        f.write('\n'.join(pdbqt_lines))

    print(f"Written {output_pdbqt}")
    return center

if __name__ == "__main__":
    center = prep_protein("data/5O3L.pdb", "data/receptor.pdbqt")
    with open("data/config.txt", "w") as f:
        f.write(f"{center[0]},{center[1]},{center[2]}")

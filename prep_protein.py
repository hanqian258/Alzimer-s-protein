import os
from rdkit import Chem
from rdkit.Chem import AllChem

def prep_protein(input_pdb, output_pdbqt):
    """
    Extracts the first model from the PDB file and converts it to PDBQT format.
    Also calculates the center of the 'KLVFF' region (residues 16-20).
    Note: A proper PDBQT conversion usually requires adding charges and atom types.
    Here we use a simplified approach or RDKit if possible, but RDKit's PDBQT writer is limited.
    However, for the purpose of this project and without 'prepare_receptor4.py' (MGLTools),
    we will do a best-effort conversion.

    Actually, 1IYT is an NMR structure with Hydrogens already.
    Vina requires PDBQT.

    We will try to use Meeko if it supports receptor preparation, or OpenBabel if available.
    Since we only installed 'meeko', we check if it has PDBQT writing capabilities for receptors.
    If not, we will parse the PDB and write a simplified PDBQT (keeping coordinates)
    and assuming Gasteiger charges or similar are not strictly critical for this level of demo,
    OR we try to use a minimal python PDBQT writer.

    Update: Meeko is primarily for ligands.
    Let's try to load into RDKit and write as PDB, then converting to PDBQT might be tricky.

    Strategy:
    1. Read PDB with RDKit.
    2. Extract residues 16-20 to calculate center.
    3. Write the whole protein to PDBQT.
       Since writing a valid PDBQT is complex (autodock atom types), we will use a basic mapping.
    """

    # Load PDB
    mol = Chem.MolFromPDBFile(input_pdb, removeHs=False)
    if mol is None:
        raise ValueError(f"Could not load {input_pdb}")

    # 1IYT contains multiple models. RDKit MolFromPDBFile reads the first one by default?
    # Actually it reads the first model.

    # Calculate Center of KLVFF (Residues 16, 17, 18, 19, 20)
    # 1IYT sequence: DAEFRHDSGYEVHHQKLVFFAEDVGSNKGAIIGLMVGGVVIA
    # K=16, L=17, V=18, F=19, F=20

    coords = []
    conf = mol.GetConformer()

    # PDB residue numbering usually matches.
    for atom in mol.GetAtoms():
        res_info = atom.GetPDBResidueInfo()
        if res_info:
            res_num = res_info.GetResidueNumber()
            if 16 <= res_num <= 20:
                pos = conf.GetAtomPosition(atom.GetIdx())
                coords.append((pos.x, pos.y, pos.z))

    if not coords:
        print("Warning: Residues 16-20 not found. Using Center of Mass of protein.")
        coords = [conf.GetAtomPosition(i) for i in range(mol.GetNumAtoms())]

    xs = [c[0] for c in coords]
    ys = [c[1] for c in coords]
    zs = [c[2] for c in coords]

    center = (
        sum(xs) / len(xs),
        sum(ys) / len(ys),
        sum(zs) / len(zs)
    )

    print(f"Grid Center (X, Y, Z): {center}")

    # Writing PDBQT
    # Since we lack MGLTools, we will format the PDB data into PDBQT format.
    # The critical part for Vina is the ATOM lines with partial charges (optional) and Atom Types.
    # We will assume a simple mapping for Atom Types based on element.
    # This is a heuristic approach.

    pdb_block = Chem.MolToPDBBlock(mol)
    pdbqt_lines = []

    atom_types = {
        'C': 'C', 'N': 'N', 'O': 'O', 'S': 'S', 'H': 'H', 'F': 'F',
        'Cl': 'Cl', 'Br': 'Br', 'I': 'I', 'P': 'P'
    }

    for line in pdb_block.split('\n'):
        if line.startswith("ATOM") or line.startswith("HETATM"):
            # PDB format:
            # 0-6: Record name
            # 30-54: Coordinates
            # 76-78: Element

            # PDBQT extends PDB:
            # 70-76: Charge (often left as 0.00 for rigid receptor if not calculated)
            # 77-79: Atom Type (Autodock Type)

            element = line[76:78].strip()
            if not element:
                 # Try to infer from atom name (12-16)
                 atom_name = line[12:16].strip()
                 element = atom_name[0]

            ad_type = atom_types.get(element, 'A') # A for aromatic? Or just C/N/O

            # Simple AD4 typing:
            # If Carbon, check if aromatic (not easy from text line).
            # We'll stick to element types which Vina handles reasonably well for scoring.

            # Construct PDBQT line
            # Keep original line up to 66? PDB ends coords around 54. Temp factor to 66.
            # Vina needs cols 31-54 (coords).

            # We copy the line and append charge and type
            # Charge 0.00, Type

            # Python string formatting is strict for PDB.
            # line is likely 80 chars.
            base = line[:60] # Coords + Occupancy + Temp
            # Pad to right
            base = f"{base:<70}"
            # Charge
            base += " 0.000"
            # Type
            base += f" {ad_type:<2}"

            pdbqt_lines.append(base)
        elif line.startswith("TER") or line.startswith("END"):
            pdbqt_lines.append(line)

    with open(output_pdbqt, 'w') as f:
        f.write('\n'.join(pdbqt_lines))

    print(f"Written {output_pdbqt}")
    return center

if __name__ == "__main__":
    center = prep_protein("data/1IYT.pdb", "data/receptor.pdbqt")
    # Save center to a file for easy reading later?
    with open("data/config.txt", "w") as f:
        f.write(f"{center[0]},{center[1]},{center[2]}")

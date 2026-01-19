import os
import requests
import json
from rdkit import Chem
from rdkit.Chem import AllChem
from meeko import RDKitMoleculeSetup, PDBQTWriterLegacy

def write_pdbqt_manual(mol, output_path):
    """
    Fallback PDBQT writer that preserves Gasteiger charges.
    """
    lines = []
    atom_types = {
        'C': 'C', 'N': 'N', 'O': 'O', 'S': 'S', 'H': 'H', 'F': 'F',
        'Cl': 'Cl', 'Br': 'Br', 'I': 'I', 'P': 'P'
    }

    # Ensure charges are computed
    if not mol.HasProp('gasteiger_computed'):
        try:
            AllChem.ComputeGasteigerCharges(mol)
        except:
            pass

    conf = mol.GetConformer()

    for i, atom in enumerate(mol.GetAtoms()):
        pos = conf.GetAtomPosition(i)

        # Get Charge
        charge = 0.0
        if atom.HasProp('_GasteigerCharge'):
            try:
                charge = float(atom.GetProp('_GasteigerCharge'))
            except:
                pass

        # Get Element and Type
        symbol = atom.GetSymbol()
        ad_type = atom_types.get(symbol, 'A')

        # PDBQT Atom Line
        # ATOM      1  N   MET A   1      -3.327   6.162  16.037  1.00  0.00    -0.215 N
        # Formats vary, but generally:
        # 0-6: Record name "ATOM  "
        # 6-11: Atom serial
        # 12-16: Atom name
        # 17-20: Residue name
        # 21: Chain ID
        # 22-26: Res seq
        # 30-38: X
        # 38-46: Y
        # 46-54: Z
        # 54-60: Occ
        # 60-66: Temp
        # 70-76: Partial Charge
        # 77-79: Atom Type

        res_info = atom.GetPDBResidueInfo()
        atom_name = "ATOM"
        res_name = "UNK"
        chain = "A"
        res_seq = 1

        if res_info:
            atom_name = res_info.GetName().strip()
            res_name = res_info.GetResidueName().strip()
            chain = res_info.GetChainId().strip()
            res_seq = res_info.GetResidueNumber()
        else:
            atom_name = symbol

        line = f"ATOM  {i+1:>5} {atom_name:<4} {res_name:<3} {chain:1}{res_seq:>4}    {pos.x:>8.3f}{pos.y:>8.3f}{pos.z:>8.3f}  1.00  0.00    {charge:>6.3f} {ad_type:<2}"
        lines.append(line)

    with open(output_path, 'w') as f:
        f.write('\n'.join(lines))
    print(f"Manually wrote {output_path} with charges.")

def download_pdb(pdb_id, output_path):
    """
    Downloads PDB file from RCSB.
    """
    if os.path.exists(output_path):
        print(f"{output_path} already exists.")
        return

    url = f"https://files.rcsb.org/download/{pdb_id}.pdb"
    print(f"Downloading {pdb_id} from {url}...")
    try:
        response = requests.get(url)
        response.raise_for_status()
        with open(output_path, "w") as f:
            f.write(response.text)
        print("Download successful.")
    except Exception as e:
        print(f"Failed to download PDB: {e}")

def get_vqivyk_center(mol):
    """
    Calculates the geometric center of residues 306-311 (VQIVYK).
    Since 2MZ7 is the MTBR, the residue numbering might be consistent with the full tau
    or local. We will look for the VQIVYK sequence specifically or rely on the known numbering
    for this PDB structure if possible.

    In 2MZ7, the sequence is typically residues 297-391 of Tau.
    VQIVYK is 306-311.
    """
    target_res_nums = set([306, 307, 308, 309, 310, 311])

    coords = []
    conf = mol.GetConformer()

    for atom in mol.GetAtoms():
        res = atom.GetPDBResidueInfo()
        if res:
            if res.GetResidueNumber() in target_res_nums:
                pos = conf.GetAtomPosition(atom.GetIdx())
                coords.append((pos.x, pos.y, pos.z))

    if not coords:
        print("Warning: VQIVYK residues not found by ID. Using Center of Mass of entire model.")
        coords = [conf.GetAtomPosition(i) for i in range(mol.GetNumAtoms())]

    xs = [c[0] for c in coords]
    ys = [c[1] for c in coords]
    zs = [c[2] for c in coords]

    return (sum(xs)/len(xs), sum(ys)/len(ys), sum(zs)/len(zs))

def split_and_prep(pdb_path, output_dir):
    """
    Splits PDB into models and converts each to PDBQT.
    """
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    with open(pdb_path, "r") as f:
        content = f.read()

    # Split by MODEL/ENDMDL
    # Note: PDB files might not have MODEL tags if only 1 model, but 2MZ7 should.
    models = []
    if "MODEL" in content:
        parts = content.split("MODEL")
        # Part 0 is header, rest are models
        for i, part in enumerate(parts[1:]):
            # Find ENDMDL
            if "ENDMDL" in part:
                model_content = "MODEL" + part.split("ENDMDL")[0] + "ENDMDL"
                # We need to preserve ATOM lines.
                # Actually, RDKit might fail if we just feed it a chunk without header?
                # Usually RDKit handles PDB blocks fine.
                models.append(model_content)
    else:
        models.append(content)

    print(f"Found {len(models)} models.")

    centers = {}

    for i, model_block in enumerate(models):
        model_num = i + 1
        print(f"Processing Model {model_num}...")

        # Load into RDKit
        mol = Chem.MolFromPDBBlock(model_block, removeHs=False)
        if mol is None:
            print(f"Failed to load Model {model_num}")
            continue

        # Add Hydrogens explicitly if missing (though removeHs=False should keep them,
        # but sometimes PDBs lack them or RDKit needs them marked explicit)
        mol = Chem.AddHs(mol, addCoords=True)

        # Calculate Center
        center = get_vqivyk_center(mol)
        centers[str(model_num)] = center

        # Compute Charges
        try:
            AllChem.ComputeGasteigerCharges(mol)
        except Exception as e:
            print(f"Charge computation failed for Model {model_num}: {e}")

        # Convert to PDBQT
        try:
            # Using Meeko's PDBQTReceptor for proteins is often more reliable or using PDBQTWriterLegacy
            # But RDKitMoleculeSetup.from_mol needs a clean mol.
            # Let's try to make sure we don't have sanitization issues.
            try:
                Chem.SanitizeMol(mol)
            except:
                pass

            setup = RDKitMoleculeSetup.from_mol(mol)
            # Legacy writer often works better for proteins than the new writer if we want simple PDBQT
            pdbqt_string, success, msg = PDBQTWriterLegacy.write_string(setup)

            if success:
                out_name = os.path.join(output_dir, f"2MZ7_model_{model_num}.pdbqt")
                with open(out_name, "w") as f:
                    f.write(pdbqt_string)
                print(f"Saved {out_name}")
            else:
                print(f"Meeko failed for Model {model_num}: {msg}. Using fallback.")
                out_name = os.path.join(output_dir, f"2MZ7_model_{model_num}.pdbqt")
                write_pdbqt_manual(mol, out_name)

        except Exception as e:
            print(f"PDBQT conversion failed for Model {model_num}: {e}")
            out_name = os.path.join(output_dir, f"2MZ7_model_{model_num}.pdbqt")
            write_pdbqt_manual(mol, out_name)

    # Save centers
    with open(os.path.join(output_dir, "2MZ7_centers.json"), "w") as f:
        json.dump(centers, f, indent=2)

def main():
    data_dir = "data"
    pdb_name = "2MZ7.pdb"
    pdb_path = os.path.join(data_dir, pdb_name)

    if not os.path.exists(data_dir):
        os.makedirs(data_dir)

    download_pdb("2MZ7", pdb_path)
    split_and_prep(pdb_path, data_dir)

if __name__ == "__main__":
    main()

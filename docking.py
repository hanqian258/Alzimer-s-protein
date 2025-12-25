import os
import subprocess
from rdkit import Chem
from rdkit.Chem import AllChem
from meeko import MoleculePreparation
from meeko import PDBQTMolecule

def prep_ligand(smiles, name="ligand"):
    """
    Converts SMILES to a 3D PDBQT string.
    """
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return None

        # Add Hydrogens
        mol = Chem.AddHs(mol)

        # Generate 3D Conformer
        params = AllChem.ETKDGv3()
        params.useSmallRingTorsions = True
        AllChem.EmbedMolecule(mol, params=params)

        # Optimize Geometry (MMFF)
        try:
            AllChem.MMFFOptimizeMolecule(mol)
        except:
            pass # Sometimes fails, skip

        # Convert to PDBQT using Meeko
        preparator = MoleculePreparation()
        preparator.prepare(mol)
        # In version 0.5+, prepare returns setups but we can also use write_pdbqt_string on the preparator
        # or PDBQTWriterLegacy if needed.
        # But let's try the simplest:
        pdbqt_string = preparator.write_pdbqt_string()
        return pdbqt_string
    except Exception as e:
        print(f"Error preparing ligand {name}: {e}")
        return None

def run_docking(ligand_pdbqt, receptor_path, center, size=(20, 20, 20)):
    """
    Runs AutoDock Vina.
    """
    import tempfile

    # Write ligand to temp file
    with tempfile.NamedTemporaryFile(mode='w', suffix='.pdbqt', delete=False) as tmp_lig:
        tmp_lig.write(ligand_pdbqt)
        tmp_lig_path = tmp_lig.name

    out_pdbqt_path = tmp_lig_path.replace(".pdbqt", "_out.pdbqt")

    # Construct Vina command
    # ./vina --receptor data/receptor.pdbqt --ligand tmp.pdbqt --center_x ... --center_y ... --center_z ... --size_x ... --cpu 1

    # Determine path to vina binary. Assuming it's in the current directory or specified.
    vina_path = os.path.abspath("vina")

    cmd = [
        vina_path,
        "--receptor", receptor_path,
        "--ligand", tmp_lig_path,
        "--center_x", str(center[0]),
        "--center_y", str(center[1]),
        "--center_z", str(center[2]),
        "--size_x", str(size[0]),
        "--size_y", str(size[1]),
        "--size_z", str(size[2]),
        "--cpu", "1",
        "--out", out_pdbqt_path,
        "--verbosity", "0"
    ]

    try:
        # Run Vina
        subprocess.run(cmd, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)

        # Parse output for score
        best_affinity = 0.0

        if os.path.exists(out_pdbqt_path):
            with open(out_pdbqt_path, 'r') as f:
                content = f.read()
                # Find line: REMARK VINA RESULT:   -6.5      0.000      0.000
                for line in content.split('\n'):
                    if "REMARK VINA RESULT:" in line:
                        parts = line.split()
                        if len(parts) >= 4:
                            best_affinity = float(parts[3])
                            break

            # Read the docked poses
            docked_pdbqt = content

            # Clean up
            os.remove(tmp_lig_path)
            os.remove(out_pdbqt_path)

            return best_affinity, docked_pdbqt

        else:
            print("Vina output file not found.")
            return None, None

    except subprocess.CalledProcessError as e:
        print(f"Vina failed: {e}")
        return None, None
    except Exception as e:
        print(f"Docking error: {e}")
        return None, None

if __name__ == "__main__":
    # Test block
    pass

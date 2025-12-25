import os
import subprocess
import shutil
from rdkit import Chem
from rdkit.Chem import AllChem
from meeko import MoleculePreparation

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

        # Handle fragments/salts (keep largest)
        frags = Chem.GetMolFrags(mol, asMols=True, sanitizeFrags=False)
        if len(frags) > 1:
            mol = max(frags, key=lambda m: m.GetNumAtoms())

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
        pdbqt_string = preparator.write_pdbqt_string()
        return pdbqt_string
    except Exception as e:
        print(f"Error preparing ligand {name}: {e}")
        return None

def get_vina_path():
    """
    Attempts to locate the AutoDock Vina executable.
    Prioritizes system PATH, then checks local directory.
    """
    # Check system PATH
    system_vina = shutil.which("vina")
    if system_vina:
        return system_vina

    # Check local directory (legacy/linux support)
    # Also check for common names people might give it after download
    for potential_name in ["vina", "vina_mac", "vina_linux", "vina_1.2.5_mac_x86_64"]:
        local_vina = os.path.abspath(potential_name)
        if os.path.exists(local_vina) and os.access(local_vina, os.X_OK):
            return local_vina

    return None

def run_docking(ligand_pdbqt, receptor_path, center, size=(20, 20, 20)):
    """
    Runs AutoDock Vina.
    """
    import tempfile

    # Determine path to vina binary
    vina_path = get_vina_path()

    if not vina_path:
        print("Error: AutoDock Vina executable not found.")
        print("Please download the Vina executable for your OS and place it in this folder named 'vina'.")
        print("Download link: https://github.com/ccsb-scripps/AutoDock-Vina/releases")
        return None, None

    # Write ligand to temp file
    with tempfile.NamedTemporaryFile(mode='w', suffix='.pdbqt', delete=False) as tmp_lig:
        tmp_lig.write(ligand_pdbqt)
        tmp_lig_path = tmp_lig.name

    out_pdbqt_path = tmp_lig_path.replace(".pdbqt", "_out.pdbqt")

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
        docked_pdbqt = None

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

            docked_pdbqt = content
            os.remove(out_pdbqt_path)
        else:
            print("Vina output file not found.")

        # Clean up input temp file
        if os.path.exists(tmp_lig_path):
            os.remove(tmp_lig_path)

        return best_affinity, docked_pdbqt

    except subprocess.CalledProcessError as e:
        print(f"Vina failed: {e}")
        # Clean up
        if os.path.exists(tmp_lig_path): os.remove(tmp_lig_path)
        return None, None
    except Exception as e:
        print(f"Docking error: {e}")
        if os.path.exists(tmp_lig_path): os.remove(tmp_lig_path)
        return None, None

if __name__ == "__main__":
    pass

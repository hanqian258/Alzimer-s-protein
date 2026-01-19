import os
import subprocess
import shutil
import platform
import stat
import tempfile
import requests
from rdkit import Chem
from rdkit.Chem import AllChem
from meeko import MoleculePreparation
from typing import Optional, Tuple


def prep_ligand(smiles: str, name: str = "ligand") -> Optional[str]:
    """
    Converts SMILES to a 3D PDBQT string.

    Args:
        smiles (str): The SMILES string of the ligand.
        name (str): The name of the ligand for logging/error reporting.

    Returns:
        Optional[str]: The PDBQT string if successful, None otherwise.
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
        except Exception:
            pass  # Sometimes fails, skip

        # Convert to PDBQT using Meeko
        preparator = MoleculePreparation()
        preparator.prepare(mol)
        pdbqt_string = preparator.write_pdbqt_string()
        return pdbqt_string
    except Exception as e:
        print(f"Error preparing ligand {name}: {e}")
        return None


def download_vina() -> Optional[str]:
    """
    Downloads the AutoDock Vina executable for the current OS.

    Returns:
        Optional[str]: The absolute path to the downloaded executable, or None if failed.
    """
    system = platform.system()
    machine = platform.machine()

    base_url = "https://github.com/ccsb-scripps/AutoDock-Vina/releases/download/v1.2.5"
    filename = "vina"
    url = ""

    print(f"Detected system: {system} {machine}")

    if system == "Linux":
        url = f"{base_url}/vina_1.2.5_linux_x86_64"
    elif system == "Darwin":
        url = f"{base_url}/vina_1.2.5_mac_x86_64"
    elif system == "Windows":
        url = f"{base_url}/vina_1.2.5_win.exe"
        filename = "vina.exe"
    else:
        print(f"Unsupported system: {system}")
        return None

    print(f"Downloading Vina from {url}...")
    try:
        response = requests.get(url, stream=True)
        response.raise_for_status()
        with open(filename, 'wb') as f:
            for chunk in response.iter_content(chunk_size=8192):
                f.write(chunk)

        if system != "Windows":
            st_mode = os.stat(filename)
            os.chmod(filename, st_mode.st_mode | stat.S_IEXEC)

        print("Vina downloaded successfully.")
        return os.path.abspath(filename)
    except Exception as e:
        print(f"Failed to download Vina: {e}")
        return None


def get_vina_path() -> Optional[str]:
    """
    Attempts to locate the AutoDock Vina executable.
    Prioritizes system PATH, then checks local directory.
    If not found, attempts to download it.

    Returns:
        Optional[str]: The path to the Vina executable.
    """
    # Check system PATH
    system_vina = shutil.which("vina")
    if system_vina:
        return system_vina

    # Check local directory (legacy/linux support)
    # Also check for common names people might give it after download
    for potential_name in ["vina", "vina_mac", "vina_linux", "vina_1.2.5_mac_x86_64", "vina.exe"]:
        local_vina = os.path.abspath(potential_name)
        if os.path.exists(local_vina) and os.access(local_vina, os.X_OK):
            return local_vina

    # Try downloading if not found
    return download_vina()


def run_docking(
    ligand_pdbqt: str,
    receptor_path: str,
    center: Tuple[float, float, float],
    size: Tuple[float, float, float] = (20, 20, 20)
) -> Tuple[Optional[float], Optional[str]]:
    """
    Runs AutoDock Vina to dock the ligand against the receptor.

    Args:
        ligand_pdbqt (str): The ligand structure in PDBQT format.
        receptor_path (str): The path to the receptor PDBQT file.
        center (Tuple[float, float, float]): The center of the docking grid (x, y, z).
        size (Tuple[float, float, float]): The size of the docking grid (x, y, z).

    Returns:
        Tuple[Optional[float], Optional[str]]: A tuple containing the best binding affinity (score)
                                              and the docked PDBQT string.
                                              Returns (None, None) on failure.
    """
    # Determine path to vina binary
    vina_path = get_vina_path()

    if not vina_path:
        print("Error: AutoDock Vina executable not found and download failed.")
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
        if os.path.exists(tmp_lig_path):
            os.remove(tmp_lig_path)
        return None, None
    except Exception as e:
        print(f"Docking error: {e}")
        if os.path.exists(tmp_lig_path):
            os.remove(tmp_lig_path)
        return None, None

def run_docking_on_list(ligand_pdbqt, receptor_configs):
    """
    Runs docking against multiple receptor configurations and returns the average score.

    receptor_configs: List of dicts, e.g.,
    [
      {'path': 'data/rec1.pdbqt', 'center': (0,0,0)},
      {'path': 'data/rec2.pdbqt', 'center': (1,1,1)}
    ]
    """
    scores = []
    poses = []

    for config in receptor_configs:
        path = config.get('path')
        center = config.get('center')
        if not path or not center:
            continue

        score, pose = run_docking(ligand_pdbqt, path, center)
        if score is not None:
            scores.append(score)
            poses.append(pose)

    if not scores:
        return None, None

    avg_score = sum(scores) / len(scores)
    # Return the best pose (lowest energy) as the representative pose
    best_idx = scores.index(min(scores))
    best_pose = poses[best_idx]

    return avg_score, best_pose

if __name__ == "__main__":
    pass

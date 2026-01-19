import math
from typing import Tuple, List, Optional

def calculate_kd(delta_g: float, temperature: float = 298.15) -> float:
    """
    Calculates the dissociation constant (Kd) from Gibbs free energy (delta_g).

    Args:
        delta_g (float): Binding affinity in kcal/mol.
        temperature (float): Temperature in Kelvin. Default is 298.15 K.

    Returns:
        float: Kd in micromolar (uM).
    """
    R = 0.0019872042586493  # Gas constant in kcal/(K*mol)
    try:
        # Kd = exp(delta_g / (R * T))
        kd_molar = math.exp(delta_g / (R * temperature))
        kd_micromolar = kd_molar * 1e6
        return kd_micromolar
    except OverflowError:
        return float('inf')


def interpret_affinity(delta_g: float) -> str:
    """
    Provides a qualitative interpretation of the binding affinity.

    Args:
        delta_g (float): Binding affinity in kcal/mol.

    Returns:
        str: Qualitative description (e.g., "Strong", "Moderate", "Weak").
    """
    if delta_g <= -9.0:
        return "Very Strong"
    elif delta_g <= -7.0:
        return "Strong"
    elif delta_g <= -5.5:
        return "Moderate"
    else:
        return "Weak"


def get_score_color(delta_g: float) -> str:
    """
    Returns a hex color code based on the binding affinity.

    Args:
        delta_g (float): Binding affinity in kcal/mol.

    Returns:
        str: Hex color string.
    """
    if delta_g <= -9.0:
        return "#006400"  # DarkGreen
    elif delta_g <= -7.0:
        return "#228B22"  # ForestGreen
    elif delta_g <= -5.5:
        return "#DAA520"  # GoldenRod
    else:
        return "#B22222"  # FireBrick

def get_ligand_centroid(pdbqt_string: str) -> Optional[Tuple[float, float, float]]:
    """
    Parses a PDBQT string and calculates the geometric centroid of the ligand.

    Args:
        pdbqt_string (str): The content of the PDBQT file.

    Returns:
        Optional[Tuple[float, float, float]]: The (x, y, z) centroid, or None if no atoms found.
    """
    x_sum = 0.0
    y_sum = 0.0
    z_sum = 0.0
    count = 0

    lines = pdbqt_string.splitlines()
    for line in lines:
        if line.startswith("ATOM") or line.startswith("HETATM"):
            try:
                # PDB format: X is 30-38, Y is 38-46, Z is 46-54
                x = float(line[30:38])
                y = float(line[38:46])
                z = float(line[46:54])

                x_sum += x
                y_sum += y
                z_sum += z
                count += 1
            except ValueError:
                continue

    if count == 0:
        return None

    return (x_sum / count, y_sum / count, z_sum / count)

def calculate_distance(p1: Tuple[float, float, float], p2: Tuple[float, float, float]) -> float:
    """
    Calculates Euclidean distance between two 3D points.

    Args:
        p1: (x, y, z) tuple
        p2: (x, y, z) tuple

    Returns:
        float: The distance.
    """
    return math.sqrt((p1[0] - p2[0])**2 + (p1[1] - p2[1])**2 + (p1[2] - p2[2])**2)

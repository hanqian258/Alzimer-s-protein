import math

def calculate_kd(delta_g, temperature=298.15):
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

def interpret_affinity(delta_g):
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

def get_score_color(delta_g):
    """
    Returns a hex color code based on the binding affinity.

    Args:
        delta_g (float): Binding affinity in kcal/mol.

    Returns:
        str: Hex color string.
    """
    if delta_g <= -9.0:
        return "#006400" # DarkGreen
    elif delta_g <= -7.0:
        return "#228B22" # ForestGreen
    elif delta_g <= -5.5:
        return "#DAA520" # GoldenRod
    else:
        return "#B22222" # FireBrick

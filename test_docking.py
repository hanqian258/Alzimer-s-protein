import docking
import os

def test():
    # Load config
    with open("data/config.txt") as f:
        c = f.read().strip().split(',')
        center = (float(c[0]), float(c[1]), float(c[2]))

    print(f"Testing docking with center: {center}")

    # Simple molecule: Benzene
    smiles = "c1ccccc1"

    print("Preparing ligand...")
    pdbqt = docking.prep_ligand(smiles, "Benzene")

    if pdbqt:
        print("Ligand prepared. Running Vina...")
        score, docked_pose = docking.run_docking(pdbqt, "data/receptor.pdbqt", center)
        print(f"Docking Score: {score}")
        if score < 0:
            print("SUCCESS: Docking ran and produced a negative binding energy.")
        else:
            print("FAILURE: Docking score suspicious (>=0).")
    else:
        print("Ligand preparation failed.")

if __name__ == "__main__":
    test()

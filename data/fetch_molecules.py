import pandas as pd
import requests
import time
import os

# Define lists of molecules
fda_drugs = [
    "Donepezil", "Memantine", "Galantamine", "Rivastigmine", "Methylene Blue",
    "Aspirin", "Ibuprofen", "Acetaminophen", "Metformin", "Atorvastatin",
    "Lisinopril", "Amlodipine", "Metoprolol", "Omeprazole", "Losartan",
    "Simvastatin", "Levothyroxine", "Gabapentin", "Hydrochlorothiazide", "Sertraline",
    "Furosemide", "Fluticasone", "Amoxicillin", "Alprazolam", "Atenolol",
    "Citalopram", "Escitalopram", "Bupropion", "Trazodone", "Tramadol",
    "Duloxetine", "Prednisone", "Montelukast", "Rosuvastatin", "Pravastatin",
    "Carvedilol", "Pantoprazole", "Meloxicam", "Clopidogrel", "Cyclobenzaprine",
    "Methylphenidate", "Zolpidem", "Warfarin", "Venlafaxine", "Clonazepam",
    "Fluoxetine", "Lorazepam", "Celecoxib", "Naproxen", "Doxycycline",
    "Allopurinol", "Amitriptyline", "Risperidone", "Levodopa", "Carbidopa",
    "Entacapone", "Selegiline", "Rasagiline", "Ropinirole", "Pramipexole"
]

flavonoids = [
    "Quercetin", "Kaempferol", "Myricetin", "Fisetin", "Galangin",
    "Isorhamnetin", "Luteolin", "Apigenin", "Baicalein", "Chrysin",
    "Scutellarein", "Wogonin", "Hesperetin", "Naringenin", "Eriodictyol",
    "Homoeriodictyol", "Taxifolin", "Dihydrokaempferol", "Dihydroquercetin", "Dihydromyricetin",
    "Catechin", "Epicatechin", "Gallocatechin", "Epigallocatechin", "Epicatechin gallate",
    "Epigallocatechin gallate", "Theaflavin", "Procyanidin B1", "Procyanidin B2", "Genistein",
    "Daidzein", "Glycitein", "Formononetin", "Biochanin A", "Cyanidin",
    "Delphinidin", "Pelargonidin", "Malvidin", "Peonidin", "Petunidin",
    "Phloretin", "Isoliquiritigenin", "Butein", "Tangeretin", "Nobiletin",
    "Sinensetin", "Rutin", "Hesperidin", "Naringin", "Morin",
    "Baicalin", "Wogonoside", "Liquiritigenin", "Glycyrrhizin", "Silibinin"
]

known_ligands = [
    "Methylene Blue", "Epothilone D", "PBB3", "T807", "Anandamide",
    "Surfen", "Levistolide A", "Leucomethylene blue", "C004019", "RI-AG03",
    "AChE/BChE-IN-29", "QC-01-175", "TTBK1-IN-1", "Verbenalin", "Ceperognastat",
    "PhosTAC7", "Simufilam", "Curcumin", "EGCG", "TRx0237",
    "LMTM", "Tideglusib", "Lithium carbonate", "Valproate", "Tannic acid",
    "Oleocanthal", "Resveratrol", "Berberine", "Palmatine", "Jatrorrhizine",
    "Coptisine", "Epiberberine", "Emodin", "Daunorubicin", "Mitoxantrone",
    "Pixantrone", "Epalrestat", "Troglitazone", "C11"
]

# Manual SMILES for compounds not easily found by name in PUG REST or needing specific forms
manual_smiles = {
    "AQ2S": "C1=CC=C2C(=C1)C(=O)C3=C(C2=O)C=C(C=C3)S(=O)(=O)O"
}

# Add manual entries to known_ligands if they are keys in manual_smiles
known_ligands.extend(list(manual_smiles.keys()))

# Map to categories in the CSV
categories = {
    "FDA Drug": fda_drugs,
    "Flavonoid": flavonoids,
    "Known Ligand": known_ligands
}

def get_smiles(name):
    """Fetches Isomeric SMILES from PubChem PUG REST API."""
    if name in manual_smiles:
        return manual_smiles[name]

    try:
        # Handle specific cases like C11 which is a common name
        search_name = name
        if name == "C11":
             search_name = "3,3'-Diethyl-9-methylthiacarbocyanine iodide"

        url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{search_name}/property/IsomericSMILES/TXT"
        response = requests.get(url, timeout=10)
        if response.status_code == 200:
            return response.text.strip()
        else:
            # Try removing spaces or specific salt forms if 404?
            # For now just return None
            print(f"Failed to fetch SMILES for {name}: Status {response.status_code}")
            return None
    except Exception as e:
        print(f"Error fetching SMILES for {name}: {e}")
        return None

def main():
    csv_path = "data/molecules.csv"

    # Load existing data
    if os.path.exists(csv_path):
        df = pd.read_csv(csv_path)
    else:
        df = pd.DataFrame(columns=["name", "category", "smiles", "citation"])

    # Pre-clean T807 duplication if exists
    # Keep "T807 (AV-1451)" and remove "T807" if both present
    if "T807 (AV-1451)" in df['name'].values and "T807" in df['name'].values:
         print("Removing duplicate T807 entry...")
         df = df[df['name'] != "T807"]

    existing_names = set(df['name'].str.lower())

    # Special handling for T807
    if "t807 (av-1451)" in existing_names:
        existing_names.add("t807")

    new_rows = []

    for category, names in categories.items():
        print(f"Processing category: {category}")
        for name in names:
            # Check for duplicates (case-insensitive)
            if name.lower() in existing_names:
                print(f"Skipping {name} (already exists)")
                continue

            # Fetch SMILES
            print(f"Fetching SMILES for {name}...")
            smiles = get_smiles(name)

            if smiles:
                new_rows.append({
                    "name": name,
                    "category": category,
                    "smiles": smiles,
                    "citation": f"https://pubchem.ncbi.nlm.nih.gov/compound/{name}"
                })
                # Respect API rate limits
                time.sleep(0.5)
            else:
                print(f"Could not find SMILES for {name}")

    if new_rows:
        new_df = pd.DataFrame(new_rows)
        # Append to existing
        df = pd.concat([df, new_df], ignore_index=True)

        # Final Deduplication on name and category just in case
        df = df.drop_duplicates(subset=['name', 'category'])

        df.to_csv(csv_path, index=False)
        print(f"Successfully added {len(new_rows)} new molecules to {csv_path}")
    else:
        print("No new molecules added.")

if __name__ == "__main__":
    main()

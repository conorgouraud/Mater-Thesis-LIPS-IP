import pandas as pd
import ast
import sys 

class KMDReferenceEntry:
    def __init__(self, general_abbreviation, kmd, sub_class):
        self.general_abbreviation = general_abbreviation
        self.kmd = kmd
        self.sub_class = sub_class

def calculate_kendrick_mass(mass, base_mass=14.01565):
    """Calculate the Kendrick Mass (KM) for a given mass."""
    return (mass / base_mass) * 14

def kendrick_mass_defect(kendrick_mass):
    """Calculate the Kendrick Mass Defect (KMD) for a given KM."""
    return kendrick_mass - int(kendrick_mass)

def sanitize_mass_string(mass_string):
    """Sanitize the mass string to remove unwanted characters and convert to float."""
    try:
        mass_string = mass_string.strip().replace("'", "").replace(")", "").replace("(", "")
        return float(mass_string)
    except ValueError:
        print(f"Could not convert to float: {mass_string}")
        return None

def load_reference_database(reference_database_path):
    """Load the reference database into a list of KMDReferenceEntry objects."""
    df = pd.read_csv(reference_database_path)
    reference_entries = []

    for _, row in df.iterrows():
        entry = KMDReferenceEntry(
            general_abbreviation=row['GENERAL_ABBREVIATION'],
            kmd=row['KMD'],
            sub_class=row['SUB_CLASS']
        )
        reference_entries.append(entry)

    return reference_entries

def find_kmd_matches_actual(mass_list, reference_entries):
    """Find KMD matches for a given list of masses against the reference entries."""
    results = {}

    for entry in mass_list:
        mass_part = entry.split('(')[-1].strip(')')
        mass = sanitize_mass_string(mass_part)

        if mass is None:
            continue

        km = calculate_kendrick_mass(mass)
        kmd = kendrick_mass_defect(km)

        matches = []

        for ref_entry in reference_entries:
            difference = kmd - ref_entry.kmd
            result = difference / 0.0134
            closest_integer = round(result)

            if -16 < closest_integer < 0.5 and abs(result - closest_integer) <= 0.0075:
                matches.append((ref_entry.general_abbreviation, abs(closest_integer), ref_entry.sub_class))

        results[entry] = matches

    return results

def find_kmd_matches(row, index, reference_entries):
    """Convert the string representation of a list to an actual list and find KMD matches."""
    #try:
    formula_list = ast.literal_eval(row)
    '''except (ValueError, SyntaxError):
        print(f"Could not parse row:")
        return {}'''

    result = find_kmd_matches_actual(formula_list, reference_entries)
    
    return result







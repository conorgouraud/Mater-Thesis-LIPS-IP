import pandas as pd
import ast
import re


# Function to parse the headgroup, DB, and additional O
def parse_entry(a, b):
    if not isinstance(a, str):
        return None, None, None
    headgroup_match = re.search(r'^(.*?O-)' if 'O-' in a else r'^(.*?)(?=\s*X|:|;|$)', a)
    headgroup = headgroup_match.group(1).strip() if headgroup_match else ''
    db_match = re.search(r':(\d+)', a)
    db = int(db_match.group(1)) + b if db_match else b
    additional_o_match = re.search(r';O(\d*)', a)
    additional_o = int(additional_o_match.group(1)) if additional_o_match and additional_o_match.group(1) else 0 if additional_o_match else 0

    return headgroup, db, additional_o

# Function to convert a chemical formula to a dictionary
def formula_to_dict(formula):
    matches = re.findall(r'([A-Z][a-z]*)(\d*)', formula)
    return {element: int(count) if count else 1 for (element, count) in matches}

def subtract_formula_and_validate(key_formula, chem_formula, chains, db):
    key_dict, chem_dict = formula_to_dict(key_formula), formula_to_dict(chem_formula)
    for elem in ['N', 'O', 'P', 'S']:
        if key_dict.get(elem, 0) != chem_dict.get(elem, 0):
            return None, elem
    if key_dict['C'] >= chem_dict.get('C', 0) and key_dict['H'] >= chem_dict.get('H', 0):
        remaining_C, remaining_H = key_dict['C'] - chem_dict.get('C', 0), key_dict['H'] - chem_dict.get('H', 0) - chains
        if remaining_H == 2 * remaining_C - 2 * db:
            return {'C': remaining_C, 'H': remaining_H}, None
        else:
            return {'C': remaining_C, 'H': remaining_H}, 'formula_mismatch'
    else:
        return None, 'invalid_subtraction'

# Function to subtract the formula, chains, and validate C and H
def check_heads(df, headgroup_db):

    df['Non-Empty Keys Count'] = 0
    df['Accepted Keys Count'] = 0
    df['Unique Lipids Count'] = 0
    df['Unique Lipids List'] = None
    
    for index, row in df.iterrows():
        row_dict = row['Formula:subclass combos']
        if not row_dict:
            continue
        
        non_empty_keys_count = sum(1 for k, v in row_dict.items() if v)
        accepted_keys = set()
        unique_lipids = set()
        
        for key, entries in row_dict.items():
            if entries:
                for entry in entries:
                    a, b, c = entry
                    headgroup, db, additional_o = parse_entry(a, b)
                    if headgroup is None:
                        continue

                    match = headgroup_db[(headgroup_db['class'] == headgroup) & (headgroup_db['additional O'] == additional_o)]
                    if not match.empty:
                        chem_formula, chains = match.iloc[0]['chem'], match.iloc[0]['chains']
                        remaining_elements, issue = subtract_formula_and_validate(key, chem_formula, chains, db)
                        if remaining_elements and issue is None:
                            accepted_keys.add(key)
                            final_c = remaining_elements['C'] + 1
                            final_lipid = f"{headgroup} {final_c}:{db};O{additional_o}" if additional_o > 0 else f"{headgroup} {final_c}:{db}"
                            unique_lipids.add(final_lipid)
        
        df.at[index, 'Non-Empty Keys Count'] = non_empty_keys_count
        df.at[index, 'Accepted Keys Count'] = len(accepted_keys)
        df.at[index, 'Unique Lipids Count'] = len(unique_lipids)
        df.at[index, 'Unique Lipids List'] = list(unique_lipids)

    return df




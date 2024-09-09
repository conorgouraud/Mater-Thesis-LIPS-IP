from nominal_mass_pred import *

# Function to extract element count from the formula
def extract_element_count(formula, element):
    element_index = formula.find(element)
    if element_index == -1:
        return 0
    number_start = element_index + len(element)
    number_end = number_start

    while number_end < len(formula) and formula[number_end].isdigit():
        number_end += 1

    return int(formula[number_start:number_end] or 0)


def apply_rules(row, formula_column):
    nominal_mass = row['nominal_mass']
    formulas = row[formula_column]
    valid_formulas = []

    for formula in formulas:
        h_count = extract_element_count(formula, 'H')
        br_count = extract_element_count(formula, 'Br')
        cl_count = extract_element_count(formula, 'Cl')
        n_count = extract_element_count(formula, 'N')
        p_count = extract_element_count(formula, 'P')

        '''# Check the first nitrogen rule
        if (nominal_mass % 2 == 0 and n_count % 2 != 0) or (nominal_mass % 2 != 0 and n_count % 2 == 0):
            continue  # Skip this formula if it doesn't satisfy the nitrogen rule'''
        
        result = (31 * p_count + h_count + 35 * cl_count + 79 * br_count) % 4 + nominal_mass % 4
        
        if result % 4 == 0:
            valid_formulas.append(formula)
    
    return valid_formulas

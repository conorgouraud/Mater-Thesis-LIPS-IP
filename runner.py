import numpy as np
from pyopenms import MSExperiment, MzMLFile, MzXMLFile
import pandas as pd
import os
import time
import pickle
from tqdm import tqdm
import csv
from rpy2.robjects import pandas2ri, r
from rpy2.robjects.packages import importr
import rpy2.robjects as ro
import subprocess

from peak_picking import generate_peak_lists_for_all_rt, peaks_to_networks_to_chains, is_in_peak, Peak
from centwave import centwave_on_raw_data
from nominal_mass_pred import PredNomMass, k, MSE, corr_fact, diff_mass_sum
from element_range_pred import *
from refined_hydrogen_rule import *
from Lipid_class_predictor import *

start_time = time.time()

################################################################################################################################
##### Main Parameters #####
###########################

raw_file_path = 'LIPIDS3_20231027_011.mzXML'

pred_interval = 0.995 # for the RR prediction
charge = 1 # +1 or -1 for positive or negative charge adduct
mass_error = 5  # in ppm (+/- mass*ppm/1e6)

################################################################################################################################

################################
##### STEP 1: Peak Picking #####
################################

exp = MSExperiment()
MzXMLFile().load(raw_file_path, exp)

### Getting triples from raw data ###
all_peak_lists = generate_peak_lists_for_all_rt(exp)

all_triples_dictionary = {}
for peak_list in tqdm(all_peak_lists, desc="Mining Triples...", unit="list"):
    if peak_list:
        rt = peak_list[0]['rtime']
        rt = rt

        arrays_of_triples = peaks_to_networks_to_chains(peak_list, mass_error)
        all_triples_dictionary[rt] = [array for array in arrays_of_triples if array.size > 0]

# Saving a dictionary with all isotope patterns before centwave filtering
with open(f'trackable_outputs/all_triples_dictionary_{raw_file_path}.pkl', 'wb') as f:
    pickle.dump(all_triples_dictionary, f)

print('There are', sum(len(items) for items in all_triples_dictionary.values()), 'triples before centwave.')

### Centwave part on raw data to get EICs####
if __name__ == "__main__":
    print('Getting all EICs using centWave...')
    output_file = "trackable_outputs/centwave_output.csv"
    centwave_on_raw_data(raw_file_path, output_file)
    print('Finished using Centwave.')

### Checking if mono isotopic mass of triple lies in a EIC from centWave (bad isotope patterns removed) ###
centwave_df = pd.read_csv('trackable_outputs/centwave_output.csv') # not numeric for some reason so change like this for the moment
centwave_df['rt'] = pd.to_numeric(centwave_df['rt'], errors='coerce') 
centwave_df['rtmin'] = pd.to_numeric(centwave_df['rtmin'], errors='coerce') 
centwave_df['rtmax'] = pd.to_numeric(centwave_df['rtmax'], errors='coerce') 
centwave_df['mzmin'] = pd.to_numeric(centwave_df['mzmin'], errors='coerce')
centwave_df['mzmax'] = pd.to_numeric(centwave_df['mzmax'], errors='coerce')
centwave_df['mz'] = pd.to_numeric(centwave_df['mz'], errors='coerce')

peaks = [Peak(row['rt'], row['rtmin'], row['rtmax'], row['mzmin'], row['mzmax'], row['mz'], idx) for idx, row in centwave_df.iterrows()]

not_good_triples = []
good_triple_dict = {}
for rt, triples in tqdm(all_triples_dictionary.items(), desc="Checking if triples lie in a good quality EIC..."):
    good_triple_dict[rt] = []
    for triple in triples:
        peak_id, peak_mz, peak_rt_min, peak_rt_max, peak_rt = is_in_peak(rt, triple[0][0], peaks)  # Check each triple
        if peak_id is not None:
            good_triple_dict[rt].append((triple, peak_id, peak_mz, peak_rt_min, peak_rt_max, peak_rt)) 
        else:
            not_good_triples.append((rt, triple[0][0]))

# Saving the 'bad' triples as list, these are not used anymore
with open(f'trackable_outputs/not_good_triples_dictionary_{raw_file_path}.pkl', 'wb') as file:
    pickle.dump(not_good_triples, file)

# Saving the 'good' triples as dictionary, these are what we will continue with
with open(f'trackable_outputs/good_triples_dictionary_{raw_file_path}.pkl', 'wb') as file:
    pickle.dump(good_triple_dict, file)

print('Good triples found and saved with corresponding EIC IDs.')
print('There are', sum(len(items) for items in good_triple_dict.values()), 'triples whose monoisotopic peak lies in a "good" peak range.')


###########################################
#### STEP 2: Relative Ratio Prediction ####
###########################################

print('Predicting relative ratios and element ranges...')

# Extracting all necessary features from triples
triples_and_their_features = []
for rt, triples in good_triple_dict.items():
    for triple in triples:
        # isotope pattern features
        mass1 = triple[0][0][0] + charge * 1.00727647  # FOR [M+H]+ or [M-H]-
        mass2 = triple[0][1][0] + charge * 1.00727647
        mass3 = triple[0][2][0] + charge * 1.00727647

        wmean_mz = triple[2] + charge * 1.00727647

        m21 = mass2 - mass1
        m32 = mass3 - mass2
        m31 = mass3 - mass1
        abu1 = triple[0][0][1]
        abu2 = triple[0][1][1]
        abu3 = triple[0][2][1]

        nom_mass = PredNomMass(mass1, k, MSE, corr_fact, diff_mass_sum)[0]
        frac_mass = mass1 - nom_mass

        triples_and_their_features.append({
            'RT': rt, 'mass1': mass1, 'mass2': mass2, 'mass3': mass3, 
            'm21': m21, 'm32': m32, 'm31': m31, 'abu1': abu1, 'abu2': abu2, 
            'abu3': abu3, 'nom_mass': nom_mass, 'frac_mass': frac_mass, 
            'EIC_ID': triple[1], 'wmean_mz': wmean_mz, 'rt_min': triple[3], 
            'rt_max': triple[4], 'rt': triple[5]
        })

iso_pattern_df = pd.DataFrame(triples_and_their_features)
iso_pattern_df = iso_pattern_df.sort_values(by='mass1')
iso_pattern_df['m0'] = iso_pattern_df['mass1'] / 1000  # divide by 1000 as model fitted using this
iso_pattern_df.to_csv('trackable_outputs/good_triples_and_their_features.csv.zip', compression='zip', index=False)

RR_lookup_table = pd.read_csv("database_folder/RR_lookup_table.csv") # RR reference table

pandas2ri.activate()  # switch easily between R and pandas

# Predicting C ranges
element = "C"
print("Predicting", element, "ranges...")
readRDS = ro.r['readRDS']
fModel = readRDS("relative_ratios_no_fracmass/resultRR21.noS.rds")

toPred = ro.r['data.frame'](
    m0=iso_pattern_df['m0'],
    m21=iso_pattern_df['m21'],
    m32=iso_pattern_df['m32']
)

resRR = predRR21_agg_nfc(x=toPred, modelOutput=fModel, conf_int=True, conf_level=pred_interval, S=False)
resRR_df = pandas2ri.rpy2py(resRR)


resRR_df = calculate_element_min_max(resRR_df, RR_lookup_table, element, pred_interval)

iso_pattern_df = iso_pattern_df.reset_index(drop=True)  # was some indexing issues so leave these for the moment
resRR_df = resRR_df.reset_index(drop=True)
iso_pattern_df = iso_pattern_df.merge(resRR_df[[f"{element}minpred_{pred_interval}", f"{element}maxpred_{pred_interval}"]], left_index=True, right_index=True)

# Predicting H ranges
element = "H"
print("Predicting", element, "ranges...")
readRDS = ro.r['readRDS']
fModel = readRDS("update_relative_ratios/resultRR21_fractionalMass.noS.rds")

toPred = ro.r['data.frame'](  # used as input to model for predictions
    m0=iso_pattern_df['m0'],
    m21=iso_pattern_df['m21'],
    m32=iso_pattern_df['m32'],
    fractionalMass=iso_pattern_df['frac_mass']
)

resRR = predRR21_agg_fc(x=toPred, modelOutput=fModel, conf_int=True, conf_level=pred_interval, S=False)
resRR_df = pandas2ri.rpy2py(resRR)

resRR_df = calculate_element_min_max(resRR_df, RR_lookup_table, element, pred_interval)

iso_pattern_df = iso_pattern_df.reset_index(drop=True)  # was some indexing issues so leave these for the moment
resRR_df = resRR_df.reset_index(drop=True)
iso_pattern_df = iso_pattern_df.merge(resRR_df[[f"{element}minpred_{pred_interval}", f"{element}maxpred_{pred_interval}"]], left_index=True, right_index=True)

# Predicting N ranges
element = "N"
print("Predicting", element, "ranges...")
readRDS = ro.r['readRDS']
fModel = readRDS("relative_ratios_no_fracmass/resultRR21.noS.rds")

toPred = ro.r['data.frame'](  # used as input to model for predictions
    m0=iso_pattern_df['m0'],
    m21=iso_pattern_df['m21'],
    m32=iso_pattern_df['m32']
)

resRR = predRR21_agg_nfc(x=toPred, modelOutput=fModel, conf_int=True, conf_level=pred_interval, S=False)
resRR_df = pandas2ri.rpy2py(resRR)

resRR_df = calculate_element_min_max(resRR_df, RR_lookup_table, element, pred_interval)

iso_pattern_df = iso_pattern_df.reset_index(drop=True)  # was some indexing issues so leave these for the moment
resRR_df = resRR_df.reset_index(drop=True)
iso_pattern_df = iso_pattern_df.merge(resRR_df[[f"{element}minpred_{pred_interval}", f"{element}maxpred_{pred_interval}"]], left_index=True, right_index=True)

# Predicting O ranges
element = "O"
print("Predicting", element, "ranges...")
readRDS = ro.r['readRDS']
fModel = readRDS("update_relative_ratios/resultRR43_fractionalMass.noS.rds")

toPred = ro.r['data.frame'](  # used as input to model for predictions
    m0=iso_pattern_df['m0'],
    m21=iso_pattern_df['m21'],
    m32=iso_pattern_df['m32'],
    fractionalMass=iso_pattern_df['frac_mass']
)

resRR = predRR43_agg_fc(x=toPred, modelOutput=fModel, conf_int=True, conf_level=pred_interval, S=False)
resRR_df = pandas2ri.rpy2py(resRR)

resRR_df = calculate_element_min_max(resRR_df, RR_lookup_table, element, pred_interval)

iso_pattern_df = iso_pattern_df.reset_index(drop=True)  # was some indexing issues so leave these for the moment
resRR_df = resRR_df.reset_index(drop=True)
iso_pattern_df = iso_pattern_df.merge(resRR_df[[f"{element}minpred_{pred_interval}", f"{element}maxpred_{pred_interval}"]], left_index=True, right_index=True)

# Predicting P ranges
element = "P"
print("Predicting", element, "ranges...")
readRDS = ro.r['readRDS']
fModel = readRDS("update_relative_ratios/resultRR32_fractionalMass.noS.rds")

toPred = ro.r['data.frame'](  # used as input to model for predictions
    m0=iso_pattern_df['m0'],
    m21=iso_pattern_df['m21'],
    m32=iso_pattern_df['m32'],
    fractionalMass=iso_pattern_df['frac_mass']
)

resRR = predRR32_agg_fc(x=toPred, modelOutput=fModel, conf_int=True, conf_level=pred_interval, S=False)
resRR_df = pandas2ri.rpy2py(resRR)

resRR_df = calculate_element_min_max(resRR_df, RR_lookup_table, element, pred_interval)

iso_pattern_df = iso_pattern_df.reset_index(drop=True)  # was some indexing issues so leave these for the moment
resRR_df = resRR_df.reset_index(drop=True)
iso_pattern_df = iso_pattern_df.merge(resRR_df[[f"{element}minpred_{pred_interval}", f"{element}maxpred_{pred_interval}"]], left_index=True, right_index=True)


## Prep for mass decomposition ##
iso_pattern_df['m0'] = iso_pattern_df['m0']*1000  # masses were /1000 for RR prediction

for col in iso_pattern_df.columns[-10:]:
    iso_pattern_df[col] = iso_pattern_df[col].astype(int).astype(str)  # all the element min/max need to have no '.0' for the mass decomp command

iso_pattern_df.to_csv('trackable_outputs/element_range_predictions_per_individual_triple.csv.zip', compression='zip', index=False)
print('Element ranges predicted.')

# getting min of min and max of max for each element range for each EIC and using wmean calculated by centwWave for m/z from now on 
iso_pattern_df_grouped = iso_pattern_df.groupby('EIC_ID').agg(
        Cminpred_995_min=('Cminpred_0.995', 'min'),
        Cmaxpred_995_max=('Cmaxpred_0.995', 'max'),
        Hminpred_995_min=('Hminpred_0.995', 'min'),
        Hmaxpred_995_max=('Hmaxpred_0.995', 'max'),
        Nminpred_995_min=('Nminpred_0.995', 'min'),
        Nmaxpred_995_max=('Nmaxpred_0.995', 'max'),
        Ominpred_995_min=('Ominpred_0.995', 'min'),
        Omaxpred_995_max=('Omaxpred_0.995', 'max'),
        Pminpred_995_min=('Pminpred_0.995', 'min'),
        Pmaxpred_995_max=('Pmaxpred_0.995', 'max'),
        wmean_mz=('wmean_mz', 'first'),
        nom_mass=('nom_mass', 'first'),
        RT_min=('rt_min', 'first'),
        RT_max=('rt_max', 'first'),
        RT=('rt', 'first')
    ).reset_index()

RT_list = iso_pattern_df_grouped['RT']
RTmin_list = iso_pattern_df_grouped['RT_min']
RTmax_list = iso_pattern_df_grouped['RT_max']

iso_pattern_df_grouped.to_csv('trackable_outputs/element_range_predictions_per_EIC.csv.zip', compression='zip',index=False) # ranges saved after triples gouped by EIC

print("After grouping by EIC, there are", len(iso_pattern_df_grouped), "features.")

####################################
#### STEP 3: Mass Decomposition ####
####################################


original_working_directory = os.getcwd() # saving this working directory
os.chdir('mass_decomp_background_files/tools')  # directory with all mass decompositions files just for this bit

number_of_outputs = 999999999  # limit on number of decomps before stopping

# lists for final df
monoisotopic_masses_list = []
num_of_decomps_list = []
possible_formulas_list = []

for index, row in tqdm(iso_pattern_df_grouped.iterrows(), total=len(iso_pattern_df_grouped), desc="Mass Decompositions..."):
    mass = row['wmean_mz']
    cminpred = row['Cminpred_995_min']
    cmaxpred = row['Cmaxpred_995_max']
    hminpred = row['Hminpred_995_min']
    hmaxpred = row['Hmaxpred_995_max']
    nminpred = row['Nminpred_995_min']
    nmaxpred = row['Nmaxpred_995_max']
    ominpred = row['Ominpred_995_min']
    omaxpred = row['Omaxpred_995_max']
    #sminpred = row['Sminpred_995_min']
    #smaxpred = row['Smaxpred_995_max']
    pminpred = row['Pminpred_995_min']
    pmaxpred = row['Pmaxpred_995_max']

    command = [
        './imsdecomp',
        '-a', '/Users/conorgouraud/Documents/Thesis Code/LipidMass/mass_decomp_background_files/res/atom-mono.masses',
        '-m', f'{mass}',
        '-b', f'C{cminpred}H{hminpred}N{nminpred}O{ominpred}P{pminpred}',
        '-g', f'C{cmaxpred}H{hmaxpred}N{nmaxpred}O{omaxpred}P{pmaxpred}',
        '-o', 'mono',
        '-s',
        '-x', f'{number_of_outputs}',
        '-u', 'ppm',
        '-e', f'{mass_error}'
    ]

    result = subprocess.run(command, capture_output=True, text=True)
    output = result.stdout
    lines = output.split('\n')
    num_of_decomps = len(lines[30:len(lines)-3]) # these two lines are specific to only using CHNOPS as input elements for mass decomp
    possible_formulas = lines[30:len(lines)-3]

    monoisotopic_masses_list.append(mass)
    num_of_decomps_list.append(num_of_decomps)
    possible_formulas_list.append(possible_formulas)

os.chdir(original_working_directory) # go back to original working directory

mass_decomp_results = {
    'Mass': monoisotopic_masses_list,
    'RT': RT_list, # created earlier
    'RTmin': RTmin_list, 
    'RTmax': RTmax_list,
    'Number_of_decomps': num_of_decomps_list, 
    'Possible_formulas': possible_formulas_list
}

chemical_formulas_df = pd.DataFrame(mass_decomp_results)

print('All possible chemical formulas found!')

##########################################
#### STEP 5: Refined Hydrogen Rule ####
##########################################

chemical_formulas_df['nominal_mass'] = chemical_formulas_df['Mass'].apply(lambda mass: PredNomMass(mass, k, MSE, corr_fact, diff_mass_sum)[0]) #nominal mass predicted again because last time it was per RT spectrum

chemical_formulas_df['Possible_formulas_after_Hrule'] = chemical_formulas_df.apply(apply_rules, formula_column='Possible_formulas', axis=1)

chemical_formulas_df.to_csv('trackable_outputs/chemical_formulas_after_Hrule.csv', index=False)

###########################################
#### STEP 5: Subclass Prediction (KMD) ####
###########################################

lipidmaps_df_path = 'trackable_outputs/chemical_formulas_after_Hrule.csv' # no need to reload df but was small issue will come back to
lipidmaps_df = pd.read_csv(lipidmaps_df_path) 
mass_decomp_formula_column = 'Possible_formulas_after_Hrule'

reference_database_path = 'database_folder/KMD_reference_database.csv'
reference_entries = load_reference_database(reference_database_path)

dictionaries = []
for idx, row in tqdm(lipidmaps_df.iterrows(), total=lipidmaps_df.shape[0], desc="Finding possible subclasses..."):
    result = find_kmd_matches(row[mass_decomp_formula_column], idx, reference_entries)
    dictionaries.append(result)

lipidmaps_df['Formula:subclass combos'] = dictionaries
lipidmaps_df.to_csv(f'{raw_file_path}_output.csv.zip', compression='zip', index=False)


#############
#### FIN ####
#############

print('Done!')

end_time = time.time()
print('Total Time elapsed:', end_time - start_time, 'seconds')
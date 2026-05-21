import pickle
import argparse
import pandas as pd
import numpy as np
from rdkit import DataStructs

parser = argparse.ArgumentParser()

parser.add_argument('--prolif_tsv', '-i')
parser.add_argument('--outfile', '-o', help='Name of output .pkl file with sparse bit fingerprints for each target')

args = parser.parse_args()

def filter_rec_resi_l(df):
    # Filter the receptor residue list to only use
    # nonredundant residue names. Assumes numbering is the same

    rec_resi_l = []

    for resi_ch in df['rec_resi']:
        resi = resi_ch.split('.')[0]
        if resi not in rec_resi_l:
            rec_resi_l.append(resi)
    
    
    return rec_resi_l

def assign_fp_to_target(df, target, mapping, size=1024):
    fp_vector = np.zeros(size)

    df_filter = df[df['target'] == target]
    #print(target)
    #print(df_filter)

    for i, resi_ch in enumerate(df_filter['rec_resi']):
        resi = resi_ch.split('.')[0]
        plf_type = df_filter['ProLIF_type'].iloc[i]
        idx = mapping[resi][plf_type]
        #print(resi_ch, resi, plf_type, idx)

        fp_vector[idx] = 1

    return fp_vector

def assign_sparse_fp_to_target(df, target, mapping, size=1024):
    fp_vector = np.zeros(size)

    df_filter = df[df['target'] == target]
    #print(target)
    #print(df_filter)

    for i, resi_ch in enumerate(df_filter['rec_resi']):
        resi = resi_ch.split('.')[0]
        plf_type = df_filter['ProLIF_type'].iloc[i]
        idx = mapping[resi][plf_type]
        #print(resi_ch, resi, plf_type, idx)

        fp_vector[idx] = 1
    
    fp_sparse = np.flatnonzero(fp_vector)

    return fp_sparse

def create_resi_fp_mapping(resi_l, prolif_l):
    mapping = {}
    n_plf = len(prolif_l)
    for i, resi in enumerate(resi_l):
        mapping[resi] = {}
        for j, plf in enumerate(prolif_l):
            idx = i*n_plf + j
            mapping[resi][plf] = idx

            print(resi, plf, idx)
            
    
    return mapping

def tanimoto_numpy(array1, array2):
    # Intersection: count of common '1' bits
    intersection = np.dot(array1, array2.T)

    # Union components: total '1' bits in each array
    norm1 = np.sum(array1)
    norm2 = np.sum(array2)
    
    # Tanimoto formula: intersection / (sum_a + sum_b - intersection)
    return intersection / (norm1 + norm2 - intersection)

def main():
        
    # target  rec_resi  ProLIF_type
    df = pd.read_csv(args.prolif_tsv, delimiter='\t')
    print(df)

    rec_resi_l = filter_rec_resi_l(df)

    # Use only prolif types in the dataframe
    prolif_types = list(set(df['ProLIF_type'])) 

    target_l = list(set(df['target']))

    # Use default prolif types
    prolif_types = ['Hydrophobic', 'HBDonor', 'HBAcceptor', 'PiStacking', 'Anionic', 'Cationic', 'CationPi', 'PiCation', 'VdWContact']

    print(prolif_types)
    print(len(rec_resi_l), len(set(rec_resi_l)))

    fp_len = len(rec_resi_l) * len(prolif_types)

    bitsize = 1024 * int(fp_len%1024)
    
    # Create an index mapping from receptors+interactiosn to the bit string
    plf_fp_mapping = create_resi_fp_mapping(rec_resi_l, prolif_types)
    
    # Convert fingerprints to numpy arrays
    fp_l = []
    fp_data = {'fp_len': bitsize, 'targets': {} }
    for target in target_l:
        #fp_vec = assign_fp_to_target(df, target, plf_fp_mapping, size=bitsize)
        fp_vec = assign_sparse_fp_to_target(df, target, plf_fp_mapping, size=bitsize)
        fp_l.append(fp_vec)

        fp_data['targets'][target] = fp_vec

    
    with open(args.outfile, 'wb') as f:
        pickle.dump(fp_data, f)

if __name__=='__main__':
    main()

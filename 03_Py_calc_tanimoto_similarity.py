import json
import tqdm
import pickle
import argparse
import numpy as np
import pandas as pd
from rdkit.ML.Cluster import Butina

parser = argparse.ArgumentParser()

parser.add_argument('--fp_pkl', '-i', help='pkl file with dict for targets+fingerprints')
parser.add_argument('--clust_threshold', '-c', help='Distance threshold to use for clustering by ProLIF tanimoto similarity (default = 0.5)', default=0.5, type=float)

parser.add_argument('--filter_f', '-f', help='(Optional) provide a filter file to only consider a subset of targets. Should contain a "complex_name" column, and a "fragment_screen" column to denote which targets are fragments. If this file is provided, fragment-lead similarity will be calculated.', default=None)
parser.add_argument('--outfile', '-o', help='(Optional) Name of output file for frag-lead tanimoto (default = json_frag-lead_tanimoto.json)', default='json_frag-lead_tanimoto.json')

args = parser.parse_args()

# Calculate tanimoto similarity between sparse fingerprints
def calc_tanimoto(fp1, fp2):
    
    norm1 = len(fp1)
    norm2 = len(fp2)

    intersect = np.intersect1d(fp1, fp2)
    l_intersect = len(intersect)

    #print('fp1', fp1)
    #print('fp2', fp2)
    #print('intersect', intersect)

    tanimoto = l_intersect / (norm1 + norm2 - l_intersect)
    return tanimoto

def filter_targets(filter_f):
    df = pd.read_csv(filter_f)
    
    target_l_f = []
    frag_l = []
    lead_l = []
    for i, target in enumerate(df['complex_name']):
        target_l_f.append(target)
        is_frag = df['fragment_screen'].iloc[i]

        if is_frag:
            frag_l.append(target)
        else:
            lead_l.append(target)

    return target_l_f, frag_l, lead_l

def make_merged_frag_fp(frag_l, fp_data):
    merged_fp = []
    for frag in frag_l:
        for idx in fp_data['targets'][frag]:
            if int(idx) not in merged_fp:
                merged_fp.append(int(idx))

    merged_fp.sort()
    
    return np.array(merged_fp)

# Count how many "new" interactions a lead fp has when compared
# to all known fragment interactions
def count_new_interactions(lead_fp, merged_frag_fp):
    cnt = 0
    intersect = np.intersect1d(lead_fp, merged_frag_fp)

    new_interactions = len(lead_fp) - len(intersect)
    print(lead_fp, intersect)
    print(new_interactions)
    
    return new_interactions

def lead_to_frag_fp(frag_l, lead_l, fp_data, merged_frag_fp):
    out_data = {}

    for lead in lead_l:
        try:
            fp1 = fp_data['targets'][lead]
        except:
            continue
        

        n_new_int = count_new_interactions(fp1, merged_frag_fp)
        merged_tanimoto = calc_tanimoto(fp1, merged_frag_fp)
        max_tanimoto = 0
        max_frag = None
        for frag in frag_l:
            fp2 = fp_data['targets'][frag]

            tanimoto = calc_tanimoto(fp1, fp2)

            if tanimoto > max_tanimoto:
                max_tanimoto = tanimoto
                max_frag = frag

        print(f'{lead} -> {max_frag}: {max_tanimoto}')
        print(f'{lead} -> merged frags: {merged_tanimoto}')
        
        out_data[lead] = {'max_tanimoto': max_tanimoto, 
                          'n_new_interactions': n_new_int}

    return out_data

def main():
    with open(args.fp_pkl, 'rb') as f:
        fp_data = pickle.load(f)


    target_l = list(fp_data['targets'].keys())

    # Filter the targets if a file is provided
    if args.filter_f != None:
        target_l, frag_l, lead_l = filter_targets(args.filter_f)
        n_target = len(target_l)
        merge_frag_fp = make_merged_frag_fp(frag_l, fp_data)
        tanimoto_output = lead_to_frag_fp(frag_l, lead_l, fp_data, merge_frag_fp)

        with open('json_frag-lead_tanimoto.json', 'w') as fo:
            json.dump(tanimoto_output, fo, indent=4)

    # Calculate pairwise tanimoto similarity
    n_target = len(target_l)
    distance_mtx = np.zeros((n_target, n_target))

    for i, target in enumerate(tqdm.tqdm(target_l)):
        fp1 = fp_data['targets'][target]
        for j in range(0, i+1):
            target2 = target_l[j]
            fp2 = fp_data['targets'][target2]

            tanimoto = calc_tanimoto(fp1, fp2)
            distance = (1 - tanimoto)

            #print(i, target, j, target2, tanimoto)

            distance_mtx[i][j] = distance
            distance_mtx[j][i] = distance

    clusters = Butina.ClusterData(data=distance_mtx, nPts=n_target, distThresh=args.clust_threshold, isDistData=True)
    
    # Sort by size
    clusters = list(clusters)
    clusters.sort(key=len, reverse=True)

    print(f'{len(clusters)} clusters identified.')
    # Save target names for each cluster
    out_data = {}
    for i, clust in enumerate(clusters):
        clust_n = f'clust.{i}.{len(clust)}'
        clust_l = []
        for idx in clust:
            clust_l.append(target_l[idx])

        out_data[clust_n] = clust_l

    with open(args.outfile, 'w') as fo:
        json.dump(out_data, fo, indent=4)   
            
        
    

if __name__=='__main__':
    main()

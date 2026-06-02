import os
import json
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('--clust_json', '-cj', help='json with prolif clusters')
parser.add_argument('--fragalysis_dir', '-fd', help='Path to aligned_files/ directory for fragalysis hits')
parser.add_argument('--outdir', '-od', help='Path to output directory')

args = parser.parse_args()

def main():
    os.makedirs(args.outdir, exist_ok=True)

    with open(args.clust_json) as f:
        cj_data = json.load(f)

    for clust in cj_data:
        cj_lines = []
        for target in cj_data[clust]:
            t_lig = f'{args.fragalysis_dir}/{target}/{target}_ligand.sdf'

            with open(t_lig) as sdf:
                sdf_lines = sdf.readlines()

            cj_lines += sdf_lines

        with open(f'{args.outdir}/{clust}_ligands.sdf', 'w') as fo:
            fo.write(''.join(cj_lines))
        




if __name__=='__main__':
    main()

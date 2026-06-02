import os
import sys
import argparse
import prolif as plf
import pandas as pd
import numpy as np
from pymol import cmd
from rdkit import Chem

parser = argparse.ArgumentParser()

parser.add_argument('--fragalysis_dir', '-fd', help='aligned_files/ directory for a fragalysis target')
parser.add_argument('--outfile', '-o', help='Name of output .tsv file')

args = parser.parse_args()

def main():
    outlines = ['target\trec_resi\tProLIF_type']
    err_log = []
    for target in os.listdir(args.fragalysis_dir):
        print(target)
        t_path = f'{args.fragalysis_dir}/{target}/'

        if os.path.isdir(t_path) == False:
            continue

        rec_pdb = f'{t_path}/{target}_delig-desolv.pdb'
        lig_sdf = f'{t_path}/{target}_ligand.sdf'

        cmd.reinitialize()
        cmd.load(rec_pdb)
        cmd.h_add()
        cmd.save('tmp.pdb')
        
        try:
            rdkit_prot = Chem.MolFromPDBFile('tmp.pdb', removeHs=False)
            protein_mol = plf.Molecule(rdkit_prot)

            rdkit_lig = Chem.MolFromMolFile(lig_sdf)
            ligand_mol = plf.Molecule(rdkit_lig)
        except Exception as e:
            err_log.append(f'RDKIT_ERR: {target}\n')
            err_log.append(str(e) + '\n')
            continue

        print('\t', protein_mol)
        print('\t', ligand_mol)

        try:
            fp = plf.Fingerprint(count=True)
            ifp = fp.generate(ligand_mol, protein_mol, metadata=True)
            ifp_df = plf.to_dataframe({0: ifp}, fp.interactions, dtype=np.uint8)
            print(ifp_df)
        except Exception as e:
            err_log.append(f'PROLIF_ERR: {target}\n')
            err_log.append(str(e) + '\n')
            continue
            
        # Save output lines
        for h in ifp_df.head():
            print('\t', h)
            outline = f'{target}\t{h[1]}\t{h[2]}'
            outlines.append(outline)

    with open(args.outfile, 'w') as fo:
        fo.write('\n'.join(outlines))

    with open(f'{args.outfile[:-4]}.err', 'w') as fo:
        for l in err_log:
            fo.write(l)

if __name__=='__main__':
    main()

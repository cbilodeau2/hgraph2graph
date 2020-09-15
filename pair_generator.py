import numpy as np
import sys
import os
from chemprop.data.data import MoleculeDatapoint
from chemprop.data.data import MoleculeDataset
from chemprop.data.scaffold import scaffold_split
from chemprop.data.scaffold import *
import argparse
import pandas as pd
import time
from rdkit import Chem


from hgraph import common_atom_vocab

def get_index_smiles(index_set,data_set):
    return [data_set[i].smiles for i in index_set]
            
def drop_dots(df):
    # Removes any molecular entries contianing more than one molecule
    
    return df[df['SMILES'].apply(lambda smiles: (bool(not '.' in smiles)))]

def drop_unusual_atoms(df):
    # Removes any molecular entries containing atoms not present in the common atom vocab
    
    def check_atoms(smiles):
        for atom in Chem.MolFromSmiles(smiles).GetAtoms():
            label = (atom.GetSymbol(), atom.GetFormalCharge())
            valid = True
            if not (label in common_atom_vocab.vocab):
                valid = False
                break
        return valid
    
    return df[df['SMILES'].apply(check_atoms)]

def smile2target(smiles,data,target):
    # Takes in smiles list and dataframe containing target data and gives list of targets
    
    targets = []
    for x in smiles:
        targets.append(data[data['SMILES']==x][target].values[0])
    return(targets)

def smile2pairs(smiles, data, target='Solubility', cutoff= 0.78*2, sample_n=20):
    # smiles: list of smiles to be sorted into pairs
    # data: dataframe containing smiles/target data
    # cutoff: minimum target difference between x/y
    # number of times pairs for each x molecule to try to participate in
    
    # Assign all molecules below median target to x category:
    df = pd.DataFrame([{"Smile":x,"Target":y} for x,y in zip(smiles, smile2target(smiles,data,target))])
    med = df['Target'].median()
    df_x_sorted = df[df['Target']<med].sort_values('Target',ascending=False)
    df_y = df[df['Target']>med]
    
    # Assign y molecules for each x molecule:
    assigned_pairs = pd.DataFrame()
    for i in range(0,len(df_x_sorted)):
        smile1 = df_x_sorted.iloc[i]['Smile']
        targetx = df_x_sorted.iloc[i]['Target']

        df_y_subset = df_y[(df_y['Target']-df_x_sorted.iloc[i]['Target'])>cutoff]

        # Place a limit on permutations:
        if len(df_y_subset)>sample_n:
            df_y_subset = df_y_subset.sample(n = sample_n)

        for j in range(0,len(df_y_subset)):
            smile2 = df_y_subset.iloc[j]['Smile']
            targety = df_y_subset.iloc[j]['Target']

            assigned_pairs = assigned_pairs.append({'X':smile1,'Y':smile2,'targetx':targetx,'targety':targety},ignore_index=True)
            if len(assigned_pairs)%500==0:
                print('{} Pairs Assigned'.format(len(assigned_pairs)))
    
    return assigned_pairs

def remove_tails(df,target):
    df_sorted = df.sort_values(target)
    df_notails = df_sorted.iloc[int(len(df_sorted)/8):int(len(df_sorted)*7/8)]
    
    # Return without tails and shuffled:
    return df_notails.sample(frac=1)

def drop_invalid_smiles(df):
    for smile in df['SMILES']:
        try:
            mol = Chem.MolFromSmiles(smile)
        except:
            df = df[df['SMILES']!=smile]
    return df

def generate_pairs(data_file,outfile,molfile,target='Solubility',cutoff= 0.78*2, sample_n=20,remove_tails_flag=True):  
    data = pd.read_csv(data_file)

    # Apply data filters:
    data = drop_dots(data)
    data = drop_unusual_atoms(data)
    data= drop_invalid_smiles(data)
    data= data[['SMILES',target]]
    
    if remove_tails_flag:
        data = remove_tails(data,target)

    # Generate Dataset object:
    data_list = []
    for i in range(0,len(data)):
        data_list.append(MoleculeDatapoint(smiles = data.iloc[i].values[0], targets = data.iloc[i].values[1]))
    data_set = MoleculeDataset(data_list)

    # Generate Mapping:
    scaffold_to_indices = scaffold_to_smiles(data_set.mols(), use_indices=True)
    index_sets = list(scaffold_to_indices.values())

    # Generate Pairlist:
    g2g = pd.DataFrame(columns=['X','Y'])
    i=0
    print('{} Total Scaffolds'.format(len(index_sets)))
    for index_set in index_sets:
        i+=1
        print('Starting Scaffold {}'.format(i))
        smile_list = get_index_smiles(index_set,data_set)    
        # There must be at least two molecules in the set for a pair to be made:
        if len(smile_list)>1:
            pairs = smile2pairs(smile_list,data,target,cutoff,sample_n)
            g2g = pd.concat((g2g,pairs),axis=0)
        else:
            print('No pair assignment')
            
    print('Number of Pairs:{}'.format(len(g2g)))
    g2g['X'] = [x.strip("\r\n ") for x in g2g['X'].values]
    g2g['Y'] = [x.strip("\r\n ") for x in g2g['Y'].values]
    
    # Create output files:
    mols = pd.DataFrame(data=list(set(list(g2g['X'].values)+list(g2g['Y'].values))))
    
    g2g[['X','Y']].to_csv(outfile,index=False,header=None,sep=' ')
    mols.to_csv(molfile,index=False,header=None,sep=' ')
            
if __name__ == "__main__":
    
    start = time.time()
    
    parser = argparse.ArgumentParser()
    parser.add_argument('--infile', required=True)
    parser.add_argument('--outfile', required=True)
    parser.add_argument('--molfile', required=True)
    parser.add_argument('--target', type=str,default='Solubility')
    parser.add_argument('--cutoff', type=float, default=0.78*2)
    parser.add_argument('--sample', type=int, default=20)
    parser.add_argument('--remove_tails', type=bool, default=True)
    
    args = parser.parse_args()
    
    generate_pairs(args.infile,
                   args.outfile,
                   args.molfile,
                   args.target,
                   args.cutoff,
                   args.sample,
                   args.remove_tails)
    
    end = time.time()
    
    print('Completed. Time Elapsed:{}'.format(end-start))
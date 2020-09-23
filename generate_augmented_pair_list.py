import numpy as np
import pandas as pd
import argparse

from rdkit.Chem.Descriptors import ExactMolWt
from rdkit import Chem

from chemprop.data.scaffold import generate_scaffold
from hgraph.pairing import *

# SETTINGS:
parser = argparse.ArgumentParser()

parser.add_argument('--gen_out', type=str, required=True) # Input: Pairs generated during generation step
parser.add_argument('--gen_evaluated', type=str, required=True) # Input: Predictions generated for all molecules after generation step
parser.add_argument('--datafile', type=str, required=True) # Output: Labeled data file
parser.add_argument('--molfile', type=str, required=True) # Output: Mol file for next iteration
parser.add_argument('--pairfile', type=str, required=True) # Output: Pair file for next iteration
parser.add_argument('--threshold', type=float, default=0.78) # Molecules must improve by this much to be considered a successful translation
parser.add_argument('--pairgen_cutoff', type=float, default=0.78) # Threshold used to generate new pairs
parser.add_argument('--pairgen_sample', type=int, default=20) # Number of pairs to generate per molecule
parser.add_argument('--min_mol_wt', type=float, default=50.0) # Lowest molecular weight to include for successful translation
parser.add_argument('--target', type=str, default='Target') # Y column title (in gen_evaluated input file)

args = parser.parse_args()

paired = pd.read_csv(args.gen_out,header=None,sep=' ')
labeled = pd.read_csv(args.gen_evaluated)
paired = paired.rename(columns={0:'X',1:'Y'})

# Remove molecules that don't improve enough:
paired['Targetx']=paired['X'].apply(lambda x: labeled[labeled[labeled.columns[0]]==x][labeled.columns[1]].iloc[0])
paired['Targety']=paired['Y'].apply(lambda x: labeled[labeled[labeled.columns[0]]==x][labeled.columns[1]].iloc[0])
paired = paired[(paired['Targety']-paired['Targetx'])>args.threshold]

# Remove Y molecules with low mw:
paired['MolWt Y'] = paired['Y'].apply(lambda x: ExactMolWt(Chem.MolFromSmiles(x)))
paired = paired[paired['MolWt Y']>args.min_mol_wt]

# Remove molecules outside scaffold
paired['Scaffoldx'] = paired['X'].apply(generate_scaffold)
paired['Scaffoldy'] = paired['Y'].apply(generate_scaffold)
paired = paired[paired['Scaffoldx']==paired['Scaffoldy']]

# Make labeled dataset for input into next iteration:
x_labeled = paired[['X','Targetx']].rename(columns={'X':'SMILES','Targetx':'Target'})
y_labeled = paired[['Y','Targety']].rename(columns={'Y':'SMILES','Targety':'Target'})
labeled = pd.DataFrame(labeled.dropna().values,columns=['SMILES','Target'])

data_out = pd.concat([x_labeled,y_labeled,labeled]).drop_duplicates()

# Save data.csv file
data_out.to_csv(args.datafile,index=False)

# Generate new input files for next iteration (mols.txt, train_pairs.txt)
generate_pairs(args.datafile,args.pairfile,args.molfile,'Target',args.pairgen_cutoff, args.pairgen_sample,remove_tails_flag=False)



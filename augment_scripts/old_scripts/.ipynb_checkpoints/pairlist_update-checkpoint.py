import numpy as np
import pandas as pd
import argparse
import os

# Takes in list from decode.py
# Generates data required for next round of translation

parser = argparse.ArgumentParser()
parser.add_argument('--input', type=str, required=True) # New pairs from decode.py
parser.add_argument('--output', type=str, required=True) # New train_pairs.txt for the next round
parser.add_argument('--mol_list', type=str, required=True) # New molecule list for vocab.py
parser.add_argument('--folder', type=str, required=True) # Folder where extra files generated should be stored
parser.add_argument('--old_pairs', type=str, required=True) # Previous pairlist to append new pairlist to
parser.add_argument('--tox', type=bool, 
                   default=False) # Include toxicity constraint
parser.add_argument('--predictor', type=str, # Chemprop model to use for prediction
                   default='/data/rsg/chemistry/cbilod/chemprop_old/checkpoints/solubility/solubility-1594231605.496638/fold_0')
parser.add_argument('--tox_predictor', type=str, # Chemprop model to use for prediction
                   default='/data/rsg/chemistry/cbilod/chemprop_old/checkpoints/tox21/tox21-1596485275.340425/fold_0')
parser.add_argument('--threshold',type=float, default=0.79) # Minimum improvement to include
args = parser.parse_args()

results = pd.read_csv(args.input,sep=' ',header=None)

# Write out columns for prediction:
results[0].to_csv(os.path.join(args.folder,'col1.csv'),index=False)
results[1].to_csv(os.path.join(args.folder,'col2.csv'),index=False)

# Use Chemprop to predict solubilities:
os.system('python /data/rsg/chemistry/cbilod/chemprop/predict.py --test_path '+os.path.join(args.folder,'col1.csv')+' --checkpoint_dir '+args.predictor+' --preds_path '+os.path.join(args.folder,'preds_col1.csv'))
os.system('python /data/rsg/chemistry/cbilod/chemprop/predict.py --test_path '+os.path.join(args.folder,'col2.csv')+' --checkpoint_dir '+args.predictor+' --preds_path '+os.path.join(args.folder,'preds_col2.csv'))

# Read in and organize predictions:
preds1 = pd.read_csv(os.path.join(args.folder,'preds_col1.csv'))
preds1 = preds1.rename(columns={"0":"Mol1","Solubility":"Sol1"})
preds2 = pd.read_csv(os.path.join(args.folder,'preds_col2.csv'))
preds2 = preds2.rename(columns={"1":"Mol2","Solubility":"Sol2"})
preds_tot = pd.concat((preds1,preds2),axis=1)

# New pairs are molecules that have been improved beyond their threshold:
new_pairs = preds_tot[preds_tot['Sol2']-preds_tot['Sol1']>args.threshold][['Mol1','Mol2']]

# Toxicity constraint
if args.tox == True:
    new_pairs = apply_constraints(new_pairs,args.tox_predictor,args.folder)

# Combine new pairs with old pairs
old_pairs = pd.read_csv(args.old_pairs,delimiter=' ',header=None)
old_pairs.columns = ['Mol1','Mol2']
all_pairs = pd.concat([old_pairs,new_pairs],axis=0)
print('Previous list contained {}'.format(len(old_pairs)))
print('Adding {} pairs'.format(len(new_pairs)))

# Write pairs to csv:
all_pairs.to_csv(args.output,sep=' ',header=None,index=None)
mol_list = pd.concat([all_pairs['Mol1'],all_pairs['Mol2']]).drop_duplicates()
mol_list.to_csv(args.mol_list,header=None,index=False)
print('{} unique molecules'.format(len(mol_list)))
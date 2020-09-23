import numpy as np
import pandas as pd
import argparse
import os

# Takes in list from decode.py
# Generates data required for next round of translation


def apply_constraint(old_pairs,predictor,folder,threshold=0.0168): # This is a terrible threshold that should be fixed
    
    # Make folder to store intermediate files
    folder = os.path.join(folder,'tox_constraint')
    try:
        os.mkdir(folder)
    except:
        print("File Exists: {}".format(folder))
    
    results = pd.read_csv(args.input,sep=' ',header=None)
    
    # Write out columns for prediction:
    results[0].to_csv(os.path.join(folder,'col1.csv'),index=False)
    results[1].to_csv(os.path.join(.folder,'col2.csv'),index=False)
    
    # Use Chemprop to predict solubilities:
    os.system('python /data/rsg/chemistry/cbilod/chemprop/predict.py --test_path '+os.path.join(folder,'col1.csv')+' --checkpoint_dir '+predictor+' --preds_path '+os.path.join(args.folder,'preds_col1.csv'))
    os.system('python /data/rsg/chemistry/cbilod/chemprop/predict.py --test_path '+os.path.join(folder,'col2.csv')+' --checkpoint_dir '+predictor+' --preds_path '+os.path.join(args.folder,'preds_col2.csv'))
    
    # Read in and organize predictions:
    preds1 = pd.read_csv(os.path.join(args.folder,'preds_col1.csv'))
    preds1 = preds1.rename(columns={"0":"Mol1",preds1.columns[1]:"Tox1"})
    preds2 = pd.read_csv(os.path.join(args.folder,'preds_col2.csv'))
    preds2 = preds2.rename(columns={"1":"Mol2",preds2.columns[1]:"Tox2"})
    preds_tot = pd.concat((preds1,preds2),axis=1)
    
    preds_tot = preds_tot[preds_tot['Tox2']<threshold] # Remove molecules that translate into toxic molecules
    
    return preds_tot[['Mol1','Mol2']]
    
            
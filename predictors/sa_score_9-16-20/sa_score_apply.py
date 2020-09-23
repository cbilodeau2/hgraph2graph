import argparse
import pandas as pd
import rdkit
from rdkit import Chem
from rdkit.Chem import RDConfig
import os
import sys
sys.path.append(os.path.join(RDConfig.RDContribDir, 'SA_Score')) 
import sascorer
import time

def sa_score_apply(infile,outfile):

    data = pd.read_csv(infile)
    data['SA Score']=data['SMILES'].apply(lambda x: sascorer.calculateScore(rdkit.Chem.MolFromSmiles(x)))
    data[['SMILES','SA Score']].to_csv(outfile,index=False)
    
if __name__ == "__main__":
    
    start = time.time()
    
    parser = argparse.ArgumentParser()
    parser.add_argument('--infile', required=True)
    parser.add_argument('--outfile', required=True)
    
    args = parser.parse_args()
    
    sa_score_apply(args.infile,args.outfile)
    
    
    
import numpy as np
import pandas as pd
import argparse


parser = argparse.ArgumentParser()
parser.add_argument('--input', type=str, required=True)
parser.add_argument('--output', type=str, required=True)
parser.add_argument('--selection', type=str, default='whole') # half
args = parser.parse_args()

data = pd.read_csv(args.input)
print(data.head())
data = data.drop_duplicates() # Remove duplicates

if args.selection=='whole':
    out_data = data.sort_values('Solubility',ascending=True)[['SMILES']]
elif args.selection=='half':
    out_data = data.sort_values('Solubility',ascending=True).iloc[int(len(data)/2):][['SMILES']]
out_data.to_csv(args.output,index=False,header=None)


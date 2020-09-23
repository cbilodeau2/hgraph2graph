import numpy as np
import pandas as pd
import argparse


def check_smiles_match(data,screen):
    return (data['SMILES'].values==screen['SMILES'].values).all()

def apply_screen(data,col_name,selection_type,selection_thresh,keep):
    data = data.sort_values(col_name,ascending=True)
    if selection_type=='Fraction':
        if keep=='High':
            data = data[-int(len(data)*selection_thresh):]
        elif keep=='Low':
            data = data[0:-int(len(data)*selection_thresh)]
        else:
            print('WARNING: INVALID KEEP TYPE')
    elif selection_type=='Cutoff':
        if keep=='High':
            data = data[data[col_name]>selection_thresh]
        elif keep=='Low':
            data = data[data[col_name]<selection_thresh]
        else:
            print('WARNING: INVALID KEEP TYPE')            
    else:
        print('WARNING: INVALID SELECTION TYPE')
        
    return data


parser = argparse.ArgumentParser()
parser.add_argument('--molfile', type=str, required=True)
parser.add_argument('--outfile', type=str, required=True)

parser.add_argument('--screen_file1', type=str, default=None)
parser.add_argument('--selection_type1', type=str, default='Fraction') # Fraction or Cutoff Value
parser.add_argument('--selection_thresh1', type=float, default=0.5)
parser.add_argument('--keep1', type=str, default='High') # High or low

parser.add_argument('--screen_file2', type=str, default=None)
parser.add_argument('--selection_type2', type=str, default='Cutoff') # Fraction or Cutoff Value
parser.add_argument('--selection_thresh2', type=float, default=5.0)
parser.add_argument('--keep2', type=str, default='Low') # High or low

args = parser.parse_args()


data = pd.read_csv(args.molfile)
data = data.drop_duplicates() # Remove duplicates

if args.screen_file1 is not None:
    screen1 = pd.read_csv(args.screen_file1)
    
    # Check if smiles match:
    if not check_smiles_match(data,screen1):
        print('WARNING: SMILES LISTS DO NOT MATCH')
        
    # Add screen
    col_name1 = pd.DataFrame(screen1.columns)[[not (x =='SMILES') for x in screen1.columns]].values[0][0]
    data[col_name1]=screen1[col_name1]
    
if args.screen_file2 is not None:
    screen2 = pd.read_csv(args.screen_file2)
    
    # Check if smiles match:
    if not check_smiles_match(data,screen1):
        print('WARNING: SMILES LISTS DO NOT MATCH')
    
    # Add screen
    col_name2 = pd.DataFrame(screen2.columns)[[not (x =='SMILES') for x in screen2.columns]].values[0][0]
    data[col_name2]=screen2[col_name2]


# Apply screens
if args.screen_file1 is not None:
    data = apply_screen(data,col_name1,args.selection_type1,args.selection_thresh1,args.keep1)
if args.screen_file2 is not None:
    data = apply_screen(data,col_name2,args.selection_type2,args.selection_thresh2,args.keep2)

    
# Output data for generation
out_data = data.sample(frac=1.0) #Shuffle data
out_data[['SMILES']].to_csv(args.outfile,index=False,header=None)
    

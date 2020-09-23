import numpy as np
import pandas as pd
import argparse


parser = argparse.ArgumentParser()

parser.add_argument('--infile', type=str, required=True)
parser.add_argument('--outfile', type=str, required=True)
parser.add_argument('--old_list', type=str, required=True)

args = parser.parse_args()


data = pd.read_csv(args.infile,header=None,sep=' ')
old_list = pd.read_csv(args.old_list)
data = list(np.reshape(pd.concat((data,old_list)).to_numpy(),(-1,)))
data = set(data)
data = pd.DataFrame(data)

data.to_csv(args.outfile, sep=' ', index=False)



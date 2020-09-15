#!/bin/bash

# Starting files:
#     data/$folder/data.csv

folder=solvation
target=Solubility

# Settings:
gpu=0
n_epochs=1
timestamp=0 # If presetting timestamp, comment timestamp below

export CUDA_VISIBLE_DEVICES=$gpu

# If generating new timestap:
timestamp=`date +"%s"`

echo "=========Create Pairlist========"
python pair_generator.py --infile data/$folder/data.csv --outfile data/$folder/train_pairs.txt --molfile data/$folder/mols.txt --target Solubility

echo "=========Get Vocabulary========"
python get_vocab.py <data/$folder/mols.txt > data/$folder/vocab.txt

echo "=========Preprocess Data========"
python preprocess.py --train data/$folder/train_pairs.txt --vocab data/$folder/vocab.txt --ncpu 16 < data/$folder/train_pairs.txt

mkdir checkpoints/$folder
dir=checkpoints/$folder/${folder}_${timestamp}
mkdir $dir
mkdir $dir/train_processed
mv tensor* $dir/train_processed/

echo "=========Train Model========"
mkdir $dir/models/
python gnn_train.py --train $dir/train_processed/ --vocab data/$folder/vocab.txt --save_dir $dir/models/ --epoch $n_epochs


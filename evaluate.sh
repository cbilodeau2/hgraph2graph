#!/bin/bash

checkpoint_path=checkpoints/solvation/solvation_1600177098 # This needs to be updated
data_path=data/solvation
model=model.0

python generate_tests.py --infile $data_path/data.csv --val_path $data_path

for section in bottom med_bottom medium med_top top
do
    mkdir $checkpoint_path/results_$section
    python decode.py --test $data_path/$section.txt --vocab $data_path/vocab.txt --model $checkpoint_path/models/$model --num_decode 20 > $checkpoint_path/results_$section/results.csv
done

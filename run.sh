#!/bin/bash

# source /data/rsg/chemistry/cbilod/anaconda/bin/activate
# conda activate chemprop

folder=sol_pairs_extended_20_notails_random_cutoff_0_78 #_cutoff_0_78

# echo "=========Get Vocabulary========"
# python get_vocab.py <data/$folder/mols.txt > data/$folder/vocab.txt

# echo "=========Preprocess Data========"
# python preprocess.py --train data/$folder/train_pairs.txt --vocab data/$folder/vocab.txt --ncpu 16 < data/$folder/train_pairs.txt
# mkdir checkpoints/$folder
# mkdir checkpoints/$folder/train_processed
# mv tensor* checkpoints/$folder/train_processed/

# echo "=========Train Model========"
# mkdir checkpoints/$folder/models/
# python gnn_train.py --train checkpoints/$folder/train_processed/ --vocab data/$folder/vocab.txt --save_dir checkpoints/$folder/models/

echo "=========Decode Model========"
#section=top
for section in bottom med_bottom medium med_top top
do
    mkdir checkpoints/$folder/results_$section
    python decode.py --test data/$folder/$section.txt --vocab data/$folder/vocab.txt --model checkpoints/$folder/models/model.5 --num_decode 20 > checkpoints/$folder/results_$section/results.csv
done
# python get_vocab.py < data/qed/mols.txt > vocab.txt

# python preprocess.py --train data/qed/train_pairs.txt --vocab data/qed/vocab.txt --ncpu 16 < data/qed/train_pairs.txt
# mkdir train_processed
# mv tensor* train_processed/


# mkdir models/
# python gnn_train.py --train train_processed/ --vocab data/qed/vocab.txt --save_dir models/ 
#mkdir checkpoints/$folder/results_top
#python decode.py --test data/$folder/top.txt --vocab data/$folder/vocab.txt --model checkpoints/$folder/models/model.2 --num_decode 20 > checkpoints/$folder/results_top/results.csv
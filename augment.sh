#!/bin/bash

# Note: Duplicates should be dropped in mols.txt to speed up process
# Must start by creating folders in checkpoints and in data folder
export CUDA_VISIBLE_DEVICES=3

name=sol_augment_nomask
predictor_path=/data/rsg/chemistry/cbilod/chemprop_old/checkpoints/solubility/solubility-1594231605.496638/fold_0
n_epochs=10

# If generating new timestap:
timestamp=`date +"%s"`
folder=$name/${name}-${timestamp}

# If restarting from existing checkpoint:
folder=$name/sol_augment_nomask-1597155462 #sol_augment-1596468619

# # Round 1 =================================================================

mkdir checkpoints/$folder
dir=checkpoints/$folder/generation1
mkdir $dir

# Extract Vocabulary:
echo "=========Get Vocabulary========"
python get_vocab.py <data/$name/mols.txt > data/$name/vocab.txt

# Preprocess Data:
echo "=========Preprocess Data========"
python preprocess.py --train data/$name/train_pairs.txt --vocab data/$name/vocab.txt --ncpu 16 < data/$name/train_pairs.txt
mkdir $dir/train_processed
mv tensor* $dir/train_processed/

# Train Model:
echo "=========Train Model========"
mkdir $dir/models/
python gnn_train.py --train $dir/train_processed/ --vocab data/$name/vocab.txt --save_dir $dir/models/ --epoch $n_epochs

# Predict solubility for input dataset:
mkdir $dir/augment/
echo "=========Predict Solubility========"
rm -rf data/$name/all_mols.txt 
echo "SMILES" >> data/$name/all_mols.txt
cat data/$name/mols.txt >> data/$name/all_mols.txt

cp data/$name/train_pairs.txt $dir/augment/

python /data/rsg/chemistry/cbilod/chemprop/predict.py --test_path data/$name/all_mols.txt --checkpoint_dir $predictor_path --preds_path $dir/augment/mols_pred.csv

# Select molecules in the top half of the dataset:
echo "=========Select Molecules for Generation========"
python augment_scripts/select_mols.py --input $dir/augment/mols_pred.csv --output $dir/augment/gen_in.csv

# Generate Molecules for the top half of the dataset:
echo "=========Decode Model========"
python decode.py --test $dir/augment/gen_in.csv --vocab data/$name/vocab.txt --model $dir/models/model.9 --num_decode 20 > $dir/augment/gen_out.csv

# Prune new pairs and add to pairlist:
echo "=========Add to Pairlist========"
python augment_scripts/pairlist_update.py --input $dir/augment/gen_out.csv --output $dir/augment/pairlist.txt --folder $dir/augment --old_pairs data/$name/train_pairs.txt --mol_list $dir/augment/mols.txt

# Additional Rounds =================================================================

for i in 4 5 #2 3 4 5
do
    old_dir=checkpoints/$folder/generation$(($i-1))
    dir=checkpoints/$folder/generation$i
    mkdir $dir
    
    echo "Working Directory is:" $dir
    
    # Extract Vocabulary:
    echo "=========Update Vocabulary========"
    python get_vocab.py <$old_dir/augment/mols.txt > $old_dir/augment/vocab.txt

    # Preprocess Data Round:
    echo "=========Preprocess Data 2========"
    python preprocess.py --train $old_dir/augment/pairlist.txt --vocab $old_dir/augment/vocab.txt --ncpu 16 < $old_dir/augment/pairlist.txt
    mkdir $dir/train_processed
    mv tensor* $dir/train_processed/

    # # Train Model Round:
    echo "=========Train Model 2========"
    mkdir $dir/models/
    python gnn_train.py --train $dir/train_processed/ --vocab $old_dir/augment/vocab.txt --save_dir $dir/models/ --epoch $n_epochs


    # Predict solubility for dataset:
    mkdir $dir/augment/
    echo "=========Predict Solubility========"
    rm -rf $old_dir/augment/all_mols.txt 
    echo "SMILES" >> $old_dir/augment/all_mols.txt
    cat $old_dir/augment/mols.txt >> $old_dir/augment/all_mols.txt

    python /data/rsg/chemistry/cbilod/chemprop/predict.py --test_path $old_dir/augment/all_mols.txt --checkpoint_dir $predictor_path --preds_path $dir/augment/mols_pred.csv

    # Select molecules in the top half of the dataset:
    echo "=========Select Molecules for Generation========"
    python augment_scripts/select_mols.py --input $dir/augment/mols_pred.csv --output $dir/augment/gen_in.csv

    # Generate Molecules Round 1 from molecules in the top half of the dataset:
    echo "=========Decode Model========"
    python decode.py --test $dir/augment/gen_in.csv --vocab $old_dir/augment/vocab.txt --model $dir/models/model.9 --num_decode 20 > $dir/augment/gen_out.csv

    # Prune pairs and add to pairlist:
    echo "=========Add to Pairlist========"
    python augment_scripts/pairlist_update.py --input $dir/augment/gen_out.csv --output $dir/augment/pairlist.txt --folder $dir/augment --old_pairs $old_dir/augment/pairlist.txt --mol_list $dir/augment/mols.txt

    echo "Working Directory is:" $dir
#input: gen_out.csv # generated molecules
#output: new_pairlist.txt # Pairlist for next iteration
done


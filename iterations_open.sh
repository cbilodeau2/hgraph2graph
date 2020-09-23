#!/bin/bash

# Starting files:
#     data/$folder/data.csv

folder=solvation_open
target=Solubility
log_folder=log_files
# Predictor 1 (this is assumed to be the main predictor for determining if translation is successful)
predictor1=predictors/solubility_9-21-20
predictor1_path=`cat predictors/solubility_9-21-20/path.txt`
predictor1_exec=predict.py

sa_score_path=predictors/sa_score_9-16-20
sa_score_exec=sa_score_apply.py

# Settings=====================================================================================
gpu=2
timestamp=1600708916 # If presetting timestamp, comment timestamp below
n_epochs=10
fin_model=model.9 # Max is n_epochs-1, model will be used for generation and starting checkpoint
num_decode=20 # Number of molecules to generate (per molecule) at each generation step

export CUDA_VISIBLE_DEVICES=$gpu

# If generating new timestap===================================================================
timestamp=`date +"%s"`
log_file=$log_folder/${folder}_${timestamp}.log

# ITERATION 1 =================================================================================
echo "========Starting Iteration 1========" >> $log_file

start=$(date +'%s')
echo "=========Create Pairlist========"
echo "Create Pairlist" >> $log_file
python pair_generator.py --infile data/$folder/data.csv --outfile data/$folder/train_pairs.txt --molfile data/$folder/mols.txt --target Solubility
end=$(date +'%s')
echo Execution time was $(($end-$start)) seconds. >> $log_file


start=$(date +'%s')
echo "=========Get Vocabulary========"
echo "Get Vocabulary" >> $log_file
python get_vocab.py <data/$folder/mols.txt > data/$folder/vocab.txt
end=$(date +'%s')
echo Execution time was $(($end-$start)) seconds. >> $log_file


start=$(date +'%s')
echo "=========Preprocess Data========"
echo "Preprocess Data" >> $log_file
python preprocess.py --train data/$folder/train_pairs.txt --vocab data/$folder/vocab.txt --ncpu 16 < data/$folder/train_pairs.txt
end=$(date +'%s')
echo Execution time was $(($end-$start)) seconds. >> $log_file


mkdir checkpoints/$folder
mkdir checkpoints/$folder/${folder}_${timestamp}
dir=checkpoints/$folder/${folder}_${timestamp}/iteration1
mkdir $dir
mkdir $dir/train_processed
mv tensor* $dir/train_processed/


start=$(date +'%s')
echo "=========Train Model========"
echo "Train Model" >> $log_file
mkdir $dir/models/
python gnn_train.py --train $dir/train_processed/ --vocab data/$folder/vocab.txt --save_dir $dir/models/ --epoch $n_epochs
end=$(date +'%s')
echo Execution time was $(($end-$start)) seconds. >> $log_file


mkdir $dir/augment
rm -rf $dir/augment/mols_old.txt
echo "SMILES" >> $dir/augment/mols_old.txt
cat data/$folder/mols.txt >> $dir/augment/mols_old.txt

# start=$(date +'%s')
# echo "=========Make Predictions: Model 1========"
# echo "Make Predictions: Model 1" >> $log_file
# python $predictor1_path/$predictor1_exec --test_path $dir/augment/mols_old.txt --checkpoint_dir $predictor1 --preds_path $dir/augment/mols_pred_model1.csv
# end=$(date +'%s')
# echo Execution time was $(($end-$start)) seconds. >> $log_file

# start=$(date +'%s')
# echo "=========Make Predictions: SA Score========"
# echo "Make Predictions: SA Score" >> $log_file
# python $sa_score_path/$sa_score_exec --infile $dir/augment/mols_old.txt --outfile $dir/augment/mols_pred_sa_score.csv
# end=$(date +'%s')
# echo Execution time was $(($end-$start)) seconds. >> $log_file

#### HERE IS WHERE WE CHANGE MASKING PROCEDURE:

start=$(date +'%s')
echo "=========Select Molecules for Generation========"
echo "Select Molecules for Generation" >> $log_file
python augment_scripts/select_mols.py --molfile $dir/augment/mols_old.txt --outfile $dir/augment/gen_in.csv #--screen_file1 $dir/augment/mols_pred_model1.csv #--screen_file2 $dir/augment/mols_pred_sa_score.csv
end=$(date +'%s')
echo Execution time was $(($end-$start)) seconds. >> $log_file


start=$(date +'%s')
echo "=========Generate Molecules========"
echo "Generate Molecules" >> $log_file
python decode.py --test $dir/augment/gen_in.csv --vocab data/$folder/vocab.txt --model $dir/models/$fin_model --num_decode $num_decode > $dir/augment/gen_out.csv
python augment_scripts/generate_mol_list.py --infile $dir/augment/gen_out.csv --outfile $dir/augment/mols_between.txt --old_list $dir/augment/mols_old.txt
end=$(date +'%s')
echo Execution time was $(($end-$start)) seconds. >> $log_file


start=$(date +'%s')
echo "=========Make Predictions========"
echo "Make Predictions" >> $log_file
python $predictor1_path/$predictor1_exec --test_path $dir/augment/mols_between.txt --checkpoint_dir $predictor1 --preds_path $dir/augment/gen_evaluated.csv
end=$(date +'%s')
echo Execution time was $(($end-$start)) seconds. >> $log_file


start=$(date +'%s')
echo "=========Add Molecules to Pairlist========"
echo "Add Molecules to Pairlist" >> $log_file
python generate_augmented_pair_list.py --gen_out $dir/augment/gen_out.csv --gen_evaluated $dir/augment/gen_evaluated.csv --datafile $dir/augment/data.csv --molfile $dir/augment/mols.txt --pairfile $dir/augment/train_pairs.txt --threshold 0.78 --pairgen_cutoff 0.78 --min_mol_wt 50.0
end=$(date +'%s')
echo Execution time was $(($end-$start)) seconds. >> $log_file





# ADDITIONAL ITERATIONS=======================================================================


for i in 2 3 4 5
do
    prev_dir=$dir
    dir=checkpoints/$folder/${folder}_${timestamp}/iteration$i
    mkdir $dir
    
    echo "========Starting Iteration " $i "========" >> $log_file

    start=$(date +'%s')
    echo "=========Get Vocabulary========"
    echo "Get Vocabulary" >> $log_file
    python get_vocab.py <$prev_dir/augment/mols.txt > $prev_dir/augment/vocab.txt
    end=$(date +'%s')
    echo Execution time was $(($end-$start)) seconds. >> $log_file


    start=$(date +'%s')
    echo "=========Preprocess Data========"
    echo "Preprocess Data" >> $log_file
    python preprocess.py --train $prev_dir/augment/train_pairs.txt --vocab $prev_dir/augment/vocab.txt --ncpu 16 < $prev_dir/augment/train_pairs.txt
    end=$(date +'%s')
    echo Execution time was $(($end-$start)) seconds. >> $log_file
    
    mkdir $dir/train_processed
    mv tensor* $dir/train_processed/
    

    start=$(date +'%s')
    echo "=========Train Model========"
    echo "Train Model" >> $log_file
    mkdir $dir/models/
    python gnn_train.py --train $dir/train_processed/ --vocab $prev_dir/augment/vocab.txt --save_dir $dir/models/ --epoch $n_epochs #--load_dir $prev_dir/models/ --load_epoch $(($n_epochs-1))
    end=$(date +'%s')
    echo Execution time was $(($end-$start)) seconds. >> $log_file
    
    mkdir $dir/augment
    rm -rf $dir/augment/mols_old.txt
    echo "SMILES" >> $dir/augment/mols_old.txt
    cat $prev_dir/augment/mols.txt >> $dir/augment/mols_old.txt

#     start=$(date +'%s')
#     echo "=========Make Predictions: Model 1========"
#     echo "Make Predictions: Model 1" >> $log_file
#     python $predictor1_path/$predictor1_exec --test_path $dir/augment/mols_old.txt --checkpoint_dir $predictor1 --preds_path $dir/augment/mols_pred_model1.csv
#     end=$(date +'%s')
#     echo Execution time was $(($end-$start)) seconds. >> $log_file

#     start=$(date +'%s')
#     echo "=========Make Predictions: SA Score========"
#     echo "Make Predictions: SA Score" >> $log_file
#     python $sa_score_path/$sa_score_exec --infile $dir/augment/mols_old.txt --outfile $dir/augment/mols_pred_sa_score.csv
#     end=$(date +'%s')
#     echo Execution time was $(($end-$start)) seconds. >> $log_file
    
    start=$(date +'%s')
    echo "=========Select Molecules for Generation========"
    echo "Select Molecules for Generation" >> $log_file
    python augment_scripts/select_mols.py --molfile $dir/augment/mols_old.txt --outfile $dir/augment/gen_in.csv #--screen_file1 $dir/augment/mols_pred_model1.csv #--screen_file2 $dir/augment/mols_pred_sa_score.csv
    end=$(date +'%s')
    echo Execution time was $(($end-$start)) seconds. >> $log_file
    
    
    start=$(date +'%s')
    echo "=========Generate Molecules========"
    echo "Generate Molecules" >> $log_file
    python decode.py --test $dir/augment/gen_in.csv --vocab $prev_dir/augment/vocab.txt --model $dir/models/$fin_model --num_decode $num_decode > $dir/augment/gen_out.csv
    python augment_scripts/generate_mol_list.py --infile $dir/augment/gen_out.csv --outfile $dir/augment/mols_between.txt --old_list $dir/augment/mols_old.txt
    end=$(date +'%s')
    echo Execution time was $(($end-$start)) seconds. >> $log_file


    start=$(date +'%s')
    echo "=========Make Predictions========"
    echo "Make Predictions" >> $log_file
    python $predictor1_path/$predictor1_exec --test_path $dir/augment/mols_between.txt --checkpoint_dir $predictor1 --preds_path $dir/augment/gen_evaluated.csv
    end=$(date +'%s')
    echo Execution time was $(($end-$start)) seconds. >> $log_file


    start=$(date +'%s')
    echo "=========Add Molecules to Pairlist========"
    echo "Add Molecules to Pairlist" >> $log_file
    python generate_augmented_pair_list.py --gen_out $dir/augment/gen_out.csv --gen_evaluated $dir/augment/gen_evaluated.csv --datafile $dir/augment/data.csv --molfile $dir/augment/mols.txt --pairfile $dir/augment/train_pairs.txt --threshold 0.78 --pairgen_cutoff 0.78 --min_mol_wt 50.0
    end=$(date +'%s')
    echo Execution time was $(($end-$start)) seconds. >> $log_file


done

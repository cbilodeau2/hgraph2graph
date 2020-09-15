# Hierarchical Generation of Molecular Graphs using Structural Motifs

Our paper is at https://arxiv.org/pdf/2002.03230.pdf

## Installation
First install the dependencies via conda:
 * PyTorch >= 1.0.0
 * networkx
 * RDKit
 * numpy
 * Python >= 3.6
 * Python < 3.8

And then run `pip install .`

## Molecule Generation
The molecule generation code is in the `generation/` folder.

## Graph translation Data Format
* The training file should contain pairs of molecules (molA, molB) that are similar to each other but molB has better chemical properties. Please see `data_original/qed/train_pairs.txt`.
* The test file is a list of molecules to be optimized. Please see `data_original/qed/test.txt`.

## Graph translation training procedure
1. Extract substructure vocabulary from a given set of molecules:
```
python get_vocab.py < data_original/qed/mols.txt > vocab.txt
```
Please replace `data_original/qed/mols.txt` with your molecules data file.

2. Preprocess training data:
```
python preprocess.py --train data_original/qed/train_pairs.txt --vocab data_original/qed/vocab.txt --ncpu 16 < data_original/qed/train_pairs.txt
mkdir train_processed
mv tensor* train_processed/
```
Please replace `--train` and `--vocab` with training and vocab file.

3. Train the model:
```
mkdir models/
python gnn_train.py --train train_processed/ --vocab data_original/qed/vocab.txt --save_dir models/ 
```

4. Make prediction on your lead compounds (you can use any model checkpoint, here we use model.5 for illustration)
```
python decode.py --test data_original/qed/valid.txt --vocab data_original/qed/vocab.txt --model models/model.5 --num_decode 20 > results.csv
```

## Using the scripts/notebooks

```
sh run.sh
sh evaluate.sh
```
You can then use Test_Summary.ipynb to calculate standard statistics for the test sets. Be sure to check that the paths are correct for your usage (set in each of the shell scripts). You also need to train a separate chemprop or other predictor model.


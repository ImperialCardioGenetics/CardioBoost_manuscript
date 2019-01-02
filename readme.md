# CardioBoost
This repository includes all the datasets and scripts to reproduce the machine learning training and evalution in the manuscript [Disease-specific variant pathogenicity prediction significantly improves clinical variant interpretation in inherited cardiac conditions] (URL).

## Dependencies

## Overview of the repository

### Collect features for variant classification

We wrote a handy package to collect the features relevant for variant classification in the project. 

Usage:
 - ../script/collect_feature/maste.py cm cm_varaint.txt cm_variant.vcf

The second argument denotes the disease type: 'cm' for cardiomyopathies and 'arm' for arrhythmias.

The third argument denotes a list of variants in table format. For example: 
| CHROM |POS  |REF  |ALT|SYMBOL|pathogenic|
|--|--|--|--|--|--|
|3 | 38655272 | C | T | SCN5A | 1| .

The forth argument denotes a list of variants in [VCF format] (https://en.wikipedia.org/wiki/Variant_Call_Format). The ID field of VCF format is the concatenate genomic coordinate of the variants. For example:
| #CHROM |POS  |ID|REF  |ALT|QUAL|FILTER|INFO|
|--|--|--|--|--|--|--|--|--|
|3 | 38655272 | 3_38655272_C_T| C | T | . | .|.| .


### Variant Classification on Cardiomyopathy and Evalutation 

To **preprocess** the collected features (missing values imputation and normalization):
- ./script/cardiomyopathy/ml/preprocess.Rmd 

To **train machine learning models** using nested cross-validation: optimising hyperparameters in inner CV loop and evaluating on outer CV loop

Set up how to split dataset in nested CV (i.e., fix the usage of samples in CV loops in the training of different models)
- ./script/cardiomyopathy/ml/tunecontrol.R

Scripts for training different machine learning models:
 - GLMNET: ./script/cardiomyopathy/ml/cvglmnet.R
 - Gradient Boosting: ./script/cardiomyopathy/ml/gbm.R
 - Random Forest: ./script/cardiomyopathy/ml/rf.R 
 - AdaBoost: ./script/cardiomyopathy/ml/ada.R
 - K Nearest Neighbours: ./script/cardiomyopathy/ml/knn.R
 - XGBoost (fast implementation of Gradient Boosting):  ./script/cardiomyopathy/ml/xgboost.R
 - Bayesian Average Boosting Tree:  ./script/cardiomyopathy/ml/bartMachine.R
 - Decision Tree: ./script/cardiomyopathy/ml/cart.R
 - Support Vector Machine with radial basis function kernel: ./script/cardiomyopathy/ml/svm_RBF.R

To **evaluate** the results of training and **select** the best trained model:
- ./script/cardiomyopathy/analysis/evaluation_training.Rmd

To **evaluate** the prediction performances on **out-of-sample** datasets:
- Hold-out test data: ./script/cardiomyopathy/analysis/evaluation_holdout_test.Rmd
- Additional test data: ./script/cardiomyopathy/analysis/evaluation_additional_test.Rmd

To **generate** predictions for all possible rare missense variants on disease-associated genes:
- ./script/cardiomyopathy/ml/cm_prediction_all_possible_variants.R

To **evaluate** the etiological fraction:
- ./script/cardiomyopathy/analysis/EF_calculation_cardioboost.Rmd
- ./script/cardiomyopathy/analysis/EF_calculation_mcap.Rmd
- ./script/cardiomyopathy/analysis/EF_calculation_revel.Rmd

To **evaluate** survival outcomes:
- ./script/cardiomyopathy/analysis/survival_analysis.Rmd

To **compare** the classification performances on indirectly seen and unseen data:
- ./script/cardiomyopathy/analysis/evaluation_unseen_seen.Rmd

### Variant Classification on Arrhythmia and Evalutation 

To **preprocess** the collected features (missing values imputation and normalization):
- ./script/arrhythmia/ml/preprocess.Rmd 

To **train machine learning models** using nested cross-validation: optimising hyperparameters in inner CV loop and evaluating on outer CV loop

Set up how to split dataset in nested CV (i.e., fix the usage of samples in CV loops in the training of different models)
- ./script/arrhythmia/ml/tunecontrol.R

Scripts for training different machine learning models:
 - GLMNET: ./script/arrhythmia/ml/cvglmnet.R
 - Gradient Boosting: ./script/arrhythmia/ml/gbm.R
 - Random Forest: ./script/arrhythmia/ml/rf.R 
 - AdaBoost: ./script/arrhythmia/ml/ada.R
 - K Nearest Neighbours: ./script/arrhythmia/ml/knn.R
 - XGBoost (fast implementation of Gradient Boosting):  ./script/arrhythmia/ml/xgboost.R
 - Bayesian Average Boosting Tree:  ./script/arrhythmia/ml/bartMachine.R
 - Decision Tree: ./script/arrhythmia/ml/cart.R
 - Support Vector Machine with radial basis function kernel: ./script/arrhythmia/ml/svm_RBF.R

To **evaluate** the results of training and **select** the best trained model:
- ./script/arrhythmia/analysis/evaluation_training.Rmd

To **evaluate** the prediction performances on **out-of-sample** datasets:
- Hold-out test data: ./script/arrhythmia/analysis/evaluation_holdout_test.Rmd
- Additional test data: ./script/arrhythmia/analysis/evaluation_additional_test.Rmd

To **generate** predictions for all possible rare missense variants on disease-associated genes:
- ./script/arrhythmia/ml/cm_prediction_all_possible_variants.R

To **compare** the classification performances on indirectly seen and unseen data:
- ./script/arrhythmia/analysis/evaluation_unseen_seen.Rmd

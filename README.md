# CardioBoost
This repository includes all the datasets and scripts to reproduce the machine learning training, evaluation, and analysis in the manuscript [Disease-specific variant pathogenicity prediction significantly improves clinical variant interpretation in inherited cardiac conditions](bioRxiv:URL).

## Dependencies

To install the required R packages:
`Rscript ./script/install_requiredRpkg.R`

## Overview of the repository

### Variant Classification on Cardiomyopathy and Evaluation

To **preprocess** the collected features and data splitting (missing values imputation and normalisation):
`./script/cardiomyopathy/ml/preprocess_full.Rmd`

To **train machine learning models** using nested cross-validation (CV): optimising hyper-parameters in inner CV loop and evaluating on outer CV loop

Set up how to split dataset in nested CV 
`Rscript ./script/cardiomyopathy/ml/tunecontrol.R`

Scripts for training different machine learning models:
 - GLMNET: `Rscript ./script/cardiomyopathy/ml/cvglmnet.R`
 - Gradient Boosting: `Rscript ./script/cardiomyopathy/ml/gbm.R`
 - Random Forest: `Rscript ./script/cardiomyopathy/ml/rf.R`
 - AdaBoost: `Rscript ./script/cardiomyopathy/ml/ada.R`
 - K Nearest Neighbours: `Rscript ./script/cardiomyopathy/ml/knn.R`
 - XGBoost (fast implementation of Gradient Boosting):  `Rscript ./script/cardiomyopathy/ml/xgboost.R`
 - Bayesian Average Boosting Tree:  `Rscript ./script/cardiomyopathy/ml/bartMachine.R`
 - Decision Tree: `Rscript ./script/cardiomyopathy/ml/cart.R`
 - Support Vector Machine with radial basis function kernel: `Rscript./script/cardiomyopathy/ml/svm_RBF.R`

To **evaluate** the results of training and **select** the best trained model:
`'./script/cardiomyopathy/analysis/evaluation_training.Rmd'`

To **tune** and **train** the selected model:
`Rscript ./script/cardiomyopathy/ada_tune.R`
`Rscript ./script/cardiomyopathy/train_ada.R`

To **evaluate** the prediction performances on **out-of-sample** datasets:
- Hold-out test data:
`'./script/cardiomyopathy/analysis/evaluation_holdout_test.Rmd'`
- Additional test data:
`'./script/cardiomyopathy/analysis/evaluation_additional_test.Rmd'`

To **generate** predictions for all possible rare missense variants on disease-associated genes:
`Rscript ./script/cardiomyopathy/ml/cm_prediction_all_possible_variants.R`

To **evaluate** the etiological fraction:
`'./script/cardiomyopathy/analysis/OR_calculation_cardioboost.Rmd'`
`'./script/cardiomyopathy/analysis/OR_calculation_mcap.Rmd'`
`'./script/cardiomyopathy/analysis/OR_calculation_revel.Rmd'`

To **evaluate** survival outcomes:
`'./script/cardiomyopathy/analysis/survival_analysis.Rmd'`

To **compare** the classification performances on indirectly seen and unseen data:
`'./script/cardiomyopathy/analysis/evaluation_unseen_seen.Rmd'`

### Variant Classification on Arrhythmia and Evaluation

To **preprocess** the collected features (missing values imputation and normalisation):
` ./script/arrhythmia/ml/preprocess_full.Rmd`

To **train machine learning models** using nested cross-validation: optimising hyper-parameters in inner CV loop and evaluating on outer CV loop

Set up how to split dataset in nested CV (i.e., fix the usage of samples in CV loops in the training of different models)
`Rscript ./script/arrhythmia/ml/tunecontrol.R`

Scripts for training different machine learning models:
 - GLMNET: `Rscript ./script/arrhythmia/ml/cvglmnet.R`
 - Gradient Boosting: `Rscript ./script/arrhythmia/ml/gbm.R`
 - Random Forest: `Rscript ./script/arrhythmia/ml/rf.R`
 - AdaBoost: `Rscript ./script/arrhythmia/ml/ada.R`
 - K Nearest Neighbours: `Rscript ./script/arrhythmia/ml/knn.R`
 - XGBoost (fast implementation of Gradient Boosting):  `Rscript ./script/arrhythmia/ml/xgboost.R`
 - Bayesian Average Boosting Tree: `Rscript  ./script/arrhythmia/ml/bartMachine.R`
 - Decision Tree: `Rscript ./script/arrhythmia/ml/cart.R`
 - Support Vector Machine with radial basis function kernel: `Rscript ./script/arrhythmia/ml/svm_RBF.R`

To **evaluate** the results of training and **select** the best trained model:
`'./script/arrhythmia/analysis/evaluation_training.Rmd'`

To **evaluate** the prediction performances on **out-of-sample** datasets:
- Hold-out test data:
`'./script/arrhythmia/analysis/evaluation_holdout_test.Rmd'`
- Additional test data:
`'./script/arrhythmia/analysis/evaluation_additional_test.Rmd'`

To **generate** predictions for all possible rare missense variants on disease-associated genes:
`Rscript ./script/arrhythmia/ml/arm_prediction_all_possible_variants.R`

To **compare** the classification performances on indirectly seen and unseen data:
`' ./script/arrhythmia/analysis/evaluation_unseen_seen.Rmd'`

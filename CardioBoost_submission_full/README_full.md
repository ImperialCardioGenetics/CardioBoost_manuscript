# CardioBoost
This repository includes all the datasets and scripts to reproduce the machine learning training and evalution in the manuscript [Disease-specific variant pathogenicity prediction significantly improves clinical variant interpretation in inherited cardiac conditions](Bioaxiv:URL).

## Dependencies
To install the required Python packages:
`pip install -r ./script/requirements.txt `

To install the required R packages:
`Rscript ./script/install_requiredRpkg.R`

## Overview of the repository

### Collect features for variant classification (Private)

We wrote a handy package to collect the features relevant for variant classification in the project.

Usage:
 ` python ./script/collect_feature/maste.py cm cm_varaint.txt cm_variant.vcf`

The second argument denotes the disease type: 'cm' for cardiomyopathies and 'arm' for arrhythmias.

The third argument denotes a list of variants in table format. For example:
| CHROM |POS  |REF  |ALT|SYMBOL|pathogenic|
|--|--|--|--|--|--|
|3 | 38655272 | C | T | SCN5A | 1| .

The forth argument denotes a list of variants in [VCF format](https://en.wikipedia.org/wiki/Variant_Call_Format). The ID field of VCF format is the concatenate genomic coordinate of the variants. For example:
| #CHROM |POS  |ID|REF  |ALT|QUAL|FILTER|INFO|
|--|--|--|--|--|--|--|--|--|
|3 | 38655272 | 3_38655272_C_T| C | T | . | .|.| .


### Variant Classification on Cardiomyopathy and Evalutation

To **preprocess** the collected features (missing values imputation and normalization):
`Rscript -e "rmarkdown::render('./script/cardiomyopathy/ml/preprocess.Rmd', clean=TRUE)"`

To **train machine learning models** using nested cross-validation: optimising hyperparameters in inner CV loop and evaluating on outer CV loop

Set up how to split dataset in nested CV (i.e., fix the usage of samples in CV loops in the training of different models)
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
`Rscript -e "rmarkdown::render('./script/cardiomyopathy/analysis/evaluation_training.Rmd',clean=TRUE)"`

To **evaluate** the prediction performances on **out-of-sample** datasets:
- Hold-out test data:
`Rscript -e "rmarkdown::render('./script/cardiomyopathy/analysis/evaluation_holdout_test.Rmd',clean=TRUE)"`
- Additional test data:
`Rscript -e "rmarkdown::render('./script/cardiomyopathy/analysis/evaluation_additional_test.Rmd',clean=TRUE)"`

To **generate** predictions for all possible rare missense variants on disease-associated genes:
`Rscript ./script/cardiomyopathy/ml/cm_prediction_all_possible_variants.R`

To **evaluate** the etiological fraction:
`Rscript -e "rmarkdown::render('./script/cardiomyopathy/analysis/EF_calculation_cardioboost.Rmd',clean=TRUE)"`
`Rscript -e "rmarkdown::render('./script/cardiomyopathy/analysis/EF_calculation_mcap.Rmd',clean=TRUE)"`
`Rscript -e "rmarkdown::render('./script/cardiomyopathy/analysis/EF_calculation_revel.Rmd',clean=TRUE)"`

To **evaluate** survival outcomes:
`Rscript -e "rmarkdown::render('./script/cardiomyopathy/analysis/survival_analysis.Rmd',clean=TRUE)"`

To **compare** the classification performances on indirectly seen and unseen data:
`Rscript -e "rmarkdown::render('./script/cardiomyopathy/analysis/evaluation_unseen_seen.Rmd',clean=TRUE)"`

### Variant Classification on Arrhythmia and Evalutation

To **preprocess** the collected features (missing values imputation and normalization):
- ./script/arrhythmia/ml/preprocess.Rmd

To **train machine learning models** using nested cross-validation: optimising hyperparameters in inner CV loop and evaluating on outer CV loop

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
`Rscript -e "rmarkdown::render('./script/arrhythmia/analysis/evaluation_training.Rmd',clean=TRUE)"`

To **evaluate** the prediction performances on **out-of-sample** datasets:
- Hold-out test data:
`Rscript -e "rmarkdown::render('./script/arrhythmia/analysis/evaluation_holdout_test.Rmd',clean=TRUE)"`
- Additional test data:
`Rscript -e "rmarkdown::render('./script/arrhythmia/analysis/evaluation_additional_test.Rmd',clean=TRUE)"`

To **generate** predictions for all possible rare missense variants on disease-associated genes:
`Rscript ./script/arrhythmia/ml/cm_prediction_all_possible_variants.R`

To **compare** the classification performances on indirectly seen and unseen data:
`Rscript -e "rmarkdown::render(' ./script/arrhythmia/analysis/evaluation_unseen_seen.Rmd',clean=TRUE)"`

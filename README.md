# ReSSI

Code for **Resilient Semi-Supervised Inference with Heterogeneous Unlabeled Data**, by Mengyuan Wang, Chengde Qian, Haojie Ren, and Changliang Zou.

We provide an R implementation of ReSSI for semi-supervised estimation and inference with contaminated unlabeled data, together with simulation and real-data experiments.

## Instructions to run the experiments

Run the following commands from the repository directory:

```sh
# Mean estimation, linear regression, and logistic regression
Rscript simulation.R

# Promoter-expression data
Rscript gene_expression.R

# Census data
Rscript census.R
```
# CRISPR

Validate extension of the `dce` methodology for latent confounding adjustment on GTEx data. Here, we try to show robustness against confounding by comparing the output of the method when applied on the original data and when we add the provided confounding covariates.

* Data source:
    * https://gtexportal.org/home/index.html 
    * https://gtexportal.org/home/datasets


## Raw Data

In order to run this workflow, create folder called raw_data and download and unzip the following files:
	* gencode.v26.GRCh38.genes.gtf
	* GTEx_Analysis_v8_sQTL_covariates
	* GTEx_Analysis_v8_eQTL_expression_matrices

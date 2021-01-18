# dce

Compute differential causal effects on (biological) networks.


## Installation

Install the latest stable version from Bioconductor:
```r
BiocManager::install("dce")
```

Install the latest development version from GitHub:
```r
remotes::install_github("cbg-ethz/dce")
```


## Project structure

* `dce`: R package
* `analysis_workflows`: Snakemake workflows for all investigations in publication
    * `crispr_benchmark`: Real-life data validation
    * `ovarian_cancer`: How does Ovarian Cancer dysregulate pathways?
    * `synthetic_benchmark`: Synthetic data validation
    * `tcga_pipeline`: Compute effects for loads of data from TCGA

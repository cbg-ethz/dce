# dce

[![lint](https://github.com/cbg-ethz/dce/workflows/lint/badge.svg)](https://github.com/cbg-ethz/dce/actions)
[![check-bioc](https://github.com/cbg-ethz/dce/workflows/check-bioc/badge.svg)](https://github.com/cbg-ethz/dce/actions)
[![pkgdown](https://github.com/cbg-ethz/dce/workflows/pkgdown/badge.svg)](https://github.com/cbg-ethz/dce/actions)

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

* `.`: R package
* `inst/scripts/`: Snakemake workflows for all investigations in publication
    * `crispr_benchmark`: Real-life data validation
    * `ovarian_cancer`: How does Ovarian Cancer dysregulate pathways?
    * `synthetic_benchmark`: Synthetic data validation
    * `tcga_pipeline`: Compute effects for loads of data from TCGA


## Development notes

* Check package locally:
    * `Rscript -e "lintr::lint_package()"`
    * `Rscript -e "devtools::check(error_on = 'warning')"`
    * `R CMD BiocCheck`
* Documentation
    * Build locally: `Rscript -e "pkgdown::build_site()"`
    * Deploy: `Rscript -e "pkgdown::deploy_to_branch(new_process = FALSE)"`
* Bioconductor
    * The `bioc` branch stores changes specific to Bioconductor releases

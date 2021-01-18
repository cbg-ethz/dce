# Benchmarking

## Automated

All important benchmarks shall be defined in `config/config.yaml`.
They can then be executing as follows:
```bash
$ snakemake -pr -j 1
```

## Manual

Specific configurations can be easily tested:
```bash
$ Rscript workflow/scripts/main.R --variable node.num --values "25,50,100,200" --methods "dce,dce.lm" --replicates 10
$ Rscript workflow/scripts/plotting.R
```

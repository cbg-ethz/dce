# Benchmarking

## Automated

All important benchmarks shall be defined in `config.yaml`.
They can then be executing as follows:
```bash
$ snakemake -pr -j 1
```

## Manual

Specific configurations can be easily tested:
```bash
$ Rscript main.R --variable node.num --values "25,50,100,200" --methods "dce,dce.lm" --replicates 10
$ Rscript plotting.R
```

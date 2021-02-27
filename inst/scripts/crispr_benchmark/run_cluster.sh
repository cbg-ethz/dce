bsub \
  -N \
  -R 'rusage[mem=2000]' \
  -q normal.24h \
  -oo snake.out -eo snake.err \
snakemake \
  --profile lsf \
  -pr \
  -j 100 \
  "$@"

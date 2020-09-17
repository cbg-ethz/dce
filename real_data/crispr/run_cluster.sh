workdir="$(grep ^workdir config.yaml | cut -d ' ' -f 2 | tr -d "'")" # oh no, very bad
mkdir -p "./$workdir"
cp -v "lsf.yaml" "./$workdir"

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

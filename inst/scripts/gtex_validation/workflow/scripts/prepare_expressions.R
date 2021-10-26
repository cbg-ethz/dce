library(stringr)
print('Loading gene expressions for all tissues.')
df = read.csv(snakemake@input$all_expressions_file, skip = 2, header = TRUE, sep='\t')
df = read.csv(snakemake@input$all_expressions_file, skip = 2, header = TRUE, sep='\t')
rownames(df) = df[,1]
df = df[,-c(1,2)]

annotations = read.csv(snakemake@input$annotations_file, sep='\t')

transform_tissue_name = function(s){
  s = str_replace_all(s, '([[:punct:]])', ' ')
  s = gsub("\\s+", " ", str_trim(s))
  s = str_replace_all(s, ' ', '_')
  return(s)
}

index = sapply(colnames(df), function(sample){
  which.max(annotations$SAMPID == str_replace_all(sample, '\\.', '-'))
  })
tissues = annotations$SMTSD[index]

for(tissue in unique(annotations$SMTSD)){
  print(paste0('Processing data for tissue: ', tissue))
  data = df[, tissue == tissues]
  colnames(data) = sapply(colnames(data), function(s){paste0('GTEX.', strsplit(s, '\\.')[[1]][2])})
  
  dir.create(file.path(snakemake@output$expressions_dir), showWarnings = FALSE)
  filename = paste0(snakemake@output$expressions_dir, '/', transform_tissue_name(tissue), '.csv')
  write.table(data, file=filename, sep='\t')
}

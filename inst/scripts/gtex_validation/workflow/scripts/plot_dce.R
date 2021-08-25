library(ggplot2)
library(cowplot)

#out = mget(load(snakemake@input[[1]], envir=(NE. <- new.env())), envir=NE.)
out = mget(load("/Users/cevidd/Desktop/dce/inst/scripts/gtex_validation/results/output/Breast_Mammary_Tissue#Lung.Rdata", envir=(NE. <- new.env())), envir=NE.)
x1 = out$normal_fit_with_deconf$dce
x2 = out$extended_fit_with_deconf$dce
shared_genes = intersect(colnames(x1), colnames(x2))
x1 = x1[shared_genes, shared_genes]
x2 = x2[shared_genes, shared_genes]
x1 = c(x1)[is.na(c(x1))==FALSE]
x2 = c(x2)[is.na(c(x2))==FALSE]
y1 = out$normal_fit_no_deconf$dce
y2 = out$extended_fit_no_deconf$dce
shared_genes = intersect(colnames(y1), colnames(y2))
y1 = y1[shared_genes, shared_genes]
y2 = y2[shared_genes, shared_genes]
y1 = c(y1)[is.na(c(y1))==FALSE]
y2 = c(y2)[is.na(c(y2))==FALSE]


gg1 = ggplot(data.frame(x=x1, y=x2), aes(x=x, y=y))+
  geom_point(size=0.8, alpha=0.7)+
  geom_abline(color='red', linetype='dashed') +
  xlab('DCE for original pathway') + 
  ylab('DCE for extended pathway') +
  theme_light() +
  theme(legend.position='none', 
        axis.text.y=element_text(size=13), 
        axis.text.x=element_text(size=13), 
        axis.title.x=element_text(size=14), 
        axis.title.y=element_text(size=14))
gg2 = ggplot(data.frame(x=y1, y=y2), aes(x=x, y=y))+
  geom_point(size=0.8, alpha=0.7)+
  geom_abline(color='red', linetype='dashed') +
  xlab('DCE for original pathway') + 
  ylab('DCE for extended pathway') +
  theme_light() +
  theme(legend.position='none', 
        axis.text.y=element_text(size=13), 
        axis.text.x=element_text(size=13), 
        axis.title.x=element_text(size=14),
        axis.title.y=element_text(size=14))
gg = plot_grid(gg1, gg2, ncol=1, 
               labels=c('with deconfounding', 'no deconfounding'), label_size = 16)

ggsave(gg, filename=snakemake@output[[1]], width=6, height=6)
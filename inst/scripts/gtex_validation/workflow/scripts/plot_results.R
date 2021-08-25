library(ggplot2)
library(cowplot)
load(snakemake@input[[1]])

gg1 = ggplot(rbind(data.frame(value=mse_no_deconf^0.5, deconfounding='without deconfounding'), 
             data.frame(value=mse_with_deconf^0.5, deconfounding='with deconfounding')), aes(x=deconfounding, y=value))+
  geom_boxplot(aes(fill=deconfounding)) + theme_light() +
  theme(legend.position='none', 
        axis.text.y=element_text(size=12), 
        axis.text.x=element_text(size=18), 
        axis.title.y=element_text(size=15)) + 
  ylab('difference in DCEs for original and extended pathways') + xlab('') #+ ylim(c(0.1, 0.3)) 
ggsave(plot=gg1, filename=snakemake@output[[1]], width=7, height=6)

gg2 = ggplot(rbind(data.frame(value=cor_no_deconf, deconfounding='without deconfounding'), 
             data.frame(value=cor_with_deconf, deconfounding='with deconfounding')), aes(x=deconfounding, y=value))+
  geom_boxplot(aes(fill=deconfounding)) + theme_light() +
  theme(legend.position='none', 
        axis.text.y=element_text(size=12), 
        axis.text.x=element_text(size=18), 
        axis.title.y=element_text(size=15)) + 
  ylab('correlation of DCEs for original and extended pathways') + xlab('') #+ ylim(c(0.1, 0.3)) 
ggsave(plot=gg2, filename=snakemake@output[[2]], width=7, height=6)


ggsave(plot_grid(gg1,gg2), filename=snakemake@output[[3]], width=14, height=6)
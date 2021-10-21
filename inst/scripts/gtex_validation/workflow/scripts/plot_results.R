###
# Plot comparison of performance measures between
# original and extended pathway.
###


library(ggplot2)
library(cowplot)
load(snakemake@input[[1]])

gg1 = ggplot(rbind(data.frame(value=mse_no_deconf^0.5, deconfounding='without deconfounding'),
             data.frame(value=mse_with_deconf^0.5, deconfounding='with deconfounding')), aes(x=deconfounding, y=value))+
  geom_boxplot(aes(fill=deconfounding)) + theme_light() +
  theme(legend.position='none',
        axis.text.y=element_text(size=12),
        axis.text.x=element_text(size=18),
        axis.title.y=element_text(size=14)) +
  ylab('difference in DCEs for original and extended pathway') + xlab('') #+ ylim(c(0.1, 0.3))
ggsave(plot=gg1, filename=snakemake@output[[1]], width=8, height=5.5)

print(c('mean - with deconfounding:', mean(cor_with_deconf)))
print(c('sd - with deconfounding:', sd(cor_with_deconf)))
print(c('mean - no deconfounding:', mean(cor_no_deconf)))
print(c('sd - no deconfounding:', sd(cor_no_deconf)))
print(c('significance:', wilcox.test(cor_with_deconf, cor_no_deconf, paired = TRUE, alternative = "greater")$p.value))

gg2 = ggplot(rbind(data.frame(value=cor_no_deconf, deconfounding='without deconfounding'),
             data.frame(value=cor_with_deconf, deconfounding='with deconfounding')), aes(x=deconfounding, y=value))+
  geom_boxplot(aes(fill=deconfounding)) + theme_light() +
  theme(legend.position='none',
        axis.text.y=element_text(size=12),
        axis.text.x=element_text(size=18),
        axis.title.y=element_text(size=14)) +
  ylab('correlation of DCEs for original and extended pathway') + xlab('') #+ ylim(c(0.1, 0.3))
ggsave(plot=gg2, filename=snakemake@output[[2]], width=8, height=5.5)


ggsave(plot_grid(gg1,gg2), filename=snakemake@output[[3]], width=16, height=5.5)

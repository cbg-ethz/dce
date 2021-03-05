library(tidyverse)
library(magrittr)
library(graph)
library(naturalsort)

devtools::load_all("../../")

source("workflow/scripts/helper_functions.R")
source("workflow/scripts/models.R")
source("workflow/scripts/performance_measures.R")
source("workflow/scripts/LDGM.R")
source("workflow/scripts/FastGGM.R")

methods <- c('auto','kim','cluster')
runs <- 1
latents <- c(0,5,10)
acc <- matrix(NA,runs*length(latents),length(methods)+2)
colnames(acc) <- c(methods, 'seed', 'latent')
seed.list <- sample(seq_len(10^9), runs*length(latents))
for (lat in 1:length(latents)) {
    latent <- latents[lat]
    for (run in 1:runs) {
        idx <- (lat-1)*runs + run
        rng.seed <- seed.list[idx]
        set.seed(rng.seed)
        graphs <- generate.random.graphs(100, 1, 0.5)
        wt.graph <- graphs$wt
        mt.graph <- graphs$mt
        # generate data
        pop.size <- 10000
        wt.X <- simulate_data(wt.graph, n = 200, dist_dispersion = 1, dist_mean = 100, pop_size = pop.size, latent = latent)
        mt.X <- simulate_data(mt.graph, n = 200, dist_dispersion = 1, dist_mean = 100, pop_size = pop.size, latent = latent)
        # library size difference
        # lib.size.range <- 10
        # lib.size.mean <- (lib.size.range+1)/2
        # lib.size.sd <- lib.size.range/10*2
        # ptruncnorm <- dnorm(1:lib.size.range,lib.size.mean,lib.size.sd)/(pnorm(lib.size.range+1,lib.size.mean,lib.size.sd)-pnorm(0,lib.size.mean,lib.size.sd))
        # lib.size.gtn <- sample(1:lib.size.range,nrow(wt.X)+nrow(mt.X),replace=TRUE,prob=ptruncnorm)
        # wt.samples <- mt.samples <- 200
        # wt.X <- wt.X*lib.size.gtn[seq_len(wt.samples)]
        # mt.X <- mt.X*lib.size.gtn[(wt.samples+1):(wt.samples+mt.samples)]
        acc[idx,1] <- estimate_latent_count(compute.tpm(wt.X), compute.tpm(mt.X), 'auto') - latent
        acc[idx,2] <- estimate_latent_count(compute.tpm(wt.X), compute.tpm(mt.X), 'kim') - latent
        acc[idx,3] <- estimate_latent_count(compute.tpm(wt.X), compute.tpm(mt.X), 'cluster') - latent
        acc[idx,4] <- rng.seed
        acc[idx,5] <- latent
        cat('.')
    }
}
write.csv(acc, file = paste0('sub_', as.numeric(Sys.time()), '.csv'))

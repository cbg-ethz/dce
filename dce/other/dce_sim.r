source("dce/main.R")
source("dce/utils.R")

library(tidyverse)
library(purrr)
library(graph)
library(pcalg)
library(assertthat)

m <- numeric(2)
n <- as.numeric(commandArgs(TRUE)[1])
m[1] <- as.numeric(commandArgs(TRUE)[2])
m[2] <- as.numeric(commandArgs(TRUE)[3])
sd <- as.numeric(commandArgs(TRUE)[4])
runs <- as.numeric(commandArgs(TRUE)[5])
perturb <- as.numeric(commandArgs(TRUE)[6])
p <- as.numeric(commandArgs(TRUE)[7])

## n <- 100; m <- c(1000, 100); sd <- 1; runs <- 1; perturb <- 0;

if (is.na(runs)) {
    runs <- 100 # simulation runs
}
if (is.na(perturb)) {
    perturb <- 0
}
if (is.na(p)) {
    p <- 0.2 # edge prob of the dag
}
## uniform limits:
lB <- -1
uB <- 1
## others:
#n <- 10 # number of nodes
#m <- c(100,100) # number of samples tumor and normal
#sd <- 0.1 # standard deviation for variable distributions
## the fraction of true pos (causal effects that are differential)
truepos <- 0.9
bsruns <- 100

acc <- array(0, c(runs,7,5),
             dimnames = list(runs = paste0("run_", seq_len(runs)),
                             methods = c("dce", "random", "bootstrap",
                                         "subsample", "full linear", "subsample2", "full boot"),
                             metrics = c("correlation", "distance", "metric3",
                                         "metric4", "metric5")))
gtnfeat <- array(0, c(runs, 6, 2),
                 dimnames = list(runs = paste0("run_", seq_len(runs)),
                                 features = c("avg children", "avg parents",
                                              "childless", "parentless",
                                              "maxpathlength", "density"),
                                 gtn = c("originial", "transitive closure")))


for (i in 1:runs) {
    cat(i)
                                        #for (j in 1:1000) {
                                        #   set.seed(j)
    normal <- randomDAG(n, p, lB, uB)
    dn <- rmvDAG_2(m[2], normal, normpars = c(0,sd))

    tumor <- newWeights(normal, lB, uB, truepos) # resample edge weights
    dt <- rmvDAG_2(m[1], tumor, normpars = c(0,sd))

    cn <- trueCov(normal)
    ct <- trueCov(tumor)

    gm <- as(normal, "matrix")
    gm[which(gm != 0)] <- 1

    gtc <- mnem:::mytc(gm) # transitively closed graph as matrix

    ## save features of gtn, which might correlate with accuracy:
    gtr <- nem::transitive.reduction(gm)
    gtnfeat[i, 1, 1] <- mean(apply(gm, 1, sum))
    gtnfeat[i, 2, 1] <- mean(apply(gm, 2, sum))
    gtnfeat[i, 3, 1] <- sum(apply(gm, 1, sum) == 0)
    gtnfeat[i, 4, 1] <- sum(apply(gm, 2, sum) == 0)
    diag(gtr) <- 1
    for (j in seq_len(n)) {
        gtr <- gtr%*%gtr
        gtr[which(gtr > 1)] <- 1
        if (all(gtr == gtc)) { break() }
    }
    gtnfeat[i, 5, 1] <- j+1
    gtnfeat[i, 6, 1] <- sum(gm)

    diag(gtc) <- 0

    gtnfeat[i, 1, 2] <- mean(apply(gtc, 1, sum))
    gtnfeat[i, 2, 2] <- mean(apply(gtc, 2, sum))
    gtnfeat[i, 3, 2] <- sum(apply(gtc, 1, sum) == 0)
    gtnfeat[i, 4, 2] <- sum(apply(gtc, 2, sum) == 0)
    gtnfeat[i, 5, 2] <- j+1
    gtnfeat[i, 6, 2] <- sum(gtc)

    ## ground truth:
    dcet <- (cn - ct)*gtc # gtn for differential causal effects
    dcegtn <- list(dce = dcet, graph = normal, dcefull = dcet)
    class(dcegtn) <- "dce"

    dcetb <- dcet
    dcetb[which(abs(dcet) > 0.5)] <- 1
    dcetb[which(abs(dcet) <= 0.5)] <- 0

    ## perturb network:
    normalOrg <- normal
    adjn <- as(normal, "matrix")
    if (perturb == 2) {
        randag <- randomDAG(n, p, lB, uB)
        adjn <- as(randag, "matrix")
    }
    if (perturb < 0) {
        adjn <- as(normal, "matrix")
        remedge <- sample(which(adjn != 0), floor(sum(adjn != 0)*abs(perturb)))
        adjn[remedge] <- 0
    }
    if (perturb > 0 & perturb != 2) {
        adjn <- as(normal, "matrix")
        addedge <- sample(which(adjn == 0 & upper.tri(adjn)), floor(sum(adjn != 0)*perturb))
        adjn[addedge] <- 1
    }
    adjn[which(adjn != 0)] <- 1
    normal <- as(adjn, "graphNEL")

    ## our inference:

    ## ## bootstrap/subsample2:
    ## dcei <- compute_differential_causal_effects(
    ##     normal, dn,
    ##     tumor, dt,
    ##     bootstrap = TRUE, runs = bsruns, replace = 0, frac = 0.5,
    ##     strap = 1
    ## )
    ## dceibs2 <- dcei

    ## dceib <- dcei$dce
    ## dceib[which(abs(dceib) > 0.5)] <- 1
    ## dceib[which(abs(dceib) <= 0.5)] <- 0

    ## dcei <- dcei$dce
    ## coridx <- which(dcet != 0 | dcei != 0)
    ## acc[i, 6, 1] <- cor(as.vector(dcet[coridx]), as.vector(dcei[coridx]), method = "p")
    ## acc[i, 6, 2] <- dist(rbind(as.vector(dcet), as.vector(dcei)))
    ## tp <- sum(dcetb == 1 & dceib == 1)
    ## fp <- sum(dcetb == 0 & dceib == 1)
    ## fn <- sum(dcetb == 1 & dceib == 0)
    ## tn <- sum(dcetb == 0 & dceib == 0)
    ## acc[i, 6, 3] <- tp/(tp+fn)
    ## acc[i, 6, 4] <- tn/(tn+fp)
    ## acc[i, 6, 5] <- (tp+tn)/(tn+fp+tp+fn)

    ## ## bootstrap:
    ## dcei <- compute_differential_causal_effects(
    ##     normal, dn,
    ##     tumor, dt,
    ##     bootstrap = TRUE, runs = bsruns, replace = 0, frac = 0.5
    ## )
    ## dceibs <- dcei

    ## dceib <- dcei$dce
    ## dceib[which(abs(dceib) > 0.5)] <- 1
    ## dceib[which(abs(dceib) <= 0.5)] <- 0

    ## dcei <- dcei$dce
    ## coridx <- which(dcet != 0 | dcei != 0)
    ## acc[i, 3, 1] <- cor(as.vector(dcet[coridx]), as.vector(dcei[coridx]), method = "p")
    ## acc[i, 3, 2] <- dist(rbind(as.vector(dcet), as.vector(dcei)))
    ## tp <- sum(dcetb == 1 & dceib == 1)
    ## fp <- sum(dcetb == 0 & dceib == 1)
    ## fn <- sum(dcetb == 1 & dceib == 0)
    ## tn <- sum(dcetb == 0 & dceib == 0)
    ## acc[i, 3, 3] <- tp/(tp+fn)
    ## acc[i, 3, 4] <- tn/(tn+fp)
    ## acc[i, 3, 5] <- (tp+tn)/(tn+fp+tp+fn)

    ## ## subsampling:
    ## dcei <- compute_differential_causal_effects(
    ##     normal, dn,
    ##     tumor, dt,
    ##     bootstrap = TRUE, runs = bsruns, replace = 1, frac = 1
    ## )
    ## dceiss <- dcei

    ## dceib <- dcei$dce
    ## dceib[which(abs(dceib) > 0.5)] <- 1
    ## dceib[which(abs(dceib) <= 0.5)] <- 0

    ## dcei <- dcei$dce
    ## coridx <- which(dcet != 0 | dcei != 0)
    ## acc[i, 4, 1] <- cor(as.vector(dcet[coridx]), as.vector(dcei[coridx]), method = "p")
    ## acc[i, 4, 2] <- dist(rbind(as.vector(dcet), as.vector(dcei)))
    ## tp <- sum(dcetb == 1 & dceib == 1)
    ## fp <- sum(dcetb == 0 & dceib == 1)
    ## fn <- sum(dcetb == 1 & dceib == 0)
    ## tn <- sum(dcetb == 0 & dceib == 0)
    ## acc[i, 4, 3] <- tp/(tp+fn)
    ## acc[i, 4, 4] <- tn/(tn+fp)
    ## acc[i, 4, 5] <- (tp+tn)/(tn+fp+tp+fn)

    ## full linear model:
    dcei <- fulllin(
        normal, dn,
        tumor, dt
    )
    dceifl <- dcei

    dceib <- dcei$dce
    dceib[which(abs(dceib) > 0.5)] <- 1
    dceib[which(abs(dceib) <= 0.5)] <- 0

    dcei <- dcei$dce
    coridx <- which(dcet != 0 | dcei != 0)
    acc[i, 5, 1] <- cor(as.vector(dcet[coridx]), as.vector(dcei[coridx]), method = "p")
    acc[i, 5, 2] <- dist(rbind(as.vector(dcet), as.vector(dcei)))
    tp <- sum(dcetb == 1 & dceib == 1)
    fp <- sum(dcetb == 0 & dceib == 1)
    fn <- sum(dcetb == 1 & dceib == 0)
    tn <- sum(dcetb == 0 & dceib == 0)
    acc[i, 5, 3] <- tp/(tp+fn)
    acc[i, 5, 4] <- tn/(tn+fp)
    acc[i, 5, 5] <- (tp+tn)/(tn+fp+tp+fn)

    ## full linear model bootstrapped:
    ## dcei <- compute_differential_causal_effects(
    ##     normal, dn,
    ##     tumor, dt,
    ##     method = "full",
    ##     bootstrap = TRUE, runs = bsruns, replace = 0, frac = 0.5
    ## )
    ## dceiflbs <- dcei

    ## dceib <- dcei$dce
    ## dceib[which(abs(dceib) > 0.5)] <- 1
    ## dceib[which(abs(dceib) <= 0.5)] <- 0

    ## dcei <- dcei$dce
    ## coridx <- which(dcet != 0 | dcei != 0)
    ## acc[i, 7, 1] <- cor(as.vector(dcet[coridx]), as.vector(dcei[coridx]), method = "p")
    ## acc[i, 7, 2] <- dist(rbind(as.vector(dcet), as.vector(dcei)))
    ## tp <- sum(dcetb == 1 & dceib == 1)
    ## fp <- sum(dcetb == 0 & dceib == 1)
    ## fn <- sum(dcetb == 1 & dceib == 0)
    ## tn <- sum(dcetb == 0 & dceib == 0)
    ## acc[i, 7, 3] <- tp/(tp+fn)
    ## acc[i, 7, 4] <- tn/(tn+fp)
    ## acc[i, 7, 5] <- (tp+tn)/(tn+fp+tp+fn)

    ## normal
    dcei <- compute_differential_causal_effects(
        normal, dn,
        tumor, dt
    )
    dcein <- dcei

    dceib <- dcei$dce
    dceib[which(abs(dceib) > 0.5)] <- 1
    dceib[which(abs(dceib) <= 0.5)] <- 0

    dcei <- dcei$dce
    coridx <- which(dcet != 0 | dcei != 0)
    acc[i, 1, 1] <- cor(as.vector(dcet[coridx]), as.vector(dcei[coridx]), method = "p")
    acc[i, 1, 2] <- dist(rbind(as.vector(dcet), as.vector(dcei)))
    tp <- sum(dcetb == 1 & dceib == 1)
    fp <- sum(dcetb == 0 & dceib == 1)
    fn <- sum(dcetb == 1 & dceib == 0)
    tn <- sum(dcetb == 0 & dceib == 0)
    acc[i, 1, 3] <- tp/(tp+fn)
    acc[i, 1, 4] <- tn/(tn+fp)
    acc[i, 1, 5] <- (tp+tn)/(tn+fp+tp+fn)

    ## random base line:
    dcei <- dcet
    dcei[which(gtc != 0)] <- runif(sum(gtc != 0), lB, uB)
    dcer <- list(dce = dcei, graph = as(gtc, "graphNEL"), dcefull = dcei)
    dceir <- dcer
    class(dceir) <- "dce"

    dcerb <- dcer$dce
    dcerb[which(abs(dcerb) > 0.5)] <- 1
    dcerb[which(abs(dcerb) <= 0.5)] <- 0

    dcer <- dcer$dce
    coridx <- which(dcet != 0 | dcer != 0)
    acc[i, 2, 1] <- cor(as.vector(dcet[coridx]), as.vector(dcer[coridx]), method = "p")
    acc[i, 2, 2] <- dist(rbind(as.vector(dcet), as.vector(dcer)))
    tp <- sum(dcetb == 1 & dcerb == 1)
    fp <- sum(dcetb == 0 & dcerb == 1)
    fn <- sum(dcetb == 1 & dcerb == 0)
    tn <- sum(dcetb == 0 & dcerb == 0)
    acc[i, 2, 3] <- tp/(tp+fn)
    acc[i, 2, 4] <- tn/(tn+fp)
    acc[i, 2, 5] <- (tp+tn)/(tn+fp+tp+fn)

    ## ## simple correlation?:
    ## dcec <- cor(rbind(dn, dt))*gtc
    ## dcec <- list(dce = dcec, graph = as(gtc, "graphNEL"), dcefull = dcec)
    ## dceic <- dcec
    ## class(dceic) <- "dce"

    ## dcerb <- dcer$dce
    ## dcerb[which(abs(dcerb) > 0.5)] <- 1
    ## dcerb[which(abs(dcerb) <= 0.5)] <- 0

    ## dcer <- dcer$dce
    ## coridx <- which(dcet != 0 | dcer != 0)
    ## acc[i, 2, 1] <- cor(as.vector(dcet[coridx]), as.vector(dcer[coridx]), method = "p")
    ## acc[i, 2, 2] <- dist(rbind(as.vector(dcet), as.vector(dcer)))
    ## tp <- sum(dcetb == 1 & dcerb == 1)
    ## fp <- sum(dcetb == 0 & dcerb == 1)
    ## fn <- sum(dcetb == 1 & dcerb == 0)
    ## tn <- sum(dcetb == 0 & dcerb == 0)
    ## acc[i, 2, 3] <- tp/(tp+fn)
    ## acc[i, 2, 4] <- tn/(tn+fp)
    ## acc[i, 2, 5] <- (tp+tn)/(tn+fp+tp+fn)

                           #if (acc[i, 1, 1] > 0.9) { break() }
                                        #}

    ## acc[i, , ]

}

for (filen in 1:100) {
    if (!file.exists(paste("dce/dce", n, paste(m, collapse = "_"), sd, perturb, filen, ".rda", sep = "_"))) {
        break()
    }
}

save(acc, gtnfeat, file = paste("dce/dce", n, paste(m, collapse = "_"), sd, perturb, filen, ".rda", sep = "_"))

stop()

## euler commands (using local R version):

module load bioconductor/3.6

module load curl/7.49.1

module load gmp/5.1.3

##

ram=1000

rm error.txt

rm output.txt

rm .RData

queue=4

genes=100
perturb=0
runs=100
prob=0.01

## parameters: n, m[1], m[2], sd
bsub -M ${ram} -q normal.${queue}h -n 1 -e error.txt -o output.txt -R "rusage[mem=${ram}]" "R/bin/R --silent --no-save --args '${genes}' '1000' '100' '1' '${runs}' '${perturb}' '${prob}' < dce_sim.r"

for i in {2..10}; do
bsub -M ${ram} -q normal.${queue}h -n 1 -e error.txt -o output.txt -R "rusage[mem=${ram}]" "R/bin/R --silent --no-save --args '${genes}' '1000' '100' '1' '${runs}' '${perturb}' '${prob}' < dce_sim.r"
done

## results:

path <- "~/Mount/Euler/"

n <- 100
m <- c(1000, 100)
sd <- 1
perturb <- 0

## combine several into one matrix:

library(abind)
acc2 <- NULL
for (filen in 1:100) {
    if (file.exists(paste0(path, paste("dce/dce", n, paste(m, collapse = "_"), sd, perturb, filen, ".rda", sep = "_")))) {
        load(paste0(path, paste("dce/dce", n, paste(m, collapse = "_"), sd, perturb, filen, ".rda", sep = "_")))
        acc2 <- abind(acc2, acc, along = 1)
    }
}
acc <- acc2

## load(paste0(path, paste("dce/dce", n, paste(m, collapse = "_"), sd, ".rda", sep = "_")))

runs <- dim(acc)[1]
par(mfrow=c(2,3))
boxplot(acc[seq_len(runs), 1:7, 1], col = c(rgb(1,0,0), rgb(0.5,0.5,0.5), rgb(0,1,0), rgb(0,0,1)), main="Correlation")
boxplot(acc[seq_len(runs), 1:7, 2], col = c(rgb(1,0,0), rgb(0.5,0.5,0.5), rgb(0,1,0), rgb(0,0,1)), main="Distance")
boxplot(acc[seq_len(runs), 1:7, 3], col = c(rgb(1,0,0), rgb(0.5,0.5,0.5), rgb(0,1,0), rgb(0,0,1)), main="Sensitivity")
boxplot(acc[seq_len(runs), 1:7, 4], col = c(rgb(1,0,0), rgb(0.5,0.5,0.5), rgb(0,1,0), rgb(0,0,1)), main="Specificity")
boxplot(acc[seq_len(runs), 1:7, 5], col = c(rgb(1,0,0), rgb(0.5,0.5,0.5), rgb(0,1,0), rgb(0,0,1)), main="Accuracy")

## combine:

path <- "~/Mount/Euler/"

m <- c(1000, 100)
sd <- 1

show3 <- c(0,0.1,-0.1,2)
show2 <- c(1,5,2)
show <- c(10,50,100)
par(mfrow=c(length(show),length(show3)))
for (i in show) {
    for (j in show3) {
        acc2 <- NULL
        for (filen in 1:1000) {
            if (file.exists(paste0(path, paste("dce/dce", i, paste(m, collapse = "_"), sd, j, filen, ".rda", sep = "_")))) {
                load(paste0(path, paste("dce/dce", i, paste(m, collapse = "_"), sd, j, filen, ".rda", sep = "_")))
                acc2 <- abind(acc2, acc, along = 1)
            }
        }
        acc <- acc2
        runs <- dim(acc)[1]
        boxplot(acc[seq_len(runs), show2, 1], col = c(rgb(1,0,0), rgb(0.5,0.5,0.5), rgb(0,1,0), rgb(0,0,1)), main="Correlation", ylim = c(-1,1))
        abline(h=c(0.9,0.75,0.5,0.25), col = rgb(0,0,0,0.75), lty = 2)
    }
}

dev.print("temp.pdf", device = pdf)

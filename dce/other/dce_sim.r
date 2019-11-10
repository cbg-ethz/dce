source("dce/R/main.R")
source("dce/R/utils.R")

library(tidyverse)
library(purrr)
library(graph)
library(pcalg)
library(assertthat)
library(igraph)
library(matlib)

m <- numeric(2)
n <- as.numeric(commandArgs(TRUE)[1])
m[1] <- as.numeric(commandArgs(TRUE)[2])
m[2] <- as.numeric(commandArgs(TRUE)[3])
sd <- as.numeric(commandArgs(TRUE)[4])
runs <- as.numeric(commandArgs(TRUE)[5])
perturb <- as.numeric(commandArgs(TRUE)[6])
p <- as.numeric(commandArgs(TRUE)[7])
cormeth <- commandArgs(TRUE)[8]

## n <- 10; m <- c(100, 10); mu <- 0; sd <- 1; runs <- 100; perturb <- 0; cormeth <- "p"; p <- "rand"

if (is.na(runs)) {
    runs <- 100 # simulation runs
}
if (is.na(perturb)) {
    perturb <- 0
}
if (is.na(p)) {
    p <- "rand" # edge prob of the dag
}

print(p)

## uniform limits:
lB <- c(-1,0)
uB <- c(0,1)
## others:
## n <- 10 # number of nodes
## m <- c(100,100) # number of samples tumor and normal
## sd <- 0.1 # standard deviation for variable distributions
## the fraction of true pos (causal effects that are differential)
truepos <- 0.9 # if we sample -1 to 1 this is not necessary # aida samples only pos effects
mu <- 0

acc <- array(0, c(runs,5,1),
             dimnames = list(runs = paste0("run_", seq_len(runs)),
                             methods = c("dce", "random", "full linear", "simple correlation", "test"),
                             metrics = c("correlation")))
gtnfeat <- array(0, c(runs, 6),
                 dimnames = list(runs = paste0("run_", seq_len(runs)),
                                 features = c("avg children", "avg parents",
                                              "childless", "parentless",
                                              "maxpathlength", "density")))

for (i in 1:runs) {
                                        #for (j in 1:1000) {
                                        #   set.seed(j)
    if (p %in% "rand") {
        p2 <- runif(1, 0, 1)
    } else {
        p2 <- p
    }
    
    normal <- randomDAG_2(n, p2, lB, uB)
    tumor <- newWeights(normal, lB, uB, truepos) # resample edge weights

    dn <- rmvDAG_2(m[2], normal, normpars = c(mu,sd))
    dt <- rmvDAG_2(m[1], tumor, normpars = c(mu,sd))

    ## no confounding
    ## blanket <- as(normal, "matrix")
    ## blanket[upper.tri(blanket)] <- 0
    ## blanket[(1:(n-1)+c(0,(1:(n-2)*n))+n)] <- 1
    ## normalAdj <- as(normal, "matrix")*blanket
    ## normalVec <- normalAdj[which(normalAdj != 0)]
    ## normal <- as(blanket, "graphNEL")
    ## for (j in 1:(n-1)) {
    ##     normal@edgeData@data[[j]]$weight <- normalVec[j]
    ## }
    ## tumorAdj <- as(tumor, "matrix")*blanket
    ## tumorVec <- tumorAdj[which(tumorAdj != 0)]
    ## tumor <- as(blanket, "graphNEL")
    ## for (j in 1:(n-1)) {
    ##     tumor@edgeData@data[[j]]$weight <- tumorVec[j]
    ## }

    ## for all confounding set p <- 1

    gm <- as(normal, "matrix")
    gm[which(gm != 0)] <- 1

    cn <- trueEffects(normal)
    ct <- trueEffects(tumor)

    gtc <- mnem:::mytc(gm) # transitively closed graph as matrix

    ## save features of gtn, which might correlate with accuracy:
    gtr <- nem::transitive.reduction(gm)
    gtnfeat[i, 1] <- mean(apply(gm, 1, sum))
    gtnfeat[i, 2] <- mean(apply(gm, 2, sum))
    gtnfeat[i, 3] <- sum(apply(gm, 1, sum) == 0)
    gtnfeat[i, 4] <- sum(apply(gm, 2, sum) == 0)
    gtnfeat[i, 5] <- max(igraph::distance_table(graph_from_adjacency_matrix(gtr))$res)
    gtnfeat[i, 6] <- sum(gm)

    ## ground truth:
    dcet <- (cn - ct)*gtc # gtn for differential causal effects
    dcegtn <- list(dce = dcet, graph = normal, dcefull = dcet)
    class(dcegtn) <- "dce"

    ## perturb network:

    if (perturb < 0) {
        adjn <- as(normal, "matrix")
        remedge <- sample(which(adjn != 0), floor(sum(adjn != 0)*abs(perturb)))
        adjn[remedge] <- 0
    }
    if (perturb > 0) {
        adjn <- as(normal, "matrix")
        addedge <- sample(which(adjn == 0 & upper.tri(adjn)), floor(sum(adjn != 0)*perturb))
        adjn[addedge] <- 1
    }

    ## our inference:

    ## test method:
    dcei <- fulllin(
        normal, dn,
        tumor, dt, conf = 0
    )
    dceitest <- dcei

    dcei <- dcei$dce
    acc[i, 5, 1] <- cor(as.vector(dcet), as.vector(dcei), method = cormeth, use = "complete.obs")

    ## full linear model:
    dcei <- fulllin(
        normal, dn,
        tumor, dt
    )
    dceifl <- dcei

    dcei <- dcei$dce
    acc[i, 3, 1] <- cor(as.vector(dcet), as.vector(dcei), method = cormeth, use = "complete.obs")

    ## normal
    if (!any(m < n)) [
           dcei <- compute_differential_causal_effects(
               normal, dn,
               tumor, dt
           )
           dcein <- dcei
           
           dcei <- dcei$dce
           acc[i, 1, 1] <- cor(as.vector(dcet), as.vector(dcei), method = cormeth, use = "complete.obs")
    }

    ## simple correlation:
    dcei <- (cor(dn) - cor(dt))*gtc
    dcec <- list(dce = dcei, graph = as(gtc, "graphNEL"), dcefull = dcei,
                 cen = cor(dt)*gtc, cet = cor(dt)*gtc)
    dceic <- dcec
    class(dceic) <- "dce"

    dcec <- dcec$dce
    acc[i, 4, 1] <- cor(as.vector(dcet), as.vector(dcec), method = cormeth, use = "complete.obs")

    ## random base line:
    dcei <- dceicn <- dceict <- dcet
    dcei[which(gtc != 0)] <- runif(sum(gtc != 0), lB[1], uB[2])
    dceicn[which(gtc != 0)] <- runif(sum(gtc != 0), lB[1], uB[2])
    dceict[which(gtc != 0)] <- runif(sum(gtc != 0), lB[1], uB[2])
    dcer <- list(dce = dcei, graph = as(gtc, "graphNEL"), dcefull = dcei,
                 cen = dceicn, cet = dceict)
    dceir <- dcer
    class(dceir) <- "dce"

    dcer <- dcer$dce
    coridx <- which(dcet != 0 | dcer != 0)
    acc[i, 2, 1] <- cor(as.vector(dcet), as.vector(dcer), method = cormeth, use = "complete.obs")

    ## acc[i, , ]

    cat(paste0(i, "."))
}

for (filen in 1:100) {
    if (!file.exists(paste("dce/dce", n, paste(m, collapse = "_"), sd, perturb, p, filen, ".rda", sep = "_"))) {
        break()
    }
}

save(acc, gtnfeat, file = paste("dce/dce", n, paste(m, collapse = "_"), sd, perturb, p, filen, ".rda", sep = "_"))

stop()

## euler commands (using local R version):

module add /cluster/apps/modules/modulefiles/new

module load python/3.7.1

module load bioconductor/3.6

module load curl/7.49.1

module load gmp/5.1.3

module load star/2.5.3a

module load samtools/1.2

##

system("scp dce/other/dce_sim.r euler.ethz.ch:dce_sim.r")
system("scp dce/R/main.r euler.ethz.ch:dce/R/main.r")
system("scp dce/R/utils.r euler.ethz.ch:dce/R/utils.r")

##

ram=1000

rm error.txt

rm output.txt

rm .RData

queue=4

genes=100
perturb=0
runs=1
prob=rand
cormeth=p
tumor=100
normal=100

bsub -M ${ram} -q normal.${queue}h -n 1 -e error.txt -o output.txt -R "rusage[mem=${ram}]" "R/bin/R --silent --no-save --args '${genes}' '${tumor}' '${normal}' '1' '${runs}' '${perturb}' '${prob}' '${cormeth}' < dce_sim.r"

for i in {2..100}; do
bsub -M ${ram} -q normal.${queue}h -n 1 -e error.txt -o output.txt -R "rusage[mem=${ram}]" "R/bin/R --silent --no-save --args '${genes}' '${tumor}' '${normal}' '1' '${runs}' '${perturb}' '${prob}' '${cormeth}' < dce_sim.r"
done

## results:

path <- "~/Mount/Euler/"

n <- 100
m <- c(100,100)
sd <- 1
perturb <- 0
p <- "rand" # rand for random

## combine several into one matrix:

library(abind)
acc2 <- gtnfeat2 <- NULL
for (filen in 1:100) {
    if (file.exists(paste0(path, paste("dce/dce", n, paste(m, collapse = "_"), sd, perturb, p, filen, ".rda", sep = "_")))) {
        load(paste0(path, paste("dce/dce", n, paste(m, collapse = "_"), sd, perturb, p, filen, ".rda", sep = "_")))
        acc2 <- abind(acc2, acc, along = 1)
        gtnfeat2 <- abind(gtnfeat2, gtnfeat, along = 1)
        cat(paste0(filen, "."))
    }
}
acc <- acc2
gtnfeat <- gtnfeat2

## load(paste0(path, paste("dce/dce", n, paste(m, collapse = "_"), sd, ".rda", sep = "_")))

## differential causal effects plus gtn features

col <- 2:6
show <- 1:5
runs <- dim(acc)[1]
par(mfrow=c(1,2))
boxplot(acc[seq_len(runs), show, 1], col = col, main="Correlation")
show <- 1:6
boxplot(gtnfeat[seq_len(runs), show], col = 1:6, main="Ground truth Features", log = "y")

source("https://raw.githubusercontent.com/cbg-ethz/mnem/master/R/mnems_low.r")
col <- 2:6
show <- 1:5
runs <- dim(acc)[1]
par(mfrow=c(1,2))
myboxplot(acc[seq_len(runs), show, 1], col = col, main="Correlation")
show <- 1:6
myboxplot(gtnfeat[seq_len(runs), show], col = 1:6, main="Ground truth Features", log = "y")

## correlated with network features:
print(cor(acc[,,1], gtnfeat[,2]))

## conversion between data.frame and array
acc.df <- as.data.frame.table(acc)

acc.arr <- xtabs(Freq ~ runs + methods + metrics, acc.df)

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

## figures:

Ga <- c("A=B", "B=C", "C=D", "A=E", "E=F", "A=F", "A=D")

Ew <- round(runif(length(Ga), -1, 1), 2)
Ew2 <- round(runif(length(Ga), -1, 1), 2)
Dw <- Ew-Ew2

edgecol <- rgb(abs(Dw)/max(abs(Dw)), 0, 0)

pdf("temp.pdf", width = 15, height = 5)
par(mfrow=c(1,3))
plotDnf(Ga, edgelabel = Ew, main = "Causal effects under condition A")
plotDnf(Ga, edgelabel = Ew2, main = "Causal effects under condition B")
plotDnf(Ga, edgelabel = Dw, edgecol = edgecol, main = "Differential causal effects")
dev.off()


source("dce/R/main.r")
source("dce/R/utils.r")

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
cormeth <- commandArgs(TRUE)[8]
dmeth <- commandArgs(TRUE)[9]

## n <- 10; m <- c(1000, 100); mu <- 0; sd <- 1; runs <- 100; perturb <- 0; dmeth <- "manhattan"; cormeth <- "p"; p <- 1

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
lB <- c(-10,0)
uB <- c(0,10)
## others:
#n <- 10 # number of nodes
#m <- c(100,100) # number of samples tumor and normal
#sd <- 0.1 # standard deviation for variable distributions
## the fraction of true pos (causal effects that are differential)
truepos <- 0.9 # if we sample -1 to 1 this is not necessary # aida samples only pos effects
bsruns <- 100
mu <- 0

trueEffects <- function(g) {
    wm <- t(wgtMatrix(g))
    te <- wm
    wmExp <- wm
    while (any(wmExp != 0)) {
        wmExp <- wmExp%*%wm
        te[which(wmExp != 0)] <- wmExp[which(wmExp != 0)]
    }
    te[which(te == 0)] <- NA
    return(te)
}

acc <- array(0, c(runs,5,6),
             dimnames = list(runs = paste0("run_", seq_len(runs)),
                             methods = c("dce", "random", "full linear", "simple correlation", "test"),
                             metrics = c("correlation", "distance",
                                         "causal effects normal cor.", "causal effects tumor cor.",
                                         "causal effects normal dist.", "causal effects tumor dist.")))
gtnfeat <- array(0, c(runs, 6, 2),
                 dimnames = list(runs = paste0("run_", seq_len(runs)),
                                 features = c("avg children", "avg parents",
                                              "childless", "parentless",
                                              "maxpathlength", "density"),
                                 gtn = c("originial", "transitive closure")))
dcefeat <- array(0, c(runs, 6, 2),
                 dimnames = list(runs = paste0("run_", seq_len(runs)),
                                 features = c("avg children", "avg parents",
                                              "childless", "parentless",
                                              "maxpathlength", "density"),
                                 gtn = c("originial", "transitive closure")))

for (i in 1:runs) {
    cat(i)
                                        #for (j in 1:1000) {
                                        #   set.seed(j)
    normal <- randomDAG_2(n, p, lB, uB)
    tumor <- newWeights(normal, lB, uB, truepos) # resample edge weights

    dn <- rmvDAG_2(m[2], normal, normpars = c(mu,sd))
    dt <- rmvDAG_2(m[1], tumor, normpars = c(mu,sd))

    gm <- as(normal, "matrix")
    gm[which(gm != 0)] <- 1

    ## no confounding
    ## blanket <- gm
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

    ## for all confounding set p <- 1 (compute_differential_causal_effects gives singularity error)

    cn <- trueEffects(normal)
    ct <- trueEffects(tumor)

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

    accfun <- function(a,b) {
        a <- as.vector(a)
        b <- as.vector(b)
        z <- dist(rbind(a,b), method = dmeth)
        return(z)
    }

    ## test method::
    dcei <- fulllin(
        normal, dn,
        tumor, dt, conf = 0
    )
    dceitest <- dcei

    dcei <- dcei$dce
    coridx <- which(dcet != 0 | dcei != 0)
    acc[i, 5, 1] <- cor(as.vector(dcet[coridx]), as.vector(dcei[coridx]), method = cormeth)
    acc[i, 5, 2] <- accfun(dcet[coridx], dcei[coridx])
    dcei <- fulllin(
        normal, dn,
        tumor, dt, conf = 0, diff = 0
    )
    coridxn <- which(cn*gtc != 0 | dcei$cen != 0)
    acc[i, 5, 3] <- cor(as.vector((cn*gtc)[coridxn]), as.vector(dcei$cen[coridxn]), method = cormeth)
    acc[i, 5, 5] <- accfun((cn*gtc)[coridxn], dcei$cen[coridxn])
    coridxt <- which(ct*gtc != 0 | dcei$cet != 0)
    acc[i, 5, 4] <- cor(as.vector((ct*gtc)[coridxt]), as.vector(dcei$cet[coridxt]), method = cormeth)
    acc[i, 5, 6] <- accfun((ct*gtc)[coridxt], dcei$cet[coridxt])

    ## full linear model:
    dcei <- fulllin(
        normal, dn,
        tumor, dt
    )
    dceifl <- dcei

    dcei <- dcei$dce
    coridx <- which(dcet != 0 | dcei != 0)
    acc[i, 3, 1] <- cor(as.vector(dcet[coridx]), as.vector(dcei[coridx]), method = cormeth)
    acc[i, 3, 2] <- accfun(dcet[coridx], dcei[coridx])
    dcei <- fulllin(
        normal, dn,
        tumor, dt, diff = 0
    )
    coridxn <- which(cn*gtc != 0 | dcei$cen != 0)
    acc[i, 3, 3] <- cor(as.vector((cn*gtc)[coridxn]), as.vector(dcei$cen[coridxn]), method = cormeth)
    acc[i, 3, 5] <- accfun((cn*gtc)[coridxn], dcei$cen[coridxn])
    coridxt <- which(ct*gtc != 0 | dcei$cet != 0)
    acc[i, 3, 4] <- cor(as.vector((ct*gtc)[coridxt]), as.vector(dcei$cet[coridxt]), method = cormeth)
    acc[i, 3, 6] <- accfun((ct*gtc)[coridxt], dcei$cet[coridxt])

    ## normal
    dcei <- compute_differential_causal_effects(
        normal, dn,
        tumor, dt
    )
    dcein <- dcei

    dcei <- dcei$dce
    coridx <- which(dcet != 0 | dcei != 0)
    acc[i, 1, 1] <- cor(as.vector(dcet[coridx]), as.vector(dcei[coridx]), method = cormeth)
    acc[i, 1, 2] <- accfun(dcet[coridx], dcei[coridx])
    dcei <- dcein
    coridxn <- which(cn*gtc != 0 | dcei$cen != 0)
    acc[i, 1, 3] <- cor(as.vector((cn*gtc)[coridxn]), as.vector(dcei$cen[coridxn]), method = cormeth)
    acc[i, 1, 5] <- accfun((cn*gtc)[coridxn], dcei$cen[coridxn])
    coridxt <- which(ct*gtc != 0 | dcei$cet != 0)
    acc[i, 1, 4] <- cor(as.vector((ct*gtc)[coridxt]), as.vector(dcei$cet[coridxt]), method = cormeth)
    acc[i, 1, 6] <- accfun((ct*gtc)[coridxt], dcei$cet[coridxt])

    ## simple correlation:
    dcei <- (cor(dn) - cor(dt))*gtc
    dcec <- list(dce = dcei, graph = as(gtc, "graphNEL"), dcefull = dcei,
                 cen = cor(dt)*gtc, cet = cor(dt)*gtc)
    dceic <- dcec
    class(dceic) <- "dce"

    dcec <- dcec$dce
    coridx <- which(dcet != 0 | dcec != 0)
    acc[i, 4, 1] <- cor(as.vector(dcet[coridx]), as.vector(dcec[coridx]), method = cormeth)
    acc[i, 4, 2] <- accfun(dcet[coridx], dcei[coridx])
    dcei <- dceic
    coridxn <- which(cn*gtc != 0 | (cor(dn)*gtc) != 0)
    acc[i, 4, 3] <- cor(as.vector((cn*gtc)[coridxn]), as.vector((cor(dn)*gtc)[coridxn]), method = cormeth)
    acc[i, 4, 5] <- accfun((cn*gtc)[coridxn], dcei$cen[coridxn])
    coridxt <- which(ct*gtc != 0 | (cor(dt)*gtc) != 0)
    acc[i, 4, 4] <- cor(as.vector((ct*gtc)[coridxt]), as.vector((cor(dt)*gtc)[coridxt]), method = cormeth)
    acc[i, 4, 6] <- accfun((ct*gtc)[coridxt], (cor(dn)*gtc)[coridxt])

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
    acc[i, 2, 1] <- cor(as.vector(dcet[coridx]), as.vector(dcer[coridx]), method = cormeth)
    acc[i, 2, 2] <- accfun(dcet[coridx], dcei[coridx])
    dcei <- dceir
    coridxn <- which(cn*gtc != 0 | dcei$cen != 0)
    acc[i, 2, 3] <- cor(as.vector((cn*gtc)[coridxn]), as.vector(dcei$cen[coridxn]), method = cormeth)
    acc[i, 2, 5] <- accfun((cn*gtc)[coridxn], dcei$cen[coridxn])
    coridxt <- which(ct*gtc != 0 | dcei$cet != 0)
    acc[i, 2, 4] <- cor(as.vector((ct*gtc)[coridxt]), as.vector(dcei$cet[coridxt]), method = cormeth)
    acc[i, 2, 6] <- accfun((ct*gtc)[coridxt], dcei$cet[coridxt])

    ## ## full linear model bootstrapped:
    ## dcei <- compute_differential_causal_effects(
    ##     normal, dn,
    ##     tumor, dt,
    ##     method = "full",
    ##     bootstrap = TRUE, runs = bsruns, replace = 0, frac = 0.5
    ## )
    ## dceiflbs <- dcei

    ## dcei <- dcei$dce
    ## coridx <- which(dcet != 0 | dcei != 0)
    ## acc[i, 7, 1] <- cor(as.vector(dcet[coridx]), as.vector(dcei[coridx]), method = cormeth)
    ## acc[i, 7, 2] <- dist(rbind(as.vector(dcet), as.vector(dcei)))

    ## if (acc[i, 1, 1] > 0.9) { break() }

    ## if (acc[i, 5, 2] < acc[i, 3, 2]) { stop() }

    ## acc[i, , ]

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

genes=50
perturb=0
runs=100
prob=0.05
cormeth=p
dmeth=euclidean

bsub -M ${ram} -q normal.${queue}h -n 1 -e error.txt -o output.txt -R "rusage[mem=${ram}]" "R/bin/R --silent --no-save --args '${genes}' '1000' '100' '1' '${runs}' '${perturb}' '${prob}' '${cormeth}' '${dmeth}' < dce_sim.r"

for i in {2..100}; do
bsub -M ${ram} -q normal.${queue}h -n 1 -e error.txt -o output.txt -R "rusage[mem=${ram}]" "R/bin/R --silent --no-save --args '${genes}' '1000' '100' '1' '${runs}' '${perturb}' '${prob}' '${cormeth}' '${dmeth}' < dce_sim.r"
done

## results:

path <- "~/Mount/Euler/"

n <- 50
m <- c(1000, 100)
sd <- 1
perturb <- 0
p <- 0.2

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

## differential causal effects:
show <- 1:5
runs <- dim(acc)[1]
par(mfrow=c(1,2))
boxplot(acc[seq_len(runs), show, 1], col = c(rgb(1,0,0), rgb(0.5,0.5,0.5), rgb(0,1,0), rgb(0,0,1)), main="Correlation")

## gtn features:
show <- 1:6
runs <- dim(acc)[1]
boxplot(gtnfeat[seq_len(runs), show, 1], col = c(rgb(1,0,0), rgb(0.5,0.5,0.5), rgb(0,1,0), rgb(0,0,1)), main="Ground truth Features")

## causal effects:
show <- 1:5
runs <- dim(acc)[1]
par(mfrow=c(1,4))
boxplot(acc[seq_len(runs), show, 3], col = c(rgb(1,0,0), rgb(0.5,0.5,0.5), rgb(0,1,0), rgb(0,0,1)), main="Normal Effects Correlation")
boxplot(acc[seq_len(runs), show, 4], col = c(rgb(1,0,0), rgb(0.5,0.5,0.5), rgb(0,1,0), rgb(0,0,1)), main="Tumor Effects Correlation")
boxplot(acc[seq_len(runs), show, 5], col = c(rgb(1,0,0), rgb(0.5,0.5,0.5), rgb(0,1,0), rgb(0,0,1)), main="Normal Effects Distance")
boxplot(acc[seq_len(runs), show, 6], col = c(rgb(1,0,0), rgb(0.5,0.5,0.5), rgb(0,1,0), rgb(0,0,1)), main="Tumor Effects Distance")

## correlated with network features:
print(cor(acc[,,2], gtnfeat[,,1]))

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

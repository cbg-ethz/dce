source("dce/R/main.R")
source("dce/R/utils.R")
source("dce/R/simulations.R")

library(tidyverse)
library(purrr)
library(graph)
library(pcalg)
library(assertthat)
library(igraph)
library(matlib)
library(nem)
library(mnem)
library(expm)

m <- numeric(2)
n <- as.numeric(commandArgs(TRUE)[1])
m[1] <- as.numeric(commandArgs(TRUE)[2])
m[2] <- as.numeric(commandArgs(TRUE)[3])
sd <- as.numeric(commandArgs(TRUE)[4])
simruns <- as.numeric(commandArgs(TRUE)[5])
perturb <- as.numeric(commandArgs(TRUE)[6])
p <- as.numeric(commandArgs(TRUE)[7])
cormeth <- commandArgs(TRUE)[8]
bootstrap <- commandArgs(TRUE)[9]

## n <- 10; m <- c(20, 20); mu <- 0; sd <- 1; simruns <- 1; perturb <- 0; cormeth <- "p"; p <- "runif"; bootstrap <- "1"

if (is.na(simruns)) {
    simruns <- 100 # simulation runs
}
if (is.na(perturb)) {
    perturb <- 0
}
if (is.na(p)) {
    p <- "runif" # edge prob of the dag
}
if (is.na(bootstrap)) {
    bootstrap <- "none"
} else if (bootstrap %in% "1") {
    bootstrap <- c("basic", "full")
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

## Rprof(line.profiling = TRUE)
simRes <- simDce(n,m,simruns,mu,sd,c(lB,uB),truepos,perturb,cormeth,p,bootstrap=bootstrap,TRUE) # maybe additional bootstrap args: frac,replace,runs
## summaryRprof()$sampling.time; head(summaryRprof(lines = "show")$by.self)

for (filen in 1:100) {
    if (!file.exists(paste("dce/dce", n, paste(m, collapse = "_"), sd, perturb, p, filen, ".rda", sep = "_"))) {
        break()
    }
}

save(simRes, file = paste("dce/dce", n, paste(m, collapse = "_"), sd, perturb, p, filen, ".rda", sep = "_"))

stop("success")

plot.dceSim(simRes, col = col, border = col, showMeth = c(2,5,3,6,1,4), showFeat = c(1,2,3))

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
system("scp dce/R/main.r euler.ethz.ch:dce/R/main.R")
system("scp dce/R/utils.r euler.ethz.ch:dce/R/utils.R")
system("scp dce/R/simulations.r euler.ethz.ch:dce/R/simulations.R")

##

ram=5000

rm error.txt

rm output.txt

rm .RData

queue=4

genes=50 # 10, 50, 100
perturb=-0.5 # 0, 0.5, -0.5
runs=1
prob=runif
cormeth=p
tumor=100 # depends on genes... 10*genes, 0.5*genes
normal=100 # see above? 2*genes, 0.25*genes
bootstrap=1

bsub -M ${ram} -q normal.${queue}h -n 1 -e error.txt -o output.txt -R "rusage[mem=${ram}]" "R/bin/R --silent --no-save --args '${genes}' '${tumor}' '${normal}' '1' '${runs}' '${perturb}' '${prob}' '${cormeth}' '${bootstrap}' < dce_sim.r"

for i in {1..10}; do
bsub -M ${ram} -q normal.${queue}h -n 1 -e error.txt -o output.txt -R "rusage[mem=${ram}]" "R/bin/R --silent --no-save --args '${genes}' '${tumor}' '${normal}' '1' '${runs}' '${perturb}' '${prob}' '${cormeth}' '${bootstrap}' < dce_sim.r"
done

## results:

path <- "~/Mount/Euler/"

n <- 50
m <- c(100,100)
sd <- 1
perturb <- 0.5
p <- "runif" # rand for random

## combine several into one matrix:

library(abind)
simRes2 <- NULL
for (filen in 1:100) {
    if (file.exists(paste0(path, paste("dce/dce", n, paste(m, collapse = "_"), sd, perturb, p, filen, ".rda", sep = "_")))) {
        load(paste0(path, paste("dce/dce", n, paste(m, collapse = "_"), sd, perturb, p, filen, ".rda", sep = "_")))
        simRes2$acc <- abind(simRes2$acc, simRes$acc, along = 1)
        simRes2$gtnFeat <- abind(simRes2$gtnFeat, simRes$gtnFeat, along = 1)
        cat(paste0(filen, "."))
    }
}
simRes <- simRes2
class(simRes) <- "dceSim"

## load(paste0(path, paste("dce/dce", n, paste(m, collapse = "_"), sd, ".rda", sep = "_")))

## differential causal effects plus gtn features

source("dce/R/main.R")
col <- rgb(c(0.5,1,0,0.1),c(0.5,0,0,0.1),c(0.5,0,1,0.1),0.75)
plot.dceSim(simRes, col = col, border = col, showMeth = c(2,5,3,6,1,4), showFeat = c(1,2,3), methNames = c("basic", "basic\n bootstrap", "full", "full\n bootstrap", "random", "correlation"))

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

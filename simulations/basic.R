# assumes that I am executed from main directory
source("simulations/utils.R")
devtools::load_all("./dce")

## parameters for simulations:

runs <- 100
p <- 0.2 # edge prob
## uniform limits:
lB <- -1
uB <- 1
##
n <- 10 # nodes
m <- 100 # samples
sd <- 1 # standard deviation

## run simulations

acc <- matrix(0, runs, 4)

for (i in 1:runs) {
  normal <- randomDAG(n, p, lB, uB)
  dn <- rmvDAG_2(m, normal, normpars = c(0,sd))

  tumor <- newWeights(normal) # resample edge weights
  dt <- rmvDAG_2(m, tumor, normpars = c(0,sd))

  cn <- trueCov(normal)
  ct <- trueCov(tumor)

  gm <- as(normal, "matrix")
  gm[which(gm != 0)] <- 1

  gtc <- mnem:::mytc(gm) # transitively closed graph as matrix
  diag(gtc) <- 0

  dcet <- (cn - ct)*gtc # gtn for differential causal effects

  dcei <- dce::compute_differential_causal_effects(
      normal, dn,
      tumor, dt
  )
  dcei <- t(dcei) # ?
  # dcei <- dcet*0 # inferred
  # for (j in 1:n) {
  #   dcei[j, ] <- idaFast(j, 1:n, cov(dn), normal) - idaFast(j, 1:n, cov(dt), tumor)
  # }

  dcei <- dcei*gtc

  ## what as accuracy? Account for random chance.

  acc[i, 1] <- cor(as.vector(dcet), as.vector(dcei), method = "s")
  acc[i, 2] <- dist(matrix(c(dcet, dcei), 2))
}

## plot results

## comparative plot:

efreq <- round(dcet[which(dcet != 0)], 2)
efreqi <- round(dcei[which(dcei != 0)], 2)
efreq[abs(efreq) > 2] <- 2 # even though the gtn does not allow for dce > 2, noise can pass that limit
efreqi[abs(efreqi) > 2] <- 2

par(mfrow=c(1,2))
dnf <- mnem:::adj2dnf(gtc)
dnf <- dnf[grep("=", dnf)] # this preprocessing will be implemented in plotDnf (bug)
mnem::plotDnf(dnf, labels = efreq, edgecol = rgb(abs(efreq)/2,0,(2-abs(efreq))/2))
mnem::plotDnf(dnf, labels = efreqi, edgecol = rgb(abs(efreqi)/2,0,(2-abs(efreqi))/2))

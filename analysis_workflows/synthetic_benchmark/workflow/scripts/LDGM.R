LDGM <- function(wt.X,mt.X) {
    wt.LC <- sin((pi/2)*cor(wt.X,method="k"))
    diag(wt.LC) <- 1
    mt.LC <- sin((pi/2)*cor(mt.X,method="k"))
    diag(mt.LC) <- 1
    wt.f <- paste0("wt.LC.temp.",runif(1),".txt")
    mt.f <- paste0("mt.LC.temp.",runif(1),".txt")
    write.table(wt.LC,file=wt.f,row.names=FALSE,col.names=FALSE,sep=",")
    write.table(mt.LC,file=mt.f,row.names=FALSE,col.names=FALSE,sep=",")
    theta.f <- paste0("Theta.",runif(1),".csv")
    system(paste0("matlab -batch \"Sigma1 = readmatrix('",mt.f,"'); Sigma2 = readmatrix('",wt.f,"'); lambda = 0.26924; cd('LDGM'); Theta = differential_graph(Sigma1,Sigma2,lambda); writematrix(Theta,'",theta.f,"'); exit\""))
    Theta <- as(read.csv(paste0("LDGM/",theta.f),header=FALSE),"matrix")
    system(paste0("rm ", wt.f, " ", mt.f, " LDGM/", theta.f))
    return(Theta)
}
LDGM.perm <- function(x, y, iter = 1000) {
    z <- LDGM(x,y)
    p <- z * 0
    for (i in seq_len(iter)) {
        xp <- x[, sample(seq_len(ncol(x)), ncol(x))]
        yp <- y[, sample(seq_len(ncol(y)), ncol(y))]
        zp <- LDGM(xp,yp)
        idx <- which(abs(zp) >= abs(z))
        p[idx] <- p[idx] + 1
    }
    p <- p / iter
    return(p)
}

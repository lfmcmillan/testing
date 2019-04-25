library(clustglm)
# source("/am/miro/home/mcmilllo/Marsden/BinaryNegCorr/clustglm_LMedit.R")

construct_row_membership <- function(N,pi_r) {
    R <- length(pi_r)
    return(sample(1:R, N, prob=pi_R))
}

construct_col_membership <- function(M,kappa_c) {
    C <- length(kappa_c)
    if (C == 1) col_membership <- rep(1,times=M)
    else {
        col_membership <- vector()
        for (cc in 1:(C-1)) {
            col_membership <- c(col_membership,rep(cc,round(M*kappa_c[cc])))
        }
        remaining <- M - length(col_membership)
        col_membership <- c(col_membership, rep(C,remaining))
    }
    col_membership
}

construct_dat <- function(N,M,theta,row_membership,col_membership) {
    dat_rows <- lapply(1:N,function(i) {
        dat_cols <- sapply(1:M,function(j) rbinom(1,1,theta[row_membership[i],col_membership[j]]))
    })
    dat <- do.call(rbind,dat_rows)
}

construct_longdat <- function(dat) {
    longdat <- mat2df(dat)
    longdat$ntrials <- rep(1,nrow(longdat))
    longdat$nsucc <- longdat$Y
    longdat$nfail <- 1-longdat$Y
    longdat
}

check_row_results <- function(output) {
    row_assignments <- apply(output$pp.list$rowclust,1,which.max)
    row_percent_correct <- sum(row_assignments==row_membership)/N*100
    row_percent_confident_correct <- sum(row_assignments==row_membership & 
                                             (output$pp.list$rowclust[,1] > 0.8 | 
                                                  output$pp.list$rowclust[,1] < 0.2))/N*100
    
    list(row_percent_correct = row_percent_correct,
         row_percent_confident_correct = row_percent_confident_correct)
}

check_col_results <- function(output) {
    col_assignments <- apply(output$pp.list$colclust,1,which.max)
    col_percent_correct <- sum(col_assignments==col_membership)/N*100
    col_percent_confident_correct <- sum(col_assignments==col_membership & 
                                             (output$pp.list$colclust[,1] > 0.8 |
                                                  output$pp.list$colclust[,1] < 0.2))/N*100
    
    list(col_percent_correct = col_percent_correct,
         col_percent_confident_correct = col_percent_confident_correct)
}

## Row clustering --------------------------------------------------------------
if (F) {
    
    N <- 100
    M <- 100
    pi_r <- c(0.5,0.5)
    kappa_c <- 1
    theta <- matrix(c(0.9,0.1),nrow=length(pi_r),byrow=TRUE)
    row_membership <- construct_row_membership(N=N,pi_r=pi_r)
    col_membership <- construct_col_membership(M=M,kappa_c=kappa_c)
    dat <- construct_dat(N=N,M=M,theta=theta,row_membership=row_membership,
                         col_membership=col_membership)
    longdat <- construct_longdat(dat)
    row2clustonly.out <- clustglm(Y~rowclust, family="binomial", data=longdat,
                                  fact4clust = "ROW", nclus=2, clustfactnames = "rowclust",
                                  start.control = list(randstarts=2), verbose=2)
    
    summary(row2clustonly.out)
    round(row2clustonly.out$pp.list$rowclust,2)
    check_row_results(row2clustonly.out)
    
    ## Now flip a few of the columns
    dat2 <- dat
    dat2[,91:100] <- 1-dat2[,91:100]
    longdat2 <- mat2df(dat2)
    
    row2clustonlyflip.out <- clustglm(Y~rowclust, family="binomial", data=longdat2,
                                      fact4clust = "ROW", nclus=2, clustfactnames = "rowclust",
                                      start.control = list(randstarts=2), verbose=2)
    
    summary(row2clustonlyflip.out)
    round(row2clustonlyflip.out$pp.list$rowclust,2)
}

## Column clustering -----------------------------------------------------------
if (F) {
    
    N <- 100
    M <- 100
    pi_r <- 1
    kappa_c <- c(0.5,0.5)
    theta <- matrix(c(0.9,0.1),nrow=length(pi_r),byrow=TRUE)
    row_membership <- construct_row_membership(N=N,pi_r=pi_r)
    col_membership <- construct_col_membership(M=M,kappa_c=kappa_c)
    dat <- construct_dat(N=N,M=M,theta=theta,row_membership=row_membership,
                         col_membership=col_membership)
    longdat <- construct_longdat(dat)
    
    col2clust.out <- clustglm(Y~colclust, family="binomial", data=longdat,
                              fact4clust = "COL", nclus=2, clustfactnames = "colclust",
                              start.control = list(randstarts=2), verbose=2)
    
    summary(col2clust.out)
    round(col2clust.out$pp.list$colclust,2)
    check_col_results(col2clust.out)
    
    ## Now flip a few of the columns
    dat2 <- dat
    dat2[,91:100] <- 1-dat2[,91:100]
    longdat2 <- mat2df(dat2)
    
    col2clustonlyflip.out <- clustglm(Y~colclust, family="binomial", data=longdat2,
                                      fact4clust = "COL", nclus=2, clustfactnames = "colclust",
                                      start.control = list(randstarts=2), verbose=2)
    
    summary(col2clustonlyflip.out)
    round(col2clustonlyflip.out$pp.list$colclust,2)
}

## Biclustering ----------------------------------------------------------------
if (F) {
    
    N <- 100
    M <- 100
    pi_r <- c(0.75,0.25)
    kappa_c <- c(0.5,0.5)
    theta <- matrix(c(0.9,0.1,0.3,0.7),nrow=length(pi_r),byrow=TRUE)
    row_membership <- construct_row_membership(N=N,pi_r=pi_r)
    col_membership <- construct_col_membership(M=M,kappa_c=kappa_c)
    
    dat <- construct_dat(N=N,M=M,theta=theta,row_membership=row_membership, 
                         col_membership=col_membership)
    longdat <- construct_longdat(dat)
    
    clustonly.out <- clustglm(Y~rowclust*colclust, family="binomial", data=longdat,
                              fact4clust = c("ROW","COL"), nclus=c(2,2), 
                              clustfactnames = c("rowclust","colclust"),
                              start.control = list(randstarts=10),
                              EM.control = list(startEMcycles=10), verbose=1)
    
    summary(clustonly.out)
    round(clustonly.out$pp.list$rowclust,2)
    round(clustonly.out$pp.list$colclust,2)
    check_row_results(clustonly.out)
    check_col_results(clustonly.out)
    
    ## Now flip a few of the columns
    dat2 <- dat
    dat2[,91:100] <- 1-dat2[,91:100]
    longdat2 <- construct_longdat(dat2)
    
    clustonlyflip.out <- clustglm(Y~rowclust+colclust, family="binomial", data=longdat2,
                                  fact4clust = c("ROW","COL"), nclus=c(2,2), 
                                  clustfactnames = c("rowclust","colclust"),
                                  start.control = list(randstarts=10),
                                  EM.control = list(startEMcycles=10), verbose=1)
    
    summary(clustonlyflip.out)
    round(clustonlyflip.out$pp.list$rowclust,2)
    round(clustonlyflip.out$pp.list$colclust,2)
    check_row_results(clustonlyflip.out)
    check_col_results(clustonlyflip.out)
    
}


if (F) {
    N <- 100
    M <- 100
    pi_r <- c(0.5,0.5)
    kappa_c <- 1
    theta <- matrix(c(0.45,0.55),nrow=length(pi_r),byrow=TRUE)
    row_membership <- construct_row_membership(N=N,pi_r=pi_r)
    col_membership <- construct_col_membership(M=M,kappa_c=kappa_c)
    dat <- construct_dat(N=N,M=M,theta=theta,row_membership=row_membership,
                         col_membership=col_membership)
    
    clust2.out <- cluster_binary_rows(Y=dat, R=2)
    round(clust2.out$z_hat)
    
    N <- 100
    M <- 100
    pi_r <- c(0.3,0.3,0.4)
    kappa_c <- 1
    theta <- matrix(c(0.9,0.5,0.1),nrow=length(pi_r),byrow=TRUE)
    row_membership <- construct_row_membership(N=N,pi_r=pi_r)
    col_membership <- construct_col_membership(M=M,kappa_c=kappa_c)
    dat <- construct_dat(N=N,M=M,theta=theta,row_membership=row_membership,
                         col_membership=col_membership)
    
    clust2.out <- cluster_binary_rows(Y=dat, R=2)
    round(clust2.out$z_hat,2)

    clust3.out <- cluster_binary_rows(Y=dat, R=3)
    round(clust3.out$z_hat,2)
}

## temporary change
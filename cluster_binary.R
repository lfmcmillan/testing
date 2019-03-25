cluster_binary_rows <- function(Y, R, maxiter=100, tol=1e-4) {

    N <- nrow(Y)
    M <- ncol(Y)
    
    pi_hat <- rDirichlet(R,rep(1,times=R))
    theta_hat <- matrix(runif(R*M),nrow=R)
    
    iter <- 1
    diff <- 1
    while (iter <= maxiter & diff > tol) {
        print(paste("iter",iter))
        prev_pi_hat <- pi_hat
        prev_theta_hat <- theta_hat
        
        z_hat_raw <- sapply(1:R, function(r) {
            sapply(1:N, function(i) {
                pi_hat[r]*prod(theta_hat[r,]^Y[i,]*(1-theta_hat[r,])^(1-Y[i,]))
            })
        })
        z_hat <- z_hat_raw/rowSums(z_hat_raw)
        
        print("z_hat")
        print(z_hat)

        theta_hat <- t(z_hat)%*%Y/colSums(z_hat)
        print("theta_hat")
        print(theta_hat)
        
        # Note that optim minimizes so we make the log-likelihood negative to 
        # maximize it
        loglikefun <- function(pi_r) {
            # Move any overly small pi_r values away from 0
            pi_r[which(pi_r < 1e-8)] <- 1e-8
            
            -sum(z_hat%*%log(pi_r/sum(pi_r)))
        } 
        
        pi_optim <- tryCatch({optim(pi_hat,loglikefun)},
                             error=function(err) {
                                 message("Error in fitting pi_hat")
                                 err })
        pi_hat <- pi_optim$par/sum(pi_optim$par)
        print("pi_hat")
        print(pi_hat)
        print(paste("convergence",pi_optim$convergence))
        
        diff <- max(c(max(abs(theta_hat - prev_theta_hat)), pi_hat - prev_pi_hat))
        print(paste("diff",diff))
        
        iter <- iter+1
    }
    
    list(pi_hat=pi_hat, theta_hat=theta_hat, z_hat=z_hat, iter=iter-1)
}

rDirichlet <- function(ndraw, alphvec){
    ## See Wikipedia on the Dirichlet distribution for confirmation:
    gamdraw <- rgamma(ndraw, shape=alphvec, rate=1)
    gamdraw / sum(gamdraw)
}
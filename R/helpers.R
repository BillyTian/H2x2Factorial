

main.cluster <- function(m_bar, CV, power, delta_x, rho, pi_x, pi_z, correction, sigma2_y, a, max_n, z_a, z_b){

  omega_2 <- sigma2_y*(1+(m_bar-1)*rho)/(pi_x*(1-pi_x)*m_bar) *
    ( (1+(m_bar-1)*rho)^2/((1+(m_bar-1)*rho)^2-CV^2*m_bar*rho*(1-rho)) + pi_z*(1-rho)*(1+(m_bar-1)*rho)^2/((1-pi_z)*((1+(m_bar-2)*rho)*(1+(m_bar-1)*rho)^2 + CV^2*m_bar*rho^2*(1-rho))) )

  #Check the omega to be positive
  if (omega_2<=0)
    stop("Variance inflation factor is in an abnormal value. Please check whether the provided CV parameter is unrealistically large")

  if (correction==FALSE){
    n <- (z_a+z_b)^2*omega_2/delta_x^2
    n.final <- ceiling(n)

  } else if (correction==TRUE){
    n <- 2
    try.power <- 0
    while ((try.power < power) & (n < max_n)){
      n <- n+1
      #Noncentral t (two-sided)
      try.power <- pt(qt(1-a/2, n-2), n-2, ncp=delta_x/sqrt(omega_2/n), lower.tail = F)
      + pt(qt(a/2, n-2), n-2, ncp=delta_x/sqrt(omega_2/n))

      if (n==max_n) {
        warning("It achieves the maximum number of cluster. The current prediction does not guarantee the required power")
      }
    }
    n.final <- n
  }
  return(n.final)
}


marginal.cluster <- function(m_bar, CV, power, delta_x, rho, pi_x, correction, sigma2_y, a, max_n, z_a, z_b){

  eta2 <- m_bar*(1-rho)/(1+(m_bar-1)*rho)*(1-CV^2*m_bar*rho*(1-rho)/((1+(m_bar-1)*rho)^2))-m_bar
  omega_x <- sigma2_y*(1-rho)/((m_bar+eta2)*pi_x*(1-pi_x))

  #Check the omega to be positive
  if (omega_x<=0)
    stop("Variance inflation factor is in an abnormal value. Please check whether the provided CV parameter is unrealistically large")

  if (correction==FALSE){
    n <- (z_a+z_b)^2*omega_x/delta_x^2
    n.final <- ceiling(n)

  } else if (correction==TRUE){
    n <- 2
    try.power <- 0
    while ((try.power < power) & (n < max_n)){
      n <- n+1
      #Noncentral t (two-sided)
      try.power <- pt(qt(1-a/2, n-2), n-2, ncp=delta_x/sqrt(omega_x/n), lower.tail = F)
      + pt(qt(a/2, n-2), n-2, ncp=delta_x/sqrt(omega_x/n))

      if (n==max_n) {
        warning("It achieves the maximum number of cluster. The current prediction does not guarantee the required power")
      }
    }
    n.final <- n
  }
  return(n.final)
}


main.ind <- function(m_bar, CV, power, delta_z, rho, pi_x, pi_z, sigma2_y, z_a, z_b){
  omega_3 <- sigma2_y*(1-rho)*(1+(m_bar-1)*rho)^3/((1-pi_z)*pi_x*(1-pi_x)*m_bar*((1+(m_bar-2)*rho)*(1+(m_bar-1)*rho)^2 + CV^2*m_bar*rho^2*(1-rho)))

  #Check the omega to be positive
  if (omega_3<=0)
    stop("Variance inflation factor is in an abnormal value. Please check whether the provided CV parameter is unrealistically large")

  n <- (z_a+z_b)^2*omega_3/delta_z^2
  n.final <- ceiling(n)
  return(n.final)
}


marginal.ind <- function(m_bar, CV, power, delta_z, rho, pi_z, sigma2_y, z_a, z_b){
  eta1 <- -m_bar*rho/(1+(m_bar-1)*rho)*(1-CV^2*m_bar*rho*(1-rho)/((1+(m_bar-1)*rho)^2))
  omega_z <- sigma2_y*(1-rho)/((m_bar+eta1)*pi_z*(1-pi_z))

  #Check the omega to be positive
  if (omega_z<=0)
    stop("Variance inflation factor is in an abnormal value. Please check whether the provided CV parameter is unrealistically large")

  n <- (z_a+z_b)^2*omega_z/delta_z^2
  n.final <- ceiling(n)
  return(n.final)
}


interaction <- function(m_bar, CV, power, delta_xz, rho, pi_x, pi_z, sigma2_y, z_a, z_b){
  eta1 <- -m_bar*rho/(1+(m_bar-1)*rho)*(1-CV^2*m_bar*rho*(1-rho)/((1+(m_bar-1)*rho)^2))
  omega_xz <- sigma2_y*(1-rho)/((m_bar+eta1)*pi_x*(1-pi_x)*pi_z*(1-pi_z))

  #Check the omega to be positive
  if (omega_xz<=0)
    stop("Variance inflation factor is in an abnormal value. Please check whether the provided CV parameter is unrealistically large")

  n <- (z_a+z_b)^2*omega_xz/delta_xz^2
  n.final <- ceiling(n)
  return(n.final)
}


####### Simultaneous Testing #######

main.joint <- function(m_bar, CV, power, delta_x, delta_z, rho, pi_x, pi_z, correction, sigma2_y, a, max_n){

  omega_2 <- sigma2_y*(1+(m_bar-1)*rho)/(pi_x*(1-pi_x)*m_bar) *
    ( (1+(m_bar-1)*rho)^2/((1+(m_bar-1)*rho)^2-CV^2*m_bar*rho*(1-rho)) + pi_z*(1-rho)*(1+(m_bar-1)*rho)^2/((1-pi_z)*((1+(m_bar-2)*rho)*(1+(m_bar-1)*rho)^2 + CV^2*m_bar*rho^2*(1-rho))) )
  omega_3 <- sigma2_y*(1-rho)*(1+(m_bar-1)*rho)^3/((1-pi_z)*pi_x*(1-pi_x)*m_bar*((1+(m_bar-2)*rho)*(1+(m_bar-1)*rho)^2 + CV^2*m_bar*rho^2*(1-rho)))
  omega_23 <- pi_z*omega_3

  #Check the omega to be positive
  if (omega_2<=0 || omega_3<=0)
    stop("Variance inflation factor is in an abnormal value. Please check whether the provided CV parameter is unrealistically large")

  if (correction==FALSE){
    n <- 1
    try.power <- 0
    while ((try.power < power) & (n < max_n)){
      n <- n+1
      #non-centrality parameter (include the unknown of interest, n)
      beta23 <- c(delta_x, delta_z)
      omega <- matrix(c(omega_2, omega_23, omega_23, omega_3), nrow=2, byrow=T)

      theta <- n*t(beta23) %*% solve(omega) %*% beta23
      try.power <- pchisq(qchisq(1-a, 2), 2, ncp = theta, lower.tail = F)

      if (n==max_n) {
        warning("It achieves the maximum number of cluster. The current prediction does not guarantee the required power")
      }
    }
    n.final <- n
  } else if (correction==TRUE){
    n <- 2
    try.power <- 0
    while ((try.power < power) & (n < max_n)){
      n <- n+1
      #non-centrality parameter (include the unknown of interest, n)
      beta23 <- c(delta_x, delta_z)
      omega <- matrix(c(omega_2, omega_23, omega_23, omega_3), nrow=2, byrow=T)

      theta <- n*t(beta23) %*% solve(omega) %*% beta23
      try.power <- pf(qf(1-a, 2, n-2), 2, n-2, ncp = theta, lower.tail = F)

      if (n==max_n) {
        warning("It achieves the maximum number of cluster. The current prediction does not guarantee the required power")
      }
    }
    n.final <- n
  }
  return(n.final)
}



marginal.joint <- function(m_bar, CV, power, delta_x, delta_z, rho, pi_x, pi_z, correction, sigma2_y, a, max_n, seed_mix, size_mix){
  eta1 <- -m_bar*rho/(1+(m_bar-1)*rho)*(1-CV^2*m_bar*rho*(1-rho)/((1+(m_bar-1)*rho)^2))
  eta2 <- m_bar*(1-rho)/(1+(m_bar-1)*rho)*(1-CV^2*m_bar*rho*(1-rho)/((1+(m_bar-1)*rho)^2))-m_bar
  omega_x <- sigma2_y*(1-rho)/((m_bar+eta2)*pi_x*(1-pi_x))
  omega_z <- sigma2_y*(1-rho)/((m_bar+eta1)*pi_z*(1-pi_z))

  #Check the omega to be positive
  if (omega_x<=0 || omega_z<=0)
    stop("Variance inflation factor is in an abnormal value. Please check whether the provided CV parameter is unrealistically large")

  if (correction==FALSE){
    n <- 1
    try.power <- 0
    while ((try.power < power) & (n < max_n)){
      n <- n+1
      theta <- n*(delta_x^2/omega_x + delta_z^2/omega_z)
      try.power <- pchisq(qchisq(1-a, 2), 2, ncp = theta, lower.tail = F)

      if (n==max_n) {
        warning("It achieves the maximum number of cluster. The current prediction does not guarantee the required power")
      }
    }
    n.final <- n
  } else if (correction==TRUE){
    n <- 2
    try.power <- 0
    while ((try.power < power) & (n < max_n)){
      set.seed(seed_mix)
      n <- n+1
      #Simulate the mixed distribution (CENTRAL) to identify rejection region bound
      f.distn <- rf(size_mix, 1, n-2)
      chisq.distn <- rchisq(size_mix, 1)
      mix.distn <- f.distn + chisq.distn
      crt.value <- quantile(mix.distn, 1-a)

      #Simulate the mixed distribution (NONCENTRAL VERSION) for power calculation
      nc.f.distn <- rf(size_mix, 1, n-2, ncp = n*delta_x^2/omega_x)
      nc.chisq.distn <- rchisq(size_mix, 1, ncp = n*delta_z^2/omega_z)
      nc.mix.distn <- nc.f.distn + nc.chisq.distn

      try.power <- mean(nc.mix.distn>crt.value)

      if (n==max_n) {
        warning("It achieves the maximum number of cluster. The current prediction does not guarantee the required power")
      }
    }
    n.final <- n
  }
  return(n.final)
}



#library(mvtnorm)


main.IU <- function(m_bar, CV, power, delta_x, delta_z, rho, pi_x, pi_z, correction, sigma2_y, a, max_n){
  omega_2 <- sigma2_y*(1+(m_bar-1)*rho)/(pi_x*(1-pi_x)*m_bar) *
    ( (1+(m_bar-1)*rho)^2/((1+(m_bar-1)*rho)^2-CV^2*m_bar*rho*(1-rho)) + pi_z*(1-rho)*(1+(m_bar-1)*rho)^2/((1-pi_z)*((1+(m_bar-2)*rho)*(1+(m_bar-1)*rho)^2 + CV^2*m_bar*rho^2*(1-rho))) )
  omega_3 <- sigma2_y*(1-rho)*(1+(m_bar-1)*rho)^3/((1-pi_z)*pi_x*(1-pi_x)*m_bar*((1+(m_bar-2)*rho)*(1+(m_bar-1)*rho)^2 + CV^2*m_bar*rho^2*(1-rho)))
  omega_23 <- pi_z*omega_3

  #Check the omega to be positive
  if (omega_2<=0 || omega_3<=0)
    stop("Variance inflation factor is in an abnormal value. Please check whether the provided CV parameter is unrealistically large")

  if (correction==FALSE){
    n <- 1
    try.power <- 0
    while ((try.power < power) & (n < max_n)){
      n <- n+1
      mean_W <- c(delta_x/sqrt(omega_2/n), delta_z/sqrt(omega_3/n))
      cov_W <- (omega_23/n) / (sqrt(omega_2/n)*sqrt(omega_3/n))
      sigma_W <- matrix(c(1, cov_W, cov_W, 1), nrow=2, byrow=T)

      try.power <- pmvnorm(lower=rep(qnorm(1-a/2),2), upper=rep(Inf,2), mean=mean_W, sigma=sigma_W) +
        pmvnorm(lower=c(qnorm(1-a/2),-Inf), upper=c(Inf,qnorm(a/2)), mean=mean_W, sigma=sigma_W) +
        pmvnorm(lower=rep(-Inf,2), upper=rep(qnorm(a/2),2), mean=mean_W, sigma=sigma_W) +
        pmvnorm(lower=c(-Inf,qnorm(a/2)), upper=c(qnorm(a/2),Inf), mean=mean_W, sigma=sigma_W)

      if (n==max_n) {
        warning("It achieves the maximum number of cluster. The current prediction does not guarantee the required power")
      }
    }
    n.final <- n

  } else if (correction==TRUE){
    n <- 2
    try.power <- 0
    while ((try.power < power) & (n < max_n)){
      n <- n+1
      #non-centrality parameter (include the unknown of interest, n)
      mean_W <- c(delta_x/sqrt(omega_2/n), delta_z/sqrt(omega_3/n))
      cov_W <- (omega_23/n) / (sqrt(omega_2/n)*sqrt(omega_3/n))
      sigma_W <- matrix(c(1, cov_W, cov_W, 1), nrow=2, byrow=T)

      try.power <- pmvt(df=n-2, lower=c(qt(1-a/2, n-2), qt(1-a/2, n-2)), upper=rep(Inf,2), delta=mean_W, sigma=sigma_W) +
        pmvt(df=n-2, lower=c(qt(1-a/2, n-2),-Inf), upper=c(Inf,qt(a/2, n-2)), delta=mean_W, sigma=sigma_W) +
        pmvt(df=n-2, lower=rep(-Inf,2), upper=c(qt(a/2, n-2), qt(a/2, n-2)), delta=mean_W, sigma=sigma_W) +
        pmvt(df=n-2, lower=c(-Inf,qt(a/2, n-2)), upper=c(qt(a/2, n-2),Inf), delta=mean_W, sigma=sigma_W)

      if (n==max_n) {
        warning("It achieves the maximum number of cluster. The current prediction does not guarantee the required power")
      }
    }
    n.final <- n
  }
  return(n.final)
}


marginal.IU <- function(m_bar, CV, power, delta_x, delta_z, rho, pi_x, pi_z, correction, sigma2_y, a, max_n){
  eta1 <- -m_bar*rho/(1+(m_bar-1)*rho)*(1-CV^2*m_bar*rho*(1-rho)/((1+(m_bar-1)*rho)^2))
  eta2 <- m_bar*(1-rho)/(1+(m_bar-1)*rho)*(1-CV^2*m_bar*rho*(1-rho)/((1+(m_bar-1)*rho)^2))-m_bar
  omega_x <- sigma2_y*(1-rho)/((m_bar+eta2)*pi_x*(1-pi_x))
  omega_z <- sigma2_y*(1-rho)/((m_bar+eta1)*pi_z*(1-pi_z))

  #Check the omega to be positive
  if (omega_x<=0 || omega_z<=0)
    stop("Variance inflation factor is in an abnormal value. Please check whether the provided CV parameter is unrealistically large")

  if (correction==FALSE){
    n <- 1
    try.power <- 0
    while ((try.power < power) & (n < max_n)){
      n <- n+1
      wmean.c <- sqrt(n)*delta_x/sqrt(omega_x)
      wmean.i <- sqrt(n)*delta_z/sqrt(omega_z)
      try.power <- pnorm(qnorm(1-a/2), mean = wmean.c, lower.tail = F)*pnorm(qnorm(1-a/2), mean = wmean.i, lower.tail = F)
      + pnorm(qnorm(1-a/2), mean = wmean.c, lower.tail = F)*pnorm(qnorm(a/2), mean = wmean.i)
      + pnorm(qnorm(a/2), mean = wmean.c)*pnorm(qnorm(1-a/2), mean = wmean.i, lower.tail = F)
      + pnorm(qnorm(a/2), mean = wmean.c)*pnorm(qnorm(a/2), mean = wmean.i)

      if (n==max_n) {
        warning("It achieves the maximum number of cluster. The current prediction does not guarantee the required power")
      }
    }
    n.final <- n

  } else if (correction==TRUE){
    n <- 2
    try.power <- 0
    while ((try.power < power) & (n < max_n)){
      n <- n+1
      c.ncp <- sqrt(n)*delta_x/sqrt(omega_x)
      i.mean <- sqrt(n)*delta_z/sqrt(omega_z)
      try.power <- pt(qt(1-a/2, n-2), n-2, c.ncp, lower.tail = F)*pnorm(qnorm(1-a/2), mean = i.mean, lower.tail = F)
      + pt(qt(1-a/2, n-2), n-2, c.ncp, lower.tail = F)*pnorm(qnorm(a/2), mean = i.mean)
      + pt(qt(a/2, n-2), n-2, c.ncp)*pnorm(qnorm(1-a/2), mean = i.mean, lower.tail = F)
      + pt(qt(a/2, n-2), n-2, c.ncp)*pnorm(qnorm(a/2), mean = i.mean)

      if (n==max_n) {
        warning("It achieves the maximum number of cluster. The current prediction does not guarantee the required power")
      }
    }
    n.final <- n
  }
  return(n.final)
}


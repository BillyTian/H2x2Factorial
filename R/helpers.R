marginal.cluster <- function(m_bar, CV, power, delta_x, rho, pi_x, correction, sigma2_y, a, max_n, z_a, z_b){

  #test
  #m_bar=30
  #CV=2.2
  #power=0.8
  #delta_x=0.25
  #rho=0.1
  #pi_x=0.5
  #correction=F
  #sigma2_y=1
  #a=0.05
  #max_n=1e8
  #z_a=qnorm(1-0.05/2)
  #z_b=qnorm(0.8)

  eta2 <- m_bar*(1-rho)/(1+(m_bar-1)*rho)*(1-CV^2*m_bar*rho*(1-rho)/((1+(m_bar-1)*rho)^2))-m_bar
  omega_x <- sigma2_y*(1-rho)/((m_bar+eta2)*pi_x*(1-pi_x))

  #Check the omega to be positive
  if (omega_x<=0)
    stop("Variance inflation factor is in an abnormal value. Please check whether the provided CV parameter is unrealistically large")

  if (correction==FALSE){
    n <- (z_a+z_b)^2*omega_x/delta_x^2
    n.final <- ceiling(n)

  } else if (correction==TRUE){
    n <- 3
    try.power <- 0
    while ((try.power < power) & (n < max_n)){
      #Noncentral t (two-sided)
      try.power <- pt(qt(1-a/2, n-2), n-2, ncp=delta_x/sqrt(omega_x/n), lower.tail = F)
      + pt(qt(a/2, n-2), n-2, ncp=delta_x/sqrt(omega_x/n))
      n <- n+1
      if (n==max_n) {
        warning("It achieves the maximum number of cluster. The current prediction does not guarantee the required power")
      }
    }
    n.final <- n
  }
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


joint <- function(m_bar, CV, power, delta_x, delta_z, rho, pi_x, pi_z, correction, sigma2_y, a, max_n, seed_mix, size_mix){
  eta1 <- -m_bar*rho/(1+(m_bar-1)*rho)*(1-CV^2*m_bar*rho*(1-rho)/((1+(m_bar-1)*rho)^2))
  eta2 <- m_bar*(1-rho)/(1+(m_bar-1)*rho)*(1-CV^2*m_bar*rho*(1-rho)/((1+(m_bar-1)*rho)^2))-m_bar
  omega_x <- sigma2_y*(1-rho)/((m_bar+eta2)*pi_x*(1-pi_x))
  omega_z <- sigma2_y*(1-rho)/((m_bar+eta1)*pi_z*(1-pi_z))

  #Check the omega to be positive
  if (omega_x<=0 || omega_z<=0)
    stop("Variance inflation factor is in an abnormal value. Please check whether the provided CV parameter is unrealistically large")

  if (correction==FALSE){
    n <- 2
    try.power <- 0
    while ((try.power < power) & (n < max_n)){
      theta <- n*(delta_x^2/omega_x + delta_z^2/omega_z)
      try.power <- pchisq(qchisq(1-a, 2), 2, ncp = theta, lower.tail = F)
      n <- n+1
      if (n==max_n) {
        warning("It achieves the maximum number of cluster. The current prediction does not guarantee the required power")
      }
    }
    n.final <- n
  } else if (correction==TRUE){
    n <- 3
    try.power <- 0
    while ((try.power < power) & (n < max_n)){
      set.seed(seed_mix)
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
      n <- n+1
    }
  }
  n.final <- n
  return(n.final)
}


IU <- function(m_bar, CV, power, delta_x, delta_z, rho, pi_x, pi_z, correction, sigma2_y, a, max_n){
  eta1 <- -m_bar*rho/(1+(m_bar-1)*rho)*(1-CV^2*m_bar*rho*(1-rho)/((1+(m_bar-1)*rho)^2))
  eta2 <- m_bar*(1-rho)/(1+(m_bar-1)*rho)*(1-CV^2*m_bar*rho*(1-rho)/((1+(m_bar-1)*rho)^2))-m_bar
  omega_x <- sigma2_y*(1-rho)/((m_bar+eta2)*pi_x*(1-pi_x))
  omega_z <- sigma2_y*(1-rho)/((m_bar+eta1)*pi_z*(1-pi_z))

  #Check the omega to be positive
  if (omega_x<=0 || omega_z<=0)
    stop("Variance inflation factor is in an abnormal value. Please check whether the provided CV parameter is unrealistically large")

  if (correction==FALSE){
    n <- 2
    try.power <- 0
    while ((try.power < power) & (n < max_n)){
      wmean.c <- sqrt(n)*delta_x/sqrt(omega_x)
      wmean.i <- sqrt(n)*delta_z/sqrt(omega_z)
      try.power <- pnorm(qnorm(1-a/2), mean = wmean.c, lower.tail = F)*pnorm(qnorm(1-a/2), mean = wmean.i, lower.tail = F)
      + pnorm(qnorm(1-a/2), mean = wmean.c, lower.tail = F)*pnorm(qnorm(a/2), mean = wmean.i)
      + pnorm(qnorm(a/2), mean = wmean.c)*pnorm(qnorm(1-a/2), mean = wmean.i, lower.tail = F)
      + pnorm(qnorm(a/2), mean = wmean.c)*pnorm(qnorm(a/2), mean = wmean.i)
      n <- n+1
      if (n==max_n) {
        warning("It achieves the maximum number of cluster. The current prediction does not guarantee the required power")
      }
    }
    n.final <- n

  } else if (correction==TRUE){
    n <- 3
    try.power <- 0
    while ((try.power < power) & (n < max_n)){
      c.ncp <- sqrt(n)*delta_x/sqrt(omega_x)
      i.mean <- sqrt(n)*delta_z/sqrt(omega_z)
      try.power <- pt(qt(1-a/2, n-2), n-2, c.ncp, lower.tail = F)*pnorm(qnorm(1-a/2), mean = i.mean, lower.tail = F)
      + pt(qt(1-a/2, n-2), n-2, c.ncp, lower.tail = F)*pnorm(qnorm(a/2), mean = i.mean)
      + pt(qt(a/2, n-2), n-2, c.ncp)*pnorm(qnorm(1-a/2), mean = i.mean, lower.tail = F)
      + pt(qt(a/2, n-2), n-2, c.ncp)*pnorm(qnorm(a/2), mean = i.mean)
      n <- n+1
      if (n==max_n) {
        warning("It achieves the maximum number of cluster. The current prediction does not guarantee the required power")
      }
    }
    n.final <- n
  }
}

#' Calculate required number of clusters under the joint test
#' 
#' @param m_bar a mean cluster size.
#' @param CV a coefficient of variation of the cluster sizes.
#' @param power a power level requirement.
#' @param delta_x an (unstandardized) effect size of the marginal effect of cluster-level treatment.
#' @param delta_z an (unstandardized) effect size of the marginal effect of individual-level treatment.
#' @param rho an intraclass correlation coefficient.
#' @param pi_x a proportion of clusters randomized to the cluster-level treatment.
#' @param pi_z a proportion of individuals randomized to the individual-level treatment within each cluster.
#' @param correction a logical argument indicating whether a finite sample correction should be used.
#' @param sigma2_y total variance of outcome.
#' @param a type I error rate.
#' @param max_n an optional setting of maximum number of clusters.
#' @param seed_mix an optional setting of seed for conducting the simulation-based testing under a mixed distribution.
#' @param size_mix a size for the mixed distribution in the simulation-based procedure.
#' 
#' @return a number representing the required number of cluster.
#' @export
#' 
joint <- function(m_bar, CV, power, delta_x, delta_z, rho, pi_x, pi_z, correction, sigma2_y, a, max_n, seed_mix, size_mix){
  eta1 <- -m_bar*rho/(1+(m_bar-1)*rho)*(1-CV^2*m_bar*rho*(1-rho)/((1+(m_bar-1)*rho)^2))
  eta2 <- m_bar*(1-rho)/(1+(m_bar-1)*rho)*(1-CV^2*m_bar*rho*(1-rho)/((1+(m_bar-1)*rho)^2))-m_bar
  omega_x <- sigma2_y*(1-rho)/((m_bar+eta2)*pi_x*(1-pi_x))
  omega_z <- sigma2_y*(1-rho)/((m_bar+eta1)*pi_z*(1-pi_z))
  
  if (correction==F){
    n <- 2
    try.power <- 0
    while ((try.power < power) & (n < max_n)){
      theta <- n*(delta_x^2/omega_x + delta_z^2/omega_z)
      try.power <- pchisq(qchisq(1-a, 2), 2, ncp = theta, lower.tail = F)
      n <- n+1
      if (n==max_n) {
        warning("It achieves the maximum number of cluster. The current prediction does not guarantee the required power.")
      }
    }
    n.final <- n
  } else if (correction==T){
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
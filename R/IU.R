#' Calculate required number of clusters under the I-U test
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
#' 
#' @return a number representing the required number of cluster.
#' @export
#' 
IU <- function(m_bar, CV, power, delta_x, delta_z, rho, pi_x, pi_z, correction, sigma2_y, a, max_n){
  eta1 <- -m_bar*rho/(1+(m_bar-1)*rho)*(1-CV^2*m_bar*rho*(1-rho)/((1+(m_bar-1)*rho)^2))
  eta2 <- m_bar*(1-rho)/(1+(m_bar-1)*rho)*(1-CV^2*m_bar*rho*(1-rho)/((1+(m_bar-1)*rho)^2))-m_bar
  omega_x <- sigma2_y*(1-rho)/((m_bar+eta2)*pi_x*(1-pi_x))
  omega_z <- sigma2_y*(1-rho)/((m_bar+eta1)*pi_z*(1-pi_z))
  
  if (correction==F){
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
        warning("It achieves the maximum number of cluster. The current prediction does not guarantee the required power.")
      }
    }
    n.final <- n
    
  } else if (correction==T){
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
        warning("It achieves the maximum number of cluster. The current prediction does not guarantee the required power.")
      }
    }
    n.final <- n
  }
}
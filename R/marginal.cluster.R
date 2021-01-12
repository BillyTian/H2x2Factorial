#' Calculate required number of clusters under the cluster-level marginal test
#' 
#' @param m_bar a mean cluster size.
#' @param CV a coefficient of variation of the cluster sizes.
#' @param power a power level requirement.
#' @param delta_x an (unstandardized) effect size of the marginal effect of cluster-level treatment.
#' @param rho an intraclass correlation coefficient.
#' @param pi_x a proportion of clusters randomized to the cluster-level treatment.
#' @param correction a logical argument indicating whether a finite sample correction should be used.
#' @param sigma2_y total variance of outcome.
#' @param a type I error rate.
#' @param max_n an optional setting of maximum number of clusters.
#' @param z_a z_a needed in the closed-form sample size formula.
#' @param z_b z_b needed in the closed-form sample size formula.
#' 
#' @return a number representing the required number of cluster.
#' @export

marginal.cluster <- function(m_bar, CV, power, delta_x, rho, pi_x, correction, sigma2_y, a, max_n, z_a, z_b){
  eta2 <- m_bar*(1-rho)/(1+(m_bar-1)*rho)*(1-CV^2*m_bar*rho*(1-rho)/((1+(m_bar-1)*rho)^2))-m_bar
  omega_x <- sigma2_y*(1-rho)/((m_bar+eta2)*pi_x*(1-pi_x))
  
  if (correction==F){
    n <- (z_a+z_b)^2*omega_x/delta_x^2
    n.final <- ceiling(n)
    
  } else if (correction==T){
    n <- 3
    try.power <- 0
    while ((try.power < power) & (n < max_n)){
      #Noncentral t (two-sided)
      try.power <- pt(qt(1-a/2, n-2), n-2, ncp=delta_x/sqrt(omega_x/n), lower.tail = F) 
      + pt(qt(a/2, n-2), n-2, ncp=delta_x/sqrt(omega_x/n))
      n <- n+1
      if (n==max_n) {
        warning("It achieves the maximum number of cluster. The current prediction does not guarantee the required power.")
      }
    }
    n.final <- n
  }
  return(n.final)
}
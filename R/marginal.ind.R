#' Calculate required number of clusters under the individual-level marginal test
#' 
#' @param m_bar a mean cluster size.
#' @param CV a coefficient of variation of the cluster sizes.
#' @param power a power level requirement.
#' @param delta_z an (unstandardized) effect size of the marginal effect of individual-level treatment.
#' @param rho an intraclass correlation coefficient.
#' @param pi_z a proportion of individuals randomized to the individual-level treatment within each cluster.
#' @param sigma2_y total variance of outcome.
#' @param z_a z_a needed in the closed-form sample size formula.
#' @param z_b z_b needed in the closed-form sample size formula.
#' 
#' @return a number representing the required number of cluster.
#' @export
#' 
marginal.ind <- function(m_bar, CV, power, delta_z, rho, pi_z, sigma2_y, z_a, z_b){
  eta1 <- -m_bar*rho/(1+(m_bar-1)*rho)*(1-CV^2*m_bar*rho*(1-rho)/((1+(m_bar-1)*rho)^2))
  omega_z <- sigma2_y*(1-rho)/((m_bar+eta1)*pi_z*(1-pi_z))
  
  n <- (z_a+z_b)^2*omega_z/delta_z^2
  n.final <- ceiling(n)
  return(n.final)
}

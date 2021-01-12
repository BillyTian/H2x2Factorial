#' Calculate required number of clusters under the interaction test
#' 
#' @param m_bar a mean cluster size.
#' @param CV a coefficient of variation of the cluster sizes.
#' @param power a power level requirement.
#' @param delta_xz an (unstandardized) effect size of the treatment interaction.
#' @param rho an intraclass correlation coefficient.
#' @param pi_x a proportion of clusters randomized to the cluster-level treatment.
#' @param pi_z a proportion of individuals randomized to the individual-level treatment within each cluster.
#' @param sigma2_y total variance of outcome.
#' @param z_a z_a needed in the closed-form sample size formula.
#' @param z_b z_b needed in the closed-form sample size formula.
#' 
#' @return a number representing the required number of cluster.
#' @export
#' 
interaction <- function(m_bar, CV, power, delta_xz, rho, pi_x, pi_z, sigma2_y, z_a, z_b){
  eta1 <- -m_bar*rho/(1+(m_bar-1)*rho)*(1-CV^2*m_bar*rho*(1-rho)/((1+(m_bar-1)*rho)^2))
  omega_xz <- sigma2_y*(1-rho)/((m_bar+eta1)*pi_x*(1-pi_x)*pi_z*(1-pi_z))
  
  n <- (z_a+z_b)^2*omega_xz/delta_xz^2
  n.final <- ceiling(n)
  return(n.final)
}
#' @description Plot sample size estimations (m_bar-n combinations) under varying CV for a chosen test
#' 
#' @param m_lower a lower bound for the mean cluster sizes on the horizontal axis. Default is 10.
#' @param m_upper an upper bound for the mean cluster sizes on the horizontal axis. Default is 100.
#' @param m_step a step size for the horizontal axis for plotting the sample size combinations. Default is 2.
#' @param CV a vector of the coefficient of variation of the cluster sizes of interest. Default is c(0, 0.3, 0.6, 0.9).
#' @param color_palette a vector to specify the color choices corresponding to the lines in the plot. 
#' Default is c("#0F2080", "#85C0F9", "#DDCC77", "#F5793A", "#A95AA1").
#' @param line_width a vector to specify the widths of the lines in the plot. Default is rep(3, 5).
#' @param line_type a vector to specify the line types of the lines in the plot. Default is seq(1, 5, 1).
#' @param power a target power level. Default is 0.8.
#' @param alpha a type I error rate. Default is 0.05.
#' @param pi_x a proportion of clusters randomized to the cluster-level treatment. Default is 0.5, representing a balanced allocation.
#' @param pi_z a proportion of individuals randomized to the individual-level treatment within each cluster. Default is 0.5, representing a balanced allocation.
#' @param delta_x an effect size of the marginal effect of cluster-level treatment. Default is 0.25.
#' @param delta_z an effect size of the marginal effect of individual-level treatment. Default is 0.33.
#' @param delta_xz an effect size of the treatment interaction. Default is 0.3.
#' @param sigma2_y a total variance of the (continuous) outcome. Default is 1.
#' @param rho an intraclass correlation coefficient. Default is 0.01.
#' @param test a character argument indicating the type of hypothesis of interest. Supported choices include 
#' cluster (test of the marginal effect of cluster-level treatment), individual (test of the marginal effect of individual-level treatment),  
#' interaction (test of interaction), joint (joint test of the two marginal effects), I-U (Intersection-Union test of the two marginal effects).
#' Default is "cluster".  
#' @param correction a logical argument indicating whether a small sample correction should be used. Default is T.
#' @param max_n an optional setting of maximum number of clusters, possibly needed under test="cluster", "joint", or "I-U". Default is 1e8.
#' @param seed_mix an optional setting of seed for conducting the simulation-based testing under a mixed distribution, possibly needed under test="joint". Default is NULL.
#' @param size_mix a pre-specified size for the mixed distribution in the simulation-based procedure, possibly needed under test="joint". Default is 1e4.
#' 
#' @return a plot comparing the sample size requirements under different CV.
#' 
#' @examples
#' 
#' plot.H2x2Factorial(power=0.9, test="cluster")
#' 
#' @export



#Based on the desired test and power, give a plot with m on the x-axis, n on the y-axis. 
#Multiple lines if a vector of CV is specified. (CV vector with the length less or equal to 5 is suggested)

#Also restrict in the y-direction by the range function
#A color-blind-friendly palette was set by default and this can be changed or edited by the user

plot.H2x2Factorial <- function(m_lower=10,
                               m_upper=100,
                               m_step=2,
                               CV=c(0, 0.3, 0.6, 0.9),
                               color_palette=c("#0F2080", "#85C0F9", "#DDCC77", "#F5793A", "#A95AA1"),
                               line_width=rep(3, 5),
                               line_type=seq(1, 5, 1),
                               power=0.8,
                               alpha=0.05, pi_x=0.5, pi_z=0.5,
                               delta_x=0.25, delta_z=0.33, delta_xz=0.3,
                               sigma2_y=1,
                               rho=0.01, 
                               test="cluster",
                               correction=T, 
                               max_n=1e8,
                               seed_mix=123456,
                               size_mix=1e4){

  a <- alpha
  b <- 1-power
  z_a <- qnorm(1-a/2)
  z_b <- qnorm(1-b)
  
  if (test=="cluster") {
    v <- Vectorize(marginal.cluster, c("m_bar"))
  } else if (test=="individual") {
    v <- Vectorize(marginal.ind, c("m_bar"))
  } else if (test=="interaction") {
    v <- Vectorize(interaction, c("m_bar"))
  } else if (test=="joint") {
    v <- Vectorize(joint, c("m_bar"))
  } else if (test=="I-U") {
    v <- Vectorize(IU, c("m_bar"))
  }
  
  Mset <- seq(m_lower, m_upper, by=m_step)
  
  line.CV <- NULL

  for (i in 1:length(CV)){
    if (test=="cluster"){
      line.CV[[i]] <- v(m_bar=Mset, CV=CV[i], power=power, correction=correction, 
                        delta_x=delta_x, rho=rho, sigma2_y=sigma2_y, pi_x=pi_x, max_n=max_n, z_a=z_a, z_b=z_b, a=a)
    } else if (test=="individual"){
      line.CV[[i]] <- v(m_bar=Mset, CV=CV[i], power=power, 
                        delta_z=delta_z, rho=rho, sigma2_y=sigma2_y, pi_z=pi_z, z_a=z_a, z_b=z_b)
    } else if (test=="interaction"){
      line.CV[[i]] <- v(m_bar=Mset, CV=CV[i], power=power,
                        delta_xz=delta_xz, rho=rho, sigma2_y=sigma2_y, pi_x=pi_x, pi_z=pi_z, z_a=z_a, z_b=z_b)
    } else if (test=="joint"){
      line.CV[[i]] <- v(m_bar=Mset, CV=CV[i], power=power, correction=correction, 
                        delta_x=delta_x, delta_z=delta_z, rho=rho, sigma2_y=sigma2_y, pi_x=pi_x, pi_z=pi_z, max_n=max_n, a=a, seed_mix=seed_mix, size_mix=size_mix)
    } else if (test=="I-U"){
      line.CV[[i]] <- v(m_bar=Mset, CV=CV[i], power=power, correction=correction, 
                        delta_x=delta_x, delta_z=delta_z, rho=rho, sigma2_y=sigma2_y, pi_x=pi_x, pi_z=pi_z, max_n=max_n, a=a)
    }
  }
  
  #Compute the parameters for ylim
  calc.range <- NULL
  for (j in 1:length(CV)){
    calc.range <- c(calc.range, line.CV[[j]])
  }
  y_lower <- range(calc.range)[1]
  y_upper <- range(calc.range)[2]
  
  if (test=="cluster"){
    if (correction==T){title <- "Test for the cluster-level marginal effect \n with finite-sample correction"}
    else if (correction==F){title <- "Test for the cluster-level marginal effect \n without finite-sample correction"}
    
  } else if (test=="individual"){
    title <- "Test for the individual-level marginal effect"
    
  } else if (test=="interaction"){
    title <- "Interaction test"
    
  } else if (test=="joint"){
    if (correction==T){title <- "Joint test with finite-sample correction"}
    else if (correction==F){title <- "Joint test without finite-sample correction"}
    
  } else if (test=="I-U"){
    if (correction==T){title <- "Intersection-union test \n with finite-sample correction"}
    else if (correction==F){title <- "Intersection-union test \n without finite-sample correction"}
    
  }
  
  plot(Mset, line.CV[[1]], ylim=c(y_lower, y_upper), las=1,
       xlab=expression(bar(m)),
       main=title,
       ylab="", cex.lab=1, cex.axis=1, cex.main=1, cex=1, 
       type="l", lwd=line_width[1], col=color_palette[1], lty=line_type[1])
  mtext(expression(n),side=2,las=1,line=3,cex=1.2)
  
  for (k in 2:length(CV)){
    lines(Mset, line.CV[[k]], type="l", lwd=line_width[k], col=color_palette[k], lty=line_type[k])
  }
  
  legend_content <- NULL
  for (l in 1:length(CV)){
    legend_content <- c(legend_content, paste0("CV = ", CV[l]))
  }
  
  legend("topright", inset=0.01, legend=legend_content,
         col=color_palette[1:length(CV)], lty=rep(1,length(CV)), cex=0.8, lwd=3, box.lty=0)

}


#Example
plot.H2x2Factorial(power=0.9, test="cluster")
plot.H2x2Factorial(power=0.9, test="individual")
plot.H2x2Factorial(power=0.9, test="interaction")
plot.H2x2Factorial(power=0.9, test="joint")
plot.H2x2Factorial(power=0.9, test="I-U")
plot.H2x2Factorial(power=0.9, test="cluster", correction=F)
plot.H2x2Factorial(power=0.9, test="joint", correction=F)
plot.H2x2Factorial(power=0.9, test="I-U", correction=F)






dev.off()

# The palette with grey:
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

# The palette with black:
cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

safe_colorblind_palette <- c("#88CCEE", "#CC6677", "#DDCC77", "#117733", "#332288", "#AA4499", 
                             "#44AA99", "#999933", "#882255", "#661100", "#6699CC", "#888888")


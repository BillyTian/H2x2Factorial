#' @title Estimate the required number of clusters or the power level under a hierarchical 2x2 factorial CRT
#' 
#' @description In a hierarchical 2x2 factorial CRT with unequal cluster sizes and a continuous outcome, 
#' estimate the number of clusters based on a stated power requirement, or estimate the power level that could be achieved based on 
#' the given sample size conditions. Five types of hypothesis tests as well as their corresponding finite-sample considerations
#' could be chosen for the predictions.
#' 
#' @usage Users can input a cluster number or not, if a cluster number is given, the function will calculate the power under the chosen test 
#' and finite-sample correction, and this will ignore potential input for the power parameter; if the number of cluster is not given, the function will help 
#' predict the number of clusters based on a given power requirement or the 0.8 power level by default.
#' 
#' @param power a desired power level for sample size estmation. Default is 0.8.
#' @param n.input a number of cluster provided by the user to estimate the power that can be achieved.
#' @param alpha a type I error rate. Default is 0.05.
#' @param pi_x a proportion of clusters randomized to the cluster-level treatment. Default is 0.5, representing a balanced allocation.
#' @param pi_z a proportion of individuals randomized to the individual-level treatment within each cluster. Default is 0.5, representing a balanced allocation.
#' @param delta_x an (unstandardized) effect size of the marginal effect of the cluster-level treatment. Default is 0.25.
#' @param delta_z an (unstandardized) effect size of the marginal effect of the individual-level treatment. Default is 0.33.
#' @param delta_xz an (unstandardized) effect size of the treatment interaction. Default is 0.3.
#' @param sigma2_y a total variance of (continuous) outcome. Default is 1.
#' @param m_bar a mean cluster size. Default is 50.
#' @param CV a coefficient of variation of the cluster sizes. Default is 0.
#' @param rho an intraclass correlation coefficient characterizing the between-cluster variability. Default is 0.01.
#' @param test a character argument indicating the type of hypothesis test of interest. Supported choices include 
#' cluster (test of the cluster-level marginal effect), individual (test of the individual-level marginal effect),  
#' interaction (interaction test of the two treatments), joint (joint test of the two marginal effects), 
#' I-U (intersection-union test of the two marginal effects). Default is "cluster".  
#' @param correction a logical argument indicating whether a finite sample correction should be used. Default is T.
#' @param max_n an optional setting of a maximum number of clusters, possibly needed under test="cluster", "joint", or "I-U". Default is 1e8.
#' @param seed_mix an optional setting of a seed for conducting the simulation-based testing under a mixed distribution, possibly needed only under test="joint". Default is NULL.
#' @param size_mix a pre-specified size for the mixed distribution in the simulation-based procedure, possibly needed only under test="joint". Default is 1e4.
#' 
#' @return a integer number representing the required number of clusters or a fraction with four 
#' decimal places, with some useful messages reiterating the vital parameter choices.
#' 
#' @examples
#' 
#' H2x2Factorial(n.input=10, delta_x=0.2, delta_z=0.1, test="joint", correction=T, seed.mix=123456, CV=0.38)
#' 
#' @export
 
H2x2Factorial <- function(power=0.8, n.input=NULL,
                          alpha=0.05, pi_x=0.5, pi_z=0.5,
                          delta_x=0.25, delta_z=0.33, delta_xz=0.3,
                          sigma2_y=1,
                          m_bar=50, CV=0, rho=0.01, 
                          test="cluster",
                          correction=T, 
                          max_n=1e8,
                          seed_mix=NULL,
                          size_mix=1e4) {
  
  #Error messages
  if (!is.null(n.input)){
    if (!is.numeric(n.input) || n.input <= 0){
      stop('Inputted number of clusters must be a positive number.')
    }
  }
  
  if (!is.numeric(power) || power <= 0 || power >= 1)
    stop('Power must be numeric in (0,1).')
  
  if (!is.numeric(alpha) || alpha <= 0 || alpha >= 1)
    stop('Type I error rate must be numeric in (0,1).')

  if (!is.numeric(sigma2_y) || sigma2_y <= 0)
    stop('Total variance of outcome must be a positive number.')
  
  if (!is.numeric(m_bar) || m_bar <= 0)
    stop('Mean cluster size must be a positive number.')
  
  if (!is.numeric(CV) || CV < 0)
    stop('Coefficient of variation of the cluster sizes must be a positive number.')
  
  if (!is.numeric(rho) || rho < 0 || rho >= 1)
    stop('Intraclass correlation coefficient must be numeric in [0,1).')
  
  if ( !(test %in% c("cluster", "individual", "interaction", "joint", "I-U")) || length(test)!=1)
    stop('Type of hypothesis tests should be a single choice from "cluster", "individual", "interaction", "joint", and "I-U".')
  
  if (test=="cluster"){
    if (!is.numeric(delta_x) || delta_x <= 0)
      stop('Effect size of the cluster-level marginal treatment effect must be a positive number.')
    if (!is.numeric(pi_x) || pi_x <= 0 || pi_x >= 1)
      stop('Proportion of clusters that are randomized to the cluster-level treatment must be numeric in (0,1).')
    
  } else if (test=="individual"){
    if (!is.numeric(delta_z) || delta_z <= 0)
      stop('Effect size of the individual-level marginal treatment effect must be a positive number.')
    if (!is.numeric(pi_z) || pi_z <= 0 || pi_z >= 1)
      stop('Proportion of individuals that are randomized to the individual-level treatment must be numeric in (0,1).')
    if (correction==T)
      warning('No finite-sample correction will be done for the test of marginal effect at the individual level due to adequate degrees of freedom.')
    
  } else if (test=="interaction"){
    if (!is.numeric(delta_xz) || delta_xz <= 0)
      stop('Effect size of interaction must be a positive number.')
    if (!is.numeric(pi_x) || pi_x <= 0 || pi_x >= 1)
      stop('Proportion of clusters that are randomized to the cluster-level treatment must be numeric in (0,1).')
    if (!is.numeric(pi_z) || pi_z <= 0 || pi_z >= 1)
      stop('Proportion of individuals that are randomized to the individual-level treatment must be numeric in (0,1).')
    if (correction==T)
      warning('No finite-sample correction will be done for the interaction test due to adequate degrees of freedom.')
    
  } else if (test=="joint"){
    if (!is.numeric(delta_x) || delta_x <= 0)
      stop('Effect size of the cluster-level marginal treatment effect must be a positive number.')
    if (!is.numeric(delta_z) || delta_z <= 0)
      stop('Effect size of the individual-level marginal treatment effect must be a positive number.')
    if (!is.numeric(pi_x) || pi_x <= 0 || pi_x >= 1)
      stop('Proportion of clusters that are randomized to the cluster-level treatment must be numeric in (0,1).')
    if (!is.numeric(pi_z) || pi_z <= 0 || pi_z >= 1)
      stop('Proportion of individuals that are randomized to the individual-level treatment must be numeric in (0,1).')
    
  } else if (test=="I-U"){
    if (!is.numeric(delta_x) || delta_x <= 0)
      stop('Effect size of the cluster-level marginal treatment effect must be a positive number.')
    if (!is.numeric(delta_z) || delta_z <= 0)
      stop('Effect size of the individual-level marginal treatment effect must be a positive number.')
    if (!is.numeric(pi_x) || pi_x <= 0 || pi_x >= 1)
      stop('Proportion of clusters that are randomized to the cluster-level treatment must be numeric in (0,1).')
    if (!is.numeric(pi_z) || pi_z <= 0 || pi_z >= 1)
      stop('Proportion of individuals that are randomized to the individual-level treatment must be numeric in (0,1).') 
  }
  
  if (!is.logical(correction))
    stop('Finite sample correction indicator should be a logical argument.')
  
  
  if (!is.numeric(max_n) || max_n <= 0)
    stop('Maximum number of clusters must be a positive number.')
  
  if (!is.null(seed_mix)){
    if (!is.numeric(seed_mix))
    stop('User-defined seed under the finite-sample corrected joint test must be numeric.')
  }
  
  if (!is.numeric(size_mix) || size_mix <= 0)
    stop('Sample size for simulating the mix distribution under the finite-sample corrected joint test must be a positive number.')
  
  
  #Re-iterate the given effect sizes and the chosen test
  if (test=="cluster"){
    message('Hypothesis test: Test of the cluster-level marginal treatment effect.')
    message(paste0('Given effect size: ', delta_x, " for the cluster-level marginal effect."))
 
  } else if (test=="individual"){
    message('Hypothesis test: Test of the individual-level marginal treatment effect.')
    message(paste0('Given effect size: ', delta_z, " for the individual-level marginal effect."))

  } else if (test=="interaction"){
    message('Hypothesis test: Interaction test of the two treatment effects.')
    message(paste0('Given effect size: ', delta_xz, " for the interaction effect."))

  } else if (test=="joint"){
    message('Hypothesis test: Joint test of the marginal effects.')
    message(paste0('Given effect sizes: ', delta_x, " for the cluster-level marginal effect, ", delta_z, " for the individual-level marginal effect."))

  } else if (test=="I-U"){
    message('Hypothesis test: Intersection-union test of the marginal effects.')
    message(paste0('Given effect sizes: ', delta_x, " for the cluster-level marginal effect, ", delta_z, " for the individual-level marginal effect."))
 
  }
  
  
  a <- alpha
  b <- 1-power
  z_a <- qnorm(1-a/2)
  z_b <- qnorm(1-b)
  
  eta1 <- -m_bar*rho/(1+(m_bar-1)*rho)*(1-CV^2*m_bar*rho*(1-rho)/((1+(m_bar-1)*rho)^2))
  eta2 <- m_bar*(1-rho)/(1+(m_bar-1)*rho)*(1-CV^2*m_bar*rho*(1-rho)/((1+(m_bar-1)*rho)^2))-m_bar
  omega_x <- sigma2_y*(1-rho)/((m_bar+eta2)*pi_x*(1-pi_x))
  omega_z <- sigma2_y*(1-rho)/((m_bar+eta1)*pi_z*(1-pi_z))
  omega_xz <- sigma2_y*(1-rho)/((m_bar+eta1)*pi_x*(1-pi_x)*pi_z*(1-pi_z))
  
  
  pred.power <- function(n.input) {
    ### Test (A1): Test of the cluster-level marginal effect
    if (test=="cluster"){
      if (correction==F){
        
        pred.power <- pnorm(sqrt(n.input*delta_x^2/omega_x)-z_a)
        message("A Wald z-test is used without finite-sample correction.")
        
      } else if (correction==T){
        
        pred.power <- pt(qt(1-a/2, n.input-2), n.input-2, ncp=delta_x/sqrt(omega_x/n.input), lower.tail = F) 
        + pt(qt(a/2, n.input-2), n.input-2, ncp=delta_x/sqrt(omega_x/n.input))
        message("A t-test is used for finite-sample correction.")
        
      }
    }
    
    ### Test (A2): Test of the individual-level marginal effect
    if (test=="individual"){
      
      pred.power <- pnorm(sqrt(n.input*delta_z^2/omega_z)-z_a)
      message("A Wald z-test is automatically used.")
      
    }
    
    ### Test (B): Test of interaction
    if (test=="interaction"){
      
      pred.power <- pnorm(sqrt(n.input*delta_xz^2/omega_xz)-z_a)
      message("A Wald z-test is automatically used.")
      
    }
    
    ### Test (C): Joint test for the two marginal effects
    if (test=="joint"){
      if (correction==F){
        
        theta <- n.input*(delta_x^2/omega_x + delta_z^2/omega_z)
        pred.power <- pchisq(qchisq(1-a, 2), 2, ncp = theta, lower.tail = F)
        message("A Chi-square test is used without finite-sample correction.")
        
      } else if (correction==T){
        
        set.seed(seed_mix)
        #Simulate the mixed distribution (CENTRAL) to identify rejection region bound
        f.distn <- rf(size_mix, 1, n.input-2)
        chisq.distn <- rchisq(size_mix, 1)
        mix.distn <- f.distn + chisq.distn
        crt.value <- quantile(mix.distn, 1-a)
        
        #Simulate the mixed distribution (NONCENTRAL VERSION) for power calculation
        nc.f.distn <- rf(size_mix, 1, n.input-2, ncp = n.input*delta_x^2/omega_x)
        nc.chisq.distn <- rchisq(size_mix, 1, ncp = n.input*delta_z^2/omega_z)
        nc.mix.distn <- nc.f.distn + nc.chisq.distn
        
        pred.power <- mean(nc.mix.distn>crt.value)
        message("A simulation-based mixed F-Chi-square test is used for finite-sample correction.")
        
      }
    }
    
    ### Test (D): Intersection-union test for the two marginal effects
    if (test=="I-U"){
      if (correction==F){
        
        wmean.c <- sqrt(n.input)*delta_x/sqrt(omega_x)
        wmean.i <- sqrt(n.input)*delta_z/sqrt(omega_z)
        pred.power <- pnorm(qnorm(1-a/2), mean = wmean.c, lower.tail = F)*pnorm(qnorm(1-a/2), mean = wmean.i, lower.tail = F) 
        + pnorm(qnorm(1-a/2), mean = wmean.c, lower.tail = F)*pnorm(qnorm(a/2), mean = wmean.i) 
        + pnorm(qnorm(a/2), mean = wmean.c)*pnorm(qnorm(1-a/2), mean = wmean.i, lower.tail = F) 
        + pnorm(qnorm(a/2), mean = wmean.c)*pnorm(qnorm(a/2), mean = wmean.i)
        message("A z-based intersection-union test is used without finite-sample correction.")
        
        
      } else if (correction==T){
        
        c.ncp <- sqrt(n.input)*delta_x/sqrt(omega_x)
        i.mean <- sqrt(n.input)*delta_z/sqrt(omega_z)
        pred.power <- pt(qt(1-a/2, n.input-2), n.input-2, c.ncp, lower.tail = F)*pnorm(qnorm(1-a/2), mean = i.mean, lower.tail = F) 
        + pt(qt(1-a/2, n.input-2), n.input-2, c.ncp, lower.tail = F)*pnorm(qnorm(a/2), mean = i.mean) 
        + pt(qt(a/2, n.input-2), n.input-2, c.ncp)*pnorm(qnorm(1-a/2), mean = i.mean, lower.tail = F) 
        + pt(qt(a/2, n.input-2), n.input-2, c.ncp)*pnorm(qnorm(a/2), mean = i.mean)
        message("A mixed t- and z-based intersection-union test is used for finite-sample correction.")
        
      }
    }
    return(pred.power)
  }
  
  
  
  #Function to estimate number of clusters based on a required power level
  cluster.number <- function(power) {
    ### Test (A1): Test of cluster-level marginal effect
    if (test=="cluster"){
      if(correction==F){
        message("A Wald z-test is used without finite-sample-correction.")
      } else if (correction==T){
        message("A t-test is used for finite-sample correction.")
      }
      n.out <- marginal.cluster(m_bar, CV, power, delta_x, rho, pi_x, correction, sigma2_y, a, max_n, z_a, z_b)
    }
    
    ### Test (A2): Test of individual-level marginal effect
    if (test=="individual"){
      message("A Wald z-test is automatically used.")
      n.out <- marginal.ind(m_bar, CV, power, delta_z, rho, pi_z, sigma2_y, z_a, z_b)
    }
    
    ### Test (B): Test of interaction
    if (test=="interaction"){
      message("A Wald z-test is automatically used.")
      n.out <- interaction(m_bar, CV, power, delta_xz, rho, pi_x, pi_z, sigma2_y, z_a, z_b)
    }
    
    ### Test (C): Joint test of marginal effects on both treatment levels
    if (test=="joint"){
      if(correction==F){
        message("A Chi-square test is used without finite-sample correction.")
      } else if (correction==T){
        message("A simulation-based mixed F-Chi-square test is used for finite-sample correction.")
      }
      n.out <- joint(m_bar, CV, power, delta_x, delta_z, rho, pi_x, pi_z, correction, sigma2_y, a, max_n, seed_mix, size_mix)
    }
    
    ### Test (D): Intersection-Union test of marginal effects on both treatment levels
    if (test=="I-U"){
      if(correction==F){
        message("A z-based intersection-union test is used without finite-sample correction.")
      } else if (correction==T){
        message("A mixed t- and z-based intersection-union test is used for finite sample correction.")
      }
      n.out <- IU(m_bar, CV, power, delta_x, delta_z, rho, pi_x, pi_z, correction, sigma2_y, a, max_n)
    }
    return(n.out)
  }

  #When the user input n, the potentially inputted power will be ignored, and the formula will always give the predicted power;
  if (!is.null(n.input)){
    ans <- round(pred.power(n.input), 4)
    message(paste0("Predicted power for the provided ", n.input, " clusters: ", ans, "."))
  #When the user does not input n, the formula will give the required number of clusters to hit the power
  } else {
    ans <- cluster.number(power)
    message(paste0("Required number of clusters to achieve ", power, " power: ", ans, "."))
  }
  
  return(ans)
}


#Examples
H2x2Factorial(power=0.9, delta_x=0.2, delta_z=0.1, test="cluster", CV=0.75)
H2x2Factorial(power=0.9, delta_x=0.2, delta_z=0.1, test="cluster", correction=F, CV=0.75)
H2x2Factorial(power=0.9, delta_x=0.2, delta_z=0.1, test="individual", CV=0.75)
H2x2Factorial(power=0.9, delta_xz=0.3, test="interaction", CV=0.75)
H2x2Factorial(power=0.9, delta_x=0.2, delta_z=0.1, test="joint", correction=F, CV=0.75)
H2x2Factorial(power=0.9, delta_x=0.2, delta_z=0.1, test="joint", correction=T, seed.mix=123456, CV=0.38)
H2x2Factorial(power=0.9, delta_x=0.2, delta_z=0.1, test="I-U", correction=F, CV=0.75)
H2x2Factorial(power=0.9, delta_x=0.2, delta_z=0.1, test="I-U", correction=T, CV=0.38)

H2x2Factorial(n.input=10, delta_x=0.2, delta_z=0.1, test="cluster", CV=0.75)
H2x2Factorial(n.input=10, delta_x=0.2, delta_z=0.1, test="cluster", correction=F, CV=0.75)
H2x2Factorial(n.input=10, delta_x=0.2, delta_z=0.1, test="individual", CV=0.75)
H2x2Factorial(n.input=10, delta_xz=0.3, test="interaction", CV=0.75)
H2x2Factorial(n.input=10, delta_x=0.2, delta_z=0.1, test="joint", correction=F, CV=0.75)
H2x2Factorial(n.input=10, delta_x=0.2, delta_z=0.1, test="joint", correction=T, seed.mix=123456, CV=0.38)
H2x2Factorial(n.input=10, delta_x=0.2, delta_z=0.1, test="I-U", correction=F, CV=0.75)
H2x2Factorial(n.input=10, delta_x=0.2, delta_z=0.1, test="I-U", correction=T, CV=0.38)



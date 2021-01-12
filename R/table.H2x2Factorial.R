#' @description  A table summarizing the required number of clusters and the predicted power level 
#' based on a constellation of design parameters
#' 
#' @usage if the user want a vector of power or other parameters invoking the need 
#' for multiple tables, an external loop could be easily written using this function.
#' 
#' @param power a target power level. Default is 0.8.
#' @param alpha a type I error rate. Default is 0.05.
#' @param pi_x a proportion of clusters randomized to the cluster-level treatment. Default is 0.5, representing a balanced allocation.
#' @param pi_z a proportion of individuals randomized to the individual-level treatment within each cluster. Default is 0.5, representing a balanced allocation.
#' @param delta_x an effect size of the marginal effect of the cluster-level treatment.
#' @param delta_z an effect size of the marginal effect of the individual-level treatment.
#' @param delta_xz an effect size of the treatment interaction.
#' @param sigma2_y a total variance of (continuous) outcome. Default is 1.
#' @param m_bar a vector of mean cluster sizes.
#' @param CV a vector of the coefficients of variation of the cluster sizes.
#' @param rho a vector of the intraclass correlation coefficients.
#' @param test a character argument indicating the type of hypothesis test of interest. Supported choices include 
#' cluster (test of the cluster-level marginal effect), individual (test of the individual-level marginal effect),  
#' interaction (interaction test of the two treatments), joint (joint test of the two marginal effects), I-U (Intersection-union test of the two marginal effects).
#' Default is "cluster".  
#' @param correction a logical argument indicating whether a finite sample correction should be used. Default is T.
#' @param max_n an optional setting of a maximum number of clusters, possibly needed under test="cluster", "joint", or "I-U". Default is 1e8.
#' @param seed_mix an optional setting of a seed for conducting the simulation-based testing under a mixed distribution, possibly needed only under test="joint". Default is NULL.
#' @param size_mix a pre-specified size for the mixed distribution in the simulation-based procedure, possibly needed only under test="joint". Default is 1e4.
#' 
#' @return a table with the m_bar, rho, and CV varied in a factorial setting.
#' 
#' @examples
#' 
#' table.H2x2Factorial(delta_x=0.2, delta_z=0.1, m_bar=c(10,50,100), CV=c(0, 0.3, 0.5), rho=c(0.01, 0.1), test="cluster")
#' 
#' @export
#'
#'  
table.H2x2Factorial <- function(power=0.8, alpha=0.05, 
                                pi_x=0.5, pi_z=0.5,
                                delta_x, delta_z, delta_xz,
                                sigma2_y=1,
                                m_bar, CV, rho, 
                                test="cluster",
                                correction=T, 
                                max_n=1e8,
                                seed_mix=NULL,
                                size_mix=1e4) {

  
  #Re-iterate the given effect sizes and the chosen test
  if (test=="cluster"){
    print('Hypothesis test: Test of the cluster-level marginal treatment effect.')
    print(paste0('Given effect size: ', delta_x, " for the cluster-level marginal effect."))
    #print(paste0("Based on the given effect sizes, the number of output tables is ", length(delta_x)))
    
  } else if (test=="individual"){
    print('Hypothesis test: Test of the individual-level marginal treatment effect.')
    print(paste0('Given effect size: ', delta_z, " for the individual-level marginal effect."))
    #print(paste0("Based on the given effect sizes, the number of output tables is ", length(delta_z)))
    
  } else if (test=="interaction"){
    print('Hypothesis test: Interaction test of the two treatment effects.')
    print(paste0('Given effect size: ', delta_xz, " for the interaction effect."))
    #print(paste0("Based on the given effect sizes, the number of output tables is ", length(delta_xz)))
    
  } else if (test=="joint"){
    print('Hypothesis test: Joint test of the marginal effects.')
    print(paste0('Given effect sizes: ', delta_x, " for the cluster-level marginal effect, ", delta_z, " for the individual-level marginal effect."))
    #print(paste0("Based on the given effect sizes, the number of output tables is ", length(delta_x)*length(delta_z)))
    
  } else if (test=="I-U"){
    print('Hypothesis test: Intersection-union test of the marginal effects.')
    print(paste0('Given effect sizes: ', delta_x, " for the cluster-level marginal effect, ", delta_z, " for the individual-level marginal effect."))
    #print(paste0("Based on the given effect sizes, the number of output tables is ", length(delta_x)*length(delta_z)))
    
  }
  
  #test
  #delta_x <- c(0.2, 0.1)
  #delta_z <- c(0.1, 0.05)
  #delta_xz <- c(0.1, 0.2)
  
  #m_bar <- c(10, 50, 100)
  #CV <- c(0, 0.3, 0.6)
  #rho <- c(0.01, 0.1)
  
  #delta_x.vector <- delta_x
  #delta_z.vector <- delta_z
  #delta_xz.vector <- delta_xz
  
  m_bar.vector <- m_bar 
  CV.vector <- CV 
  rho.vector <- rho
  
  table <- NULL
  for (m_bar.i in 1:length(m_bar.vector)){
    for (rho.i in 1:length(rho.vector)){
      for (CV.i in 1:length(CV.vector)){
        table <- rbind(table, c(m_bar[m_bar.i], rho[rho.i], CV[CV.i]))
      }
    }
  }
  m_bar <- table[,1]
  rho <- table[,2]
  CV <- table[,3]
  

  a <- alpha
  b <- 1-power
  z_a <- qnorm(1-a/2)
  z_b <- qnorm(1-b)
  
  eta1 <- -m_bar*rho/(1+(m_bar-1)*rho)*(1-CV^2*m_bar*rho*(1-rho)/((1+(m_bar-1)*rho)^2))
  eta2 <- m_bar*(1-rho)/(1+(m_bar-1)*rho)*(1-CV^2*m_bar*rho*(1-rho)/((1+(m_bar-1)*rho)^2))-m_bar
  omega_x <- sigma2_y*(1-rho)/((m_bar+eta2)*pi_x*(1-pi_x))
  omega_z <- sigma2_y*(1-rho)/((m_bar+eta1)*pi_z*(1-pi_z))
  omega_xz <- sigma2_y*(1-rho)/((m_bar+eta1)*pi_x*(1-pi_x)*pi_z*(1-pi_z))
  
  ### Test (A1): Test of cluster-level marginal effect
  if (test=="cluster"){
    if(correction==F){
      print("A Wald z-test is used without finite-sample correction.")
      
      n <- (z_a+z_b)^2*omega_x/delta_x^2
      n.final <- ceiling(n)
      pred.power <- pnorm(sqrt(n.final*delta_x^2/omega_x)-z_a)
      
    } else if (correction==T){
      print("A t-test is used for finite-sample correction.")
      
      cluster.n <- function(parameter){
        m_bar <- parameter[1]
        rho <- parameter[2]
        CV <- parameter[3]
        eta2 <- m_bar*(1-rho)/(1+(m_bar-1)*rho)*(1-CV^2*m_bar*rho*(1-rho)/((1+(m_bar-1)*rho)^2))-m_bar
        omega_x <- sigma2_y*(1-rho)/((m_bar+eta2)*pi_x*(1-pi_x))
        n <- 3
        try.power <- 0
        while ((try.power < power) & (n < max_n)){
          try.power <- pt(qt(1-a/2, n-2), n-2, ncp=delta_x/sqrt(omega_x/n), lower.tail = F) + pt(qt(a/2, n-2), n-2, ncp=delta_x/sqrt(omega_x/n))
          n <- n+1
        }
        return(c(n, try.power))
      }
      
      cluster.pred <- NULL
      for (i in 1:nrow(table)){
        cluster.pred <- rbind(cluster.pred, cluster.n(parameter=unlist(table[i,])))
      }
      n.final <- cluster.pred[,1]
      pred.power <- cluster.pred[,2]
    }
  }
    
  ### Test (A2): Test of individual-level marginal effect
  if (test=="individual"){
    print("A Wald z-test is automatically used.")
    n <- (z_a+z_b)^2*omega_z/delta_z^2
    n.final <- ceiling(n)
    pred.power <- pnorm(sqrt(n.final*delta_z^2/omega_z)-z_a)
  }  
    
  ### Test (B): Interaction test
  if (test=="interaction"){
    print("A Wald z-test is automatically used.")
    n <- (z_a+z_b)^2*omega_xz/delta_xz^2
    n.final <- ceiling(n)
    pred.power <- pnorm(sqrt(n.final*delta_xz^2/omega_xz)-z_a)
  }   
    
  ### Test (C): Joint test of marginal effects on both treatment levels
  if (test=="joint"){
    if(correction==F){
      print("A Chi-square test is used without finite-sample correction.")

      joint.n <- function(parameter){
        m_bar <- parameter[1]
        rho <- parameter[2]
        CV <- parameter[3]
        eta1 <- -m_bar*rho/(1+(m_bar-1)*rho)*(1-CV^2*m_bar*rho*(1-rho)/((1+(m_bar-1)*rho)^2))
        eta2 <- m_bar*(1-rho)/(1+(m_bar-1)*rho)*(1-CV^2*m_bar*rho*(1-rho)/((1+(m_bar-1)*rho)^2))-m_bar
        omega_x <- sigma2_y*(1-rho)/((m_bar+eta2)*pi_x*(1-pi_x))
        omega_z <- sigma2_y*(1-rho)/((m_bar+eta1)*pi_z*(1-pi_z))
        n <- 2
        try.power <- 0
        while ((try.power < power) & (n < max_n)){
          theta <- n*(delta_x^2/omega_x + delta_z^2/omega_z)
          try.power <- pchisq(qchisq(1-a, 2), 2, ncp = theta, lower.tail = F)
          n <- n+1
        }
        return(c(n, try.power))
      }
      
      joint.pred <- NULL
      for (i in 1:nrow(table)){
        joint.pred <- rbind(joint.pred, joint.n(parameter=unlist(table[i,])))
      }
      
      n.final <- joint.pred[,1]
      pred.power <- joint.pred[,2]
      
    } else if (correction==T){
      print("A simulation-based mixed F-Chi-square test is used for finite-sample correction.")
      
      joint.n <- function(parameter){
        m_bar <- parameter[1]
        rho <- parameter[2]
        CV <- parameter[3]
        eta1 <- -m_bar*rho/(1+(m_bar-1)*rho)*(1-CV^2*m_bar*rho*(1-rho)/((1+(m_bar-1)*rho)^2))
        eta2 <- m_bar*(1-rho)/(1+(m_bar-1)*rho)*(1-CV^2*m_bar*rho*(1-rho)/((1+(m_bar-1)*rho)^2))-m_bar
        omega_x <- sigma2_y*(1-rho)/((m_bar+eta2)*pi_x*(1-pi_x))
        omega_z <- sigma2_y*(1-rho)/((m_bar+eta1)*pi_z*(1-pi_z))
        n <- 3
        try.power <- 0
        while ((try.power < power) & (n < max_n)){
          set.seed(seed_mix)
          f.distn <- rf(size_mix, 1, n-2)
          chisq.distn <- rchisq(size_mix, 1)
          mix.distn <- f.distn + chisq.distn
          crt.value <- quantile(mix.distn, 1-a)
          
          nc.f.distn <- rf(size_mix, 1, n-2, ncp = n*delta_x^2/omega_x)
          nc.chisq.distn <- rchisq(size_mix, 1, ncp = n*delta_z^2/omega_z)
          nc.mix.distn <- nc.f.distn + nc.chisq.distn
          try.power <- mean(nc.mix.distn>crt.value)
          n <- n+1
        }
        return(c(n, try.power))
      }
      
      joint.pred <- NULL
      for (i in 1:nrow(table)){
        joint.pred <- rbind(joint.pred, joint.n(parameter=unlist(table[i,])))
      }
      n.final <- joint.pred[,1]
      pred.power <- joint.pred[,2]
    }
  }  
    
  ### Test (D): Intersection-Union test of marginal effects on both treatment levels
  if (test=="I-U"){
    if(correction==F){
      print("A z-based intersection-union test is used without finite-sample correction.")
      
      IU.n <- function(parameter){
        m_bar <- parameter[1]
        rho <- parameter[2]
        CV <- parameter[3]
        eta1 <- -m_bar*rho/(1+(m_bar-1)*rho)*(1-CV^2*m_bar*rho*(1-rho)/((1+(m_bar-1)*rho)^2))
        eta2 <- m_bar*(1-rho)/(1+(m_bar-1)*rho)*(1-CV^2*m_bar*rho*(1-rho)/((1+(m_bar-1)*rho)^2))-m_bar
        omega_x <- sigma2_y*(1-rho)/((m_bar+eta2)*pi_x*(1-pi_x))
        omega_z <- sigma2_y*(1-rho)/((m_bar+eta1)*pi_z*(1-pi_z))
        n <- 2
        try.power <- 0
        while ((try.power < power) & (n < max_n)){
          wmean.c <- sqrt(n)*delta_x/sqrt(omega_x)
          wmean.i <- sqrt(n)*delta_z/sqrt(omega_z)
          try.power <- pnorm(qnorm(1-a/2), mean = wmean.c, lower.tail = F)*pnorm(qnorm(1-a/2), mean = wmean.i, lower.tail = F) + pnorm(qnorm(1-a/2), mean = wmean.c, lower.tail = F)*pnorm(qnorm(a/2), mean = wmean.i) + pnorm(qnorm(a/2), mean = wmean.c)*pnorm(qnorm(1-a/2), mean = wmean.i, lower.tail = F) + pnorm(qnorm(a/2), mean = wmean.c)*pnorm(qnorm(a/2), mean = wmean.i)
          n <- n+1
        }
        return(c(n, try.power))
      }
      
      IU.pred <- NULL
      for (i in 1:nrow(table)){
        IU.pred <- rbind(IU.pred, IU.n(parameter=unlist(table[i,])))
      }
      n.final <- IU.pred[,1]
      pred.power <- IU.pred[,2]
      
    } else if (correction==T){
      print("A mixed t- and z-based intersection-union test is used for finite-sample correction.")
      
      IU.n <- function(parameter){
        m_bar <- parameter[1]
        rho <- parameter[2]
        CV <- parameter[3]
        eta1 <- -m_bar*rho/(1+(m_bar-1)*rho)*(1-CV^2*m_bar*rho*(1-rho)/((1+(m_bar-1)*rho)^2))
        eta2 <- m_bar*(1-rho)/(1+(m_bar-1)*rho)*(1-CV^2*m_bar*rho*(1-rho)/((1+(m_bar-1)*rho)^2))-m_bar
        omega_x <- sigma2_y*(1-rho)/((m_bar+eta2)*pi_x*(1-pi_x))
        omega_z <- sigma2_y*(1-rho)/((m_bar+eta1)*pi_z*(1-pi_z))
        n <- 3
        try.power <- 0
        while ((try.power < power) & (n < max_n)){
          c.ncp <- sqrt(n)*delta_x/sqrt(omega_x)
          i.mean <- sqrt(n)*delta_z/sqrt(omega_z)
          try.power <- pt(qt(1-a/2, n-2), n-2, c.ncp, lower.tail = F)*pnorm(qnorm(1-a/2), mean = i.mean, lower.tail = F) + pt(qt(1-a/2, n-2), n-2, c.ncp, lower.tail = F)*pnorm(qnorm(a/2), mean = i.mean) + pt(qt(a/2, n-2), n-2, c.ncp)*pnorm(qnorm(1-a/2), mean = i.mean, lower.tail = F) + pt(qt(a/2, n-2), n-2, c.ncp)*pnorm(qnorm(a/2), mean = i.mean)
          n <- n+1
        }
        return(c(n, try.power))
      }
      
      IU.pred <- NULL
      for (i in 1:nrow(table)){
        IU.pred <- rbind(IU.pred, IU.n(parameter=unlist(table[i,])))
      }
      
      n.final <- IU.pred[,1]
      pred.power <- IU.pred[,2]
    }
  }  
  
  table <- data.frame(cbind(m_bar, rho, CV, n.final, pred.power))
  names(table) <- c("m_bar", "rho", "CV", "n", "predicted power")
  return(table)
}


#Examples
table.H2x2Factorial(delta_x=0.2, delta_z=0.1, m_bar=c(10,50,100), CV=c(0, 0.3, 0.5), rho=c(0.01, 0.1), test="cluster")
table.H2x2Factorial(delta_x=0.2, delta_z=0.1, m_bar=c(10,50,100), CV=c(0, 0.3, 0.5), rho=c(0.01, 0.1), test="cluster", correction=F)
table.H2x2Factorial(delta_x=0.2, delta_z=0.1, m_bar=c(10,50,100), CV=c(0, 0.3, 0.5), rho=c(0.01, 0.1), test="individual")
table.H2x2Factorial(delta_xz=0.2, m_bar=c(10,50,100), CV=c(0, 0.3, 0.5), rho=c(0.01, 0.1), test="interaction")
table.H2x2Factorial(delta_x=0.2, delta_z=0.1, m_bar=c(50,100), CV=c(0, 0.3, 0.6, 0.9), rho=c(0.02, 0.05, 0.1), test="joint")
table.H2x2Factorial(delta_x=0.2, delta_z=0.1, m_bar=c(10,50,100), CV=c(0, 0.3, 0.5), rho=c(0.01, 0.1), test="I-U")
table.H2x2Factorial(delta_x=0.2, delta_z=0.1, m_bar=c(50,100), CV=c(0, 0.3, 0.6, 0.9), rho=c(0.02, 0.05, 0.1), test="joint", correction=F)
table.H2x2Factorial(delta_x=0.2, delta_z=0.1, m_bar=c(10,50,100), CV=c(0, 0.3, 0.5), rho=c(0.01, 0.1), test="I-U", correction=F)





#External loop for outputting multiple tables due to inputs of effect size vectors
#Example
delta_x <- c(0.4, 0.2)
delta_z <- c(0.2, 0.1)

tables <- NULL
for (delta_x.i in 1:length(delta_x)){
  for (delta_z.i in 1:length(delta_z)){

  }
}
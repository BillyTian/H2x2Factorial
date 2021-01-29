#' @title H2x2Factorial Sample Size and Power Calculation
#'
#' @description The function \code{calc.H2x2Factorial} estimates the required number of clusters or the achieved power level under different types of
#' hypothesis tests in a hierarchical 2x2 factorial trial with unequal cluster sizes and a continuous outcome. Five types of hypothesis tests
#' as well as their corresponding finite-sample considerations could be chosen for the predictions. Users may input an optional cluster number through the
#' \code{n.input} argument. When this number is provided, the function will calculate the power under a chosen hypothesis test as well as a
#' finite-sample correction if specified, and the function will ignore the potential input for the power parameter; When the number of clusters is not
#' provided, the function will calculate the required number of clusters based on a given power threshold, which is set to 0.8 by default.
#'
#' @usage
#' calc.H2x2Factorial(power=0.8, n_input=NULL, alpha=0.05,
#'                    pi_x=0.5, pi_z=0.5,
#'                    delta_x=0.25, delta_z=0.33, delta_xz=0.3, sigma2_y=1,
#'                    m_bar=50, CV=0, rho=0,
#'                    test="cluster", correction=FALSE,
#'                    max_n=1e8, seed_mix=NULL, size_mix=1e4,
#'                    verbose=TRUE)
#'
#' @param power a numeric value between 0 and 1 as the desired power level for sample size estimation. Default is \code{0.8}.
#' @param n_input a number of cluster provided by the user to estimate the power that can be achieved. Default is \code{NULL}.
#' @param alpha a numeric value between 0 and 1 as the type I error rate. Default is \code{0.05}.
#' @param pi_x a numeric value between 0 and 1 as the proportion of clusters randomized to the cluster-level treatment. Default is \code{0.5}, representing a balanced allocation.
#' @param pi_z a numeric value between 0 and 1 as the proportion of individuals randomized to the individual-level treatment within each cluster. Default is \code{0.5}, representing a balanced allocation.
#' @param delta_x a nonzero numeric value for the (unstandardized) effect size of the marginal cluster-level treatment effect. Default is \code{0.25}, which is the hypothetical value for the example in the referenced paper.
#' @param delta_z a nonzero numeric value for the (unstandardized) effect size of the marginal individual-level treatment effect. Default is \code{0.33}, which is the hypothetical value for the example in the referenced paper.
#' @param delta_xz a nonzero numeric value for the (unstandardized) effect size of the interaction effect of the two treatments. Default is \code{0.3}, which is the hypothetical value for the example in the referenced paper.
#' @param sigma2_y a positive numeric value for the total variance of the continuous outcome. Default is \code{1}.
#' @param m_bar a numeric value larger than 2 for the mean cluster size. Default is \code{50}.
#' @param CV a positive numeric value as the coefficient of variation of the cluster sizes. Default is \code{0}, representing equal cluster sizes.
#' @param rho a numeric value between 0 and 1 as the intraclass correlation coefficient characterizing the between-cluster variability. Default is \code{0}.
#' @param test a character argument indicating the type of hypothesis test of interest. Supported values include
#' \code{"cluster"} (test for marginal cluster-level treatment effect), \code{"individual"} (test for marginal individual-level treatment effect),
#' \code{"interaction"} (interaction test for the two treatments), \code{"joint"} (joint test for the two marginal treatment effects),
#' \code{"I-U"} (intersection-union test for the two marginal effects). Default is \code{"cluster"}.
#' @param correction a logical argument indicating whether a finite sample correction should be used. Default is \code{FALSE}.
#' @param max_n an optional setting of a maximum number of clusters, which is only functional under \code{test="cluster"}, \code{"joint"}, or \code{"I-U"}. Default is \code{1e8}.
#' @param seed_mix an optional setting of a seed for conducting the simulation-based testing under a mixed distribution, which is only functional under \code{test="joint"}. Default is \code{NULL}.
#' @param size_mix a pre-specified size for the mixed distribution in the simulation-based procedure, which is only needed under \code{test="joint"}. Default is \code{1e4}.
#' @param verbose a logical argument indicating whether the parameter reiterations and supplementary messages should be presented or suppressed. Default is \code{TRUE}.
#'
#' @details
#' Given the input parameters, our method will firstly compute the variances of the effects of interest based on Generalized Least Square estimators and large-sample approximations.
#' Then, the variances are used to build up either the classic sample size formulas (for the separate tests for marginal treatment effects and the interaction test) or the power formulas (for the simultaneous tests and the corrected tests), which help to deliver both the sample size and power calculations.
#' Without finite-sample considerations, the two separate tests for marginal treatment effects as well as the interaction test use the two-sided Wald z-test, the joint test use the Chi-square test, and the intersection-union (I-U) test use also a two-sided z-based test.
#' With \code{correction=T}, finite-sample corrections are customized for the three types of tests involving the marginal cluster-level treatment effect: For the test for marginal cluster-level treatment effect, a two-sided t-test is used;
#' For the joint test, a simulation-based mixed F-chi-square test is used; For the I-U test, a two-sided mixed t- and z-based test is used.
#' For the finite-sample corrected joint test, since there does not exist the required parametric distribution, we offer a simulation-based method to generate the null and alternative distributions, and we use the simulated distributions to compute the power and required sample size.
#' A seed should be set via \code{seed_mix} for this random process to promote reproducibility, and this is only needed under the joint test with finite-sample correction.
#' The five types of \code{test} and \code{correction} are defined in Tian (under review).
#'
#' @return \code{calc.H2x2Factorial} returns an integer representing the required number of clusters or a decimal representing the power that can be achieved by the provided
#' sample size, with some useful and suppressible messages elaborating vital parameter choices and results (the power will be displayed in 4 decimal places; the messages can be suppressed via \code{verbose=FALSE}).
#'
#' @export
#'
#' @examples
#' #Predict the actual power of a joint test when the number of clusters is 10
#' joint.power <- calc.H2x2Factorial(n_input=10,
#'                                   delta_x=0.2, delta_z=0.1,
#'                                   rho=0.1, CV=0.38,
#'                                   test="joint", correction=TRUE, seed_mix=123456, verbose=FALSE)
#' print(joint.power)
#'
#' @importFrom stats qnorm pnorm qt pt qchisq pchisq rchisq qf pf rf quantile
#'
calc.H2x2Factorial <- function(power=0.8, n_input=NULL,
                               alpha=0.05, pi_x=0.5, pi_z=0.5,
                               delta_x=0.25, delta_z=0.33, delta_xz=0.3,
                               sigma2_y=1,
                               m_bar=50, CV=0, rho=0,
                               test="cluster",
                               correction=FALSE,
                               max_n=1e8,
                               seed_mix=NULL,
                               size_mix=1e4,
                               verbose=TRUE) {

  #Error messages
  if (!is.null(n_input)){
    if (!is.numeric(n_input) || n_input <= 0 || length(n_input)!=1){
      stop('Inputted number of clusters must be a single positive number')
    }
  }

  if (!is.numeric(power) || power <= 0 || power >= 1 || length(power)!=1)
    stop('Target power must be a single number in (0,1)')

  if (!is.numeric(alpha) || alpha <= 0 || alpha >= 1 || length(alpha)!=1)
    stop('Type I error rate must be a single number in (0,1)')

  if (!is.numeric(sigma2_y) || sigma2_y <= 0 || length(alpha)!=1)
    stop('Total variance of the outcome must be a single positive number')

  if (!is.numeric(m_bar) || m_bar <= 2 || length(m_bar)!=1)
    stop('Mean cluster size must be a single number larger than 2')

  if (!is.numeric(CV) || CV < 0 || length(CV)!=1)
    stop('Coefficient of variation of the cluster sizes must be a single positive number')

  if (!is.numeric(rho) || rho < 0 || rho >= 1 || length(rho)!=1)
    stop('Intraclass correlation coefficient must be a single number in [0,1)')

  if ( !(test %in% c("cluster", "individual", "interaction", "joint", "I-U")) || length(test)!=1)
    stop('Type of hypothesis tests should be a single choice from "cluster", "individual", "interaction", "joint", and "I-U"')

  if (test=="cluster"){
    if (!is.numeric(delta_x) || delta_x == 0 || length(delta_x)!=1)
      stop('Effect size of the marginal cluster-level treatment effect must be a nonzero number')
    if (!is.numeric(pi_x) || pi_x <= 0 || pi_x >= 1 || length(pi_x)!=1)
      stop('Proportion of clusters that are randomized to the cluster-level treatment arm must be a single number in (0,1)')

  } else if (test=="individual"){
    if (!is.numeric(delta_z) || delta_z == 0 || length(delta_z)!=1)
      stop('Effect size of the marginal individual-level treatment effect must be a nonzero number')
    if (!is.numeric(pi_z) || pi_z <= 0 || pi_z >= 1 || length(pi_z)!=1)
      stop('Proportion of individuals that are randomized to the individual-level treatment arm must be a single number in (0,1)')
    if (correction==TRUE)
      message('No finite-sample correction will be done for the test for marginal individual-level treatment effect due to adequate degrees of freedom')

  } else if (test=="interaction"){
    if (!is.numeric(delta_xz) || delta_xz == 0 || length(delta_xz)!=1)
      stop('Effect size of the interaction effect must be a nonzero number')
    if (!is.numeric(pi_x) || pi_x <= 0 || pi_x >= 1 || length(pi_x)!=1)
      stop('Proportion of clusters that are randomized to the cluster-level treatment arm must be a single number in (0,1)')
    if (!is.numeric(pi_z) || pi_z <= 0 || pi_z >= 1 || length(pi_z)!=1)
      stop('Proportion of individuals that are randomized to the individual-level treatment arm must be a single number in (0,1)')
    if (correction==TRUE)
      message('No finite-sample correction will be done for the interaction test due to adequate degrees of freedom')

  } else if (test=="joint"){
    if (!is.numeric(delta_x) || delta_x == 0 || length(delta_x)!=1)
      stop('Effect size of the marginal cluster-level treatment effect must be a nonzero number')
    if (!is.numeric(delta_z) || delta_z == 0 || length(delta_z)!=1)
      stop('Effect size of the marginal individual-level treatment effect must be a nonzero number')
    if (!is.numeric(pi_x) || pi_x <= 0 || pi_x >= 1 || length(pi_x)!=1)
      stop('Proportion of clusters that are randomized to the cluster-level treatment arm must be a single number in (0,1)')
    if (!is.numeric(pi_z) || pi_z <= 0 || pi_z >= 1 || length(pi_z)!=1)
      stop('Proportion of individuals that are randomized to the individual-level treatment arm must be a single number in (0,1)')

  } else if (test=="I-U"){
    if (!is.numeric(delta_x) || delta_x == 0 || length(delta_x)!=1)
      stop('Effect size of the marginal cluster-level treatment effect must be a nonzero number')
    if (!is.numeric(delta_z) || delta_z == 0 || length(delta_z)!=1)
      stop('Effect size of the marginal individual-level treatment effect must be a nonzero number')
    if (!is.numeric(pi_x) || pi_x <= 0 || pi_x >= 1 || length(pi_x)!=1)
      stop('Proportion of clusters that are randomized to the cluster-level treatment arm must be a single number in (0,1)')
    if (!is.numeric(pi_z) || pi_z <= 0 || pi_z >= 1 || length(pi_z)!=1)
      stop('Proportion of individuals that are randomized to the individual-level treatment arm must be a single number in (0,1)')
  }

  if (!is.logical(correction))
    stop('Finite sample correction indicator should be a logical argument')


  if (!is.numeric(max_n) || max_n <= 0 || length(max_n)!=1)
    stop('Maximum number of clusters must be a positive number')

  if (!is.null(seed_mix)){
    if (!is.numeric(seed_mix) || length(seed_mix)!=1)
      stop('User-defined seed under the finite-sample corrected joint test must be numeric')
  }

  if (!is.numeric(size_mix) || size_mix <= 0 || length(size_mix)!=1)
    stop('Sample size for simulating the mix distribution under the finite-sample corrected joint test must be a positive number')

  if (!is.logical(verbose))
    stop('Message presentation indicator should be a logical argument')


  #Effect sizes might be negative
  delta_x <- abs(delta_x)
  delta_z <- abs(delta_z)
  delta_xz <- abs(delta_xz)


  a <- alpha
  b <- 1-power
  z_a <- qnorm(1-a/2)
  z_b <- qnorm(1-b)

  eta1 <- -m_bar*rho/(1+(m_bar-1)*rho)*(1-CV^2*m_bar*rho*(1-rho)/((1+(m_bar-1)*rho)^2))
  eta2 <- m_bar*(1-rho)/(1+(m_bar-1)*rho)*(1-CV^2*m_bar*rho*(1-rho)/((1+(m_bar-1)*rho)^2))-m_bar
  omega_x <- sigma2_y*(1-rho)/((m_bar+eta2)*pi_x*(1-pi_x))
  omega_z <- sigma2_y*(1-rho)/((m_bar+eta1)*pi_z*(1-pi_z))
  omega_xz <- sigma2_y*(1-rho)/((m_bar+eta1)*pi_x*(1-pi_x)*pi_z*(1-pi_z))


  #Check the omega to be positive
  if (omega_x<=0 || omega_z<=0 || omega_xz<=0)
    stop("Variance inflation factor is in an abnormal value. Please check whether the provided CV parameter is unrealistically large")


  #Re-iterate the given effect sizes and the chosen test
  if (verbose==TRUE){
    if (test=="cluster"){

      cat('Type of hypothesis test:\nTest for marginal cluster-level treatment effect')
      cat(paste0('\nEffect size:\n', delta_x, " for the marginal cluster-level treatment effect"))

    } else if (test=="individual"){

      cat('Type of hypothesis test:\nTest for marginal individual-level treatment effect')
      cat(paste0('\nEffect size:\n', delta_z, " for the marginal individual-level treatment effect"))

    } else if (test=="interaction"){

      cat('Type of hypothesis test:\nInteraction test')
      cat(paste0('\nEffect size:\n', delta_xz, " for the interaction effect"))

    } else if (test=="joint"){

      cat('Type of hypothesis test:\nJoint test')
      cat(paste0('\nEffect sizes:\n', delta_x, " for the marginal cluster-level treatment effect\n", delta_z, " for the marginal individual-level treatment effect"))

    } else if (test=="I-U"){

      cat('Type of hypothesis test:\nIntersection-union test')
      cat(paste0('\nEffect sizes:\n', delta_x, " for the marginal cluster-level treatment effect\n", delta_z, " for the marginal individual-level treatment effect"))

    }
  }


  pred.power <- function(n_input) {
    ### Test (A1): Test for marginal cluster-level treatment effect
    if (test=="cluster"){
      if (correction==FALSE){

        pred.power <- pnorm(sqrt(n_input*delta_x^2/omega_x)-z_a)
        if (verbose==TRUE){
          cat("\nA Wald z-test is used without finite-sample correction")
        }
      } else if (correction==TRUE){

        pred.power <- pt(qt(1-a/2, n_input-2), n_input-2, ncp=delta_x/sqrt(omega_x/n_input), lower.tail = F)
        + pt(qt(a/2, n_input-2), n_input-2, ncp=delta_x/sqrt(omega_x/n_input))
        if (verbose==TRUE){
          cat("\nA t-test is used for finite-sample correction")
        }
      }
    }

    ### Test (A2): Test for marginal individual-level treatment effect
    if (test=="individual"){

      pred.power <- pnorm(sqrt(n_input*delta_z^2/omega_z)-z_a)
      if (verbose==TRUE){
        cat("\nA Wald z-test is automatically used")
      }

    }

    ### Test (B): Interaction test
    if (test=="interaction"){

      pred.power <- pnorm(sqrt(n_input*delta_xz^2/omega_xz)-z_a)
      if (verbose==TRUE){
        cat("\nA Wald z-test is automatically used")
      }

    }

    ### Test (C): Joint test for the two marginal effects
    if (test=="joint"){
      if (correction==FALSE){

        theta <- n_input*(delta_x^2/omega_x + delta_z^2/omega_z)
        pred.power <- pchisq(qchisq(1-a, 2), 2, ncp = theta, lower.tail = F)
        if (verbose==TRUE){
          cat("\nA Chi-square test is used without finite-sample correction")
        }

      } else if (correction==TRUE){

        set.seed(seed_mix)
        #Simulate the mixed distribution (CENTRAL) to identify rejection region bound
        f.distn <- rf(size_mix, 1, n_input-2)
        chisq.distn <- rchisq(size_mix, 1)
        mix.distn <- f.distn + chisq.distn
        crt.value <- quantile(mix.distn, 1-a)

        #Simulate the mixed distribution (NONCENTRAL VERSION) for power calculation
        nc.f.distn <- rf(size_mix, 1, n_input-2, ncp = n_input*delta_x^2/omega_x)
        nc.chisq.distn <- rchisq(size_mix, 1, ncp = n_input*delta_z^2/omega_z)
        nc.mix.distn <- nc.f.distn + nc.chisq.distn

        pred.power <- mean(nc.mix.distn>crt.value)
        if (verbose==TRUE){
          cat("\nA simulation-based mixed F-Chi-square test is used for finite-sample correction")
        }

      }
    }

    ### Test (D): Intersection-union test for the two marginal effects
    if (test=="I-U"){
      if (correction==FALSE){

        wmean.c <- sqrt(n_input)*delta_x/sqrt(omega_x)
        wmean.i <- sqrt(n_input)*delta_z/sqrt(omega_z)
        pred.power <- pnorm(qnorm(1-a/2), mean = wmean.c, lower.tail = F)*pnorm(qnorm(1-a/2), mean = wmean.i, lower.tail = F)
        + pnorm(qnorm(1-a/2), mean = wmean.c, lower.tail = F)*pnorm(qnorm(a/2), mean = wmean.i)
        + pnorm(qnorm(a/2), mean = wmean.c)*pnorm(qnorm(1-a/2), mean = wmean.i, lower.tail = F)
        + pnorm(qnorm(a/2), mean = wmean.c)*pnorm(qnorm(a/2), mean = wmean.i)
        if (verbose==TRUE){
          cat("\nA z-based intersection-union test is used without finite-sample correction")
        }

      } else if (correction==TRUE){

        c.ncp <- sqrt(n_input)*delta_x/sqrt(omega_x)
        i.mean <- sqrt(n_input)*delta_z/sqrt(omega_z)
        pred.power <- pt(qt(1-a/2, n_input-2), n_input-2, c.ncp, lower.tail = F)*pnorm(qnorm(1-a/2), mean = i.mean, lower.tail = F)
        + pt(qt(1-a/2, n_input-2), n_input-2, c.ncp, lower.tail = F)*pnorm(qnorm(a/2), mean = i.mean)
        + pt(qt(a/2, n_input-2), n_input-2, c.ncp)*pnorm(qnorm(1-a/2), mean = i.mean, lower.tail = F)
        + pt(qt(a/2, n_input-2), n_input-2, c.ncp)*pnorm(qnorm(a/2), mean = i.mean)
        if (verbose==TRUE){
          cat("\nA mixed t- and z-based intersection-union test is used for finite-sample correction")
        }

      }
    }
    return(pred.power)
  }



  #Function to estimate number of clusters based on a required power level
  cluster.number <- function(power) {
    ### Test (A1): Test for marginal cluster-level treatment effect
    if (test=="cluster"){

      if (verbose==TRUE){
        if(correction==FALSE){
          cat("\nA Wald z-test is used without finite-sample-correction")
        } else if (correction==TRUE){
          cat("\nA t-test is used for finite-sample correction")
        }
      }
      n.out <- marginal.cluster(m_bar, CV, power, delta_x, rho, pi_x, correction, sigma2_y, a, max_n, z_a, z_b)

    }

    ### Test (A2): Test for marginal individual-level treatment effect
    if (test=="individual"){
      if (verbose==TRUE){
        cat("\nA Wald z-test is automatically used")
      }
      n.out <- marginal.ind(m_bar, CV, power, delta_z, rho, pi_z, sigma2_y, z_a, z_b)
    }

    ### Test (B): Interaction test
    if (test=="interaction"){
      if (verbose==TRUE){
        cat("\nA Wald z-test is automatically used")
      }
      n.out <- interaction(m_bar, CV, power, delta_xz, rho, pi_x, pi_z, sigma2_y, z_a, z_b)
    }

    ### Test (C): Joint test of the marginal effects
    if (test=="joint"){

      if (verbose==TRUE){
        if(correction==FALSE){
          cat("\nA Chi-square test is used without finite-sample correction.")
        } else if (correction==TRUE){
          cat("\nA simulation-based mixed F-Chi-square test is used for finite-sample correction")
        }
      }
      n.out <- joint(m_bar, CV, power, delta_x, delta_z, rho, pi_x, pi_z, correction, sigma2_y, a, max_n, seed_mix, size_mix)

    }

    ### Test (D): Intersection-Union test of the marginal effects
    if (test=="I-U"){

      if (verbose==TRUE){
        if(correction==FALSE){
          cat("\nA z-based intersection-union test is used without finite-sample correction")
        } else if (correction==TRUE){
          cat("\nA mixed t- and z-based intersection-union test is used for finite sample correction")
        }
      }
      n.out <- IU(m_bar, CV, power, delta_x, delta_z, rho, pi_x, pi_z, correction, sigma2_y, a, max_n)

    }

    return(n.out)
  }

  #When the user input n, the potentially inputted power will be ignored, and the formula will always give the predicted power;
  if (!is.null(n_input)){
    ans <- pred.power(n_input)
    if (verbose==TRUE){
      cat(paste0("\n\nPredicted power for the provided ", n_input, " clusters:\n"))
      cat(paste0(round(ans, 4), "\n"))
    }
  } else {#When the user does not input n, the formula will give the required number of clusters to hit the power
    ans <- cluster.number(power)
    if (verbose==TRUE){
      cat(paste0("\n\nRequired number of clusters to achieve ", power, " power:\n"))
      cat(paste0(ans, "\n"))
    }
  }

  return(ans)
}






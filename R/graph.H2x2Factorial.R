#' @title H2x2Factorial Plot
#'
#' @description The function \code{graph.H2x2Factorial} plots the sample size estimations or combinations of mean cluster sizes and cluster numbers
#' under variable CV for a chosen test. Based on the desired test and power, the function produces a plot with mean cluster size on the x-axis and number of clusters on
#' the y-axis, with multiple lines representing the dynamic sample size constraints if a vector of CV is specified. The limits of the y-axis
#' will be automatically adjusted based on the extreme values calculated. A color-blind-friendly palette is set by default but it can be updated by users.
#'
#' @usage
#' graph.H2x2Factorial(m_lower=10, m_upper=100, m_step=2,
#'                     CV=c(0,0.3,0.6,0.9),
#'                     palette=c("#0F2080","#85C0F9","#DDCC77","#F5793A","#A95AA1"),
#'                     line_width=rep(3,5), line_type=seq(1,5,1), title=NULL,
#'                     power=0.8, alpha=0.05,
#'                     pi_x=0.5, pi_z=0.5,
#'                     delta_x=0.25, delta_z=0.33, delta_xz=0.3, sigma2_y=1, rho=0,
#'                     test="cluster", correction=FALSE,
#'                     max_n=1e8, seed_mix=NULL, size_mix=1e4,
#'                     verbose=TRUE)
#'
#' @param m_lower a numeric value larger than 2 for the lower bound of the mean cluster sizes on the horizontal axis. Default is \code{10}.
#' @param m_upper a numeric value larger than \code{m_lower} for the upper bound of the mean cluster sizes on the horizontal axis. Default is \code{100}.
#' @param m_step a positive numeric value for the step size on the horizontal axis for plotting the sample size combinations. Default is \code{2}.
#' @param CV a vector of positive numeric values for a series of coefficients of variation of the cluster sizes. The length of CV vector equals the number
#' of lines presented in the plot, so the CV vector with a length less or equal to 5 is suggested for making a clear-looking graph. Besides, a reasonable magnitude of CV is highly recommended to produce effective plots.
#' Default is \code{c(0, 0.3, 0.6, 0.9)}.
#' @param palette a vector of character values to specify the color choices corresponding to the lines in the plot.
#' Default is \code{c("#0F2080", "#85C0F9", "#DDCC77", "#F5793A", "#A95AA1")}. The order should be matched with the specification of CV and the number of elements should be no less than that for CV vector.
#' @param line_width a vector of numeric values to specify the widths of the lines in the plot. Default is \code{rep(3, 5)}. The order should be matched with the specification of CV and the number of elements should be no less than that for CV vector.
#' @param line_type a vector of numeric values to specify the line types of the lines in the plot. Default is \code{seq(1, 5, 1)}. The order should be matched with the specification of CV and the number of elements should be no less than that for CV vector.
#' @param title a user-defined title or caption for the plot. Default is \code{NULL}. By default, a formal test name will be automatically given.
#' @param power a numeric value between 0 and 1 as the desired power level for sample size estimation. Default is \code{0.8}.
#' @param alpha a numeric value between 0 and 1 as the type I error rate. Default is \code{0.05}.
#' @param pi_x a numeric value between 0 and 1 as the proportion of clusters randomized to the cluster-level treatment. Default is \code{0.5}, representing a balanced allocation.
#' @param pi_z a numeric value between 0 and 1 as the proportion of individuals randomized to the individual-level treatment within each cluster. Default is \code{0.5}, representing a balanced allocation.
#' @param delta_x a nonzero numeric value for the (unstandardized) effect size of the marginal cluster-level treatment effect. Default is \code{0.25}, which is the hypothetical value for the example in the referenced paper.
#' @param delta_z a nonzero numeric value for the (unstandardized) effect size of the marginal individual-level treatment effect. Default is \code{0.33}, which is the hypothetical value for the example in the referenced paper.
#' @param delta_xz a nonzero numeric value for the (unstandardized) effect size of the interaction effect of the two treatments. Default is \code{0.3}, which is the hypothetical value for the example in the referenced paper.
#' @param sigma2_y a positive numeric value for the total variance of the continuous outcome. Default is \code{1}.
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
#' @return \code{graph.H2x2Factorial} returns a plot comparing the sample size requirements under different CV, with some suppressible messages.
#'
#' @export
#'
#' @examples
#' #Make a plot under the test for marginal cluster-level treatment effect
#' graph.H2x2Factorial(power=0.9, test="cluster", rho=0.1, verbose=FALSE)
#'
#' @importFrom stats qnorm pnorm qt pt qchisq pchisq rchisq qf pf rf quantile
#' @importFrom graphics plot lines mtext legend
#'
graph.H2x2Factorial <- function(m_lower=10,
                                m_upper=100,
                                m_step=2,
                                CV=c(0, 0.3, 0.6, 0.9),
                                palette=c("#0F2080", "#85C0F9", "#DDCC77", "#F5793A", "#A95AA1"),
                                line_width=rep(3, 5),
                                line_type=seq(1, 5, 1),
                                title=NULL,
                                power=0.8,
                                alpha=0.05, pi_x=0.5, pi_z=0.5,
                                delta_x=0.25, delta_z=0.33, delta_xz=0.3,
                                sigma2_y=1,
                                rho=0,
                                test="cluster",
                                correction=FALSE,
                                max_n=1e8,
                                seed_mix=NULL,
                                size_mix=1e4,
                                verbose=TRUE){

  if (!is.numeric(m_lower) || m_lower <=2 || length(m_lower)!=1)
    stop('The lower limit of mean cluster sizes must be a number greater than 2')

  if (!is.numeric(m_upper) || m_upper <= m_lower || length(m_upper)!=1)
    stop('The upper limit of mean cluster sizes must be numeric and greater than the lower limit')

  if (!is.numeric(m_step) || m_step <=0 || length(m_step)!=1)
    stop('The step size of mean cluster sizes for plotting must be a single positive number')

  for (i in 1:length(CV)){
    if (!is.numeric(CV[i]) || CV[i] < 0)
      stop('Each coefficient of variation of the cluster sizes must be a positive number')
  }
  if (length(unique(CV)) != length(CV)){
    CV <- unique(CV)
    CV.list <- CV[1]
    for (i in 2:length(CV)){
      CV.list <- paste0(CV.list, ",", CV[i])
    }
    warning(paste0("Duplicated elements of the CV input is deleted. CV input is changed to:\n(", CV.list, ")\n"))
  }
  for (i in 1:length(CV)){
    if (CV[i] < 0)
    stop('Each element of the CV vector must be a nonnegative number')
  }


  if (!is.character(palette))
    stop('Palette should be specified in the formal character format')
  if (length(palette)<length(CV))
    stop('The number of the specified colors should be no less than the number of CV provided')


  if (!is.numeric(line_width))
    stop('The line width parameter must be numeric')
  for (i in 1:length(line_width)){
    if (line_width[i] <= 0)
      stop('Each line width parameter must be a positive number')
  }
  if (length(line_width)<length(CV))
    stop('The number of specified line width should be no less than the number of CV')


  if (!is.numeric(line_type))
    stop('The line type parameter must be numeric')
  for (i in 1:length(line_type)){
    if (line_type[i] <= 0)
      stop('Each line type parameter must be a positive number')
  }
  if (length(line_type)<length(CV))
    stop('The number of specified line type should be no less than the number of CV')


  if (!is.null(title)){
    if (!is.character(title) || length(title)!=1){
      stop('The specified title must be a character string')
    }
  }


  if (!is.numeric(power) || power <= 0 || power >= 1 || length(power)!=1)
    stop('Target power must be a single number in (0,1)')

  if (!is.numeric(alpha) || alpha <= 0 || alpha >= 1 || length(alpha)!=1)
    stop('Type I error rate must be a single number in (0,1)')

  if (!is.numeric(sigma2_y) || sigma2_y <= 0 || length(alpha)!=1)
    stop('Total variance of the outcome must be a single positive number')

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


  #Re-iterate the given effect sizes, the chosen test, and the theoretical test name
  if (verbose==TRUE){
    if (test=="cluster"){

      cat('Type of hypothesis test:\nTest for marginal cluster-level treatment effect')
      cat(paste0('\nEffect size:\n', delta_x, " for the marginal cluster-level treatment effect"))
      if (correction==FALSE){
        cat("\nA Wald z-test is used without finite-sample correction\n")
      } else if (correction==TRUE){
        cat("\nA t-test is used for finite-sample correction\n")
      }

    } else if (test=="individual"){

      cat('Type of hypothesis test:\nTest for marginal individual-level treatment effect')
      cat(paste0('\nEffect size:\n', delta_z, " for the marginal individual-level treatment effect"))

    } else if (test=="interaction"){

      cat('Type of hypothesis test:\nInteraction test')
      cat(paste0('\nEffect size:\n', delta_xz, " for the interaction effect"))

    } else if (test=="joint"){

      cat('Type of hypothesis test:\nJoint test')
      cat(paste0('\nEffect sizes:\n', delta_x, " for the marginal cluster-level treatment effect\n", delta_z, " for the marginal individual-level treatment effect"))
      if (correction==FALSE){
        cat("\nA Chi-square test is used without finite-sample correction\n")
      } else if (correction==TRUE){
        cat("\nA simulation-based mixed F-Chi-square test is used for finite-sample correction\n")
      }

    } else if (test=="I-U"){

      cat('Type of hypothesis test:\nIntersection-union test')
      cat(paste0('\nEffect sizes:\n', delta_x, " for the marginal cluster-level treatment effect\n", delta_z, " for the marginal individual-level treatment effect"))
      if (correction==FALSE){
        cat("\nA z-based intersection-union test is used without finite-sample correction\n")
      } else if (correction==TRUE){
        cat("\nA mixed t- and z-based intersection-union test is used for finite-sample correction\n")
      }

    }
  }


  if (is.null(title)){
    if (test=="cluster"){
      if (correction==TRUE){title <- "Test for marginal cluster-level treatment effect \n with finite-sample correction"}
      else if (correction==FALSE){title <- "Test for marginal cluster-level treatment effect \n without finite-sample correction"}

    } else if (test=="individual"){
      title <- "Test for marginal individual-level treatment effect"

    } else if (test=="interaction"){
      title <- "Interaction test"

    } else if (test=="joint"){
      if (correction==TRUE){title <- "Joint test with finite-sample correction"}
      else if (correction==FALSE){title <- "Joint test without finite-sample correction"}

    } else if (test=="I-U"){
      if (correction==TRUE){title <- "Intersection-union test \n with finite-sample correction"}
      else if (correction==FALSE){title <- "Intersection-union test \n without finite-sample correction"}

    }
  }

  plot(Mset, line.CV[[1]], ylim=c(y_lower, y_upper), las=1,
       xlab=expression(bar(m)),
       main=title,
       ylab="", cex.lab=1, cex.axis=1, cex.main=1, cex=1,
       type="l", lwd=line_width[1], col=palette[1], lty=line_type[1])
  mtext(expression(n),side=2,las=1,line=3,cex=1.2)

  if (length(CV)>1){
    for (k in 2:length(CV)){
    lines(Mset, line.CV[[k]], type="l", lwd=line_width[k], col=palette[k], lty=line_type[k])
    }
  }

  legend_content <- NULL
  for (l in 1:length(CV)){
    legend_content <- c(legend_content, paste0("CV = ", CV[l]))
  }

  legend("topright", inset=0.01, legend=legend_content,
         col=palette[1:length(CV)], lty=rep(1,length(CV)), cex=0.8, lwd=3, box.lty=0)

}



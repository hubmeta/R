#'Morris weigth function
#'
#' @param correlations Correlations to meta-analyze.
#' @param sample_sizes Sample sizes aligned to `correlations`.
#' @param reliability_of_x Reliability estimates for the x variable in each
#'   study.
#' @param reliability_of_y Reliability estimates for the y variable in each
#'   study.
#' @param significance_levels Numeric vector giving interval levels, for example
#'   `c(0.95, 0.80)`.
#' @param default_reliability_of_x Default replacement for missing x
#'   reliabilities. If `NULL`, the mean of the observed reliabilities is used.
#' @param default_reliability_of_y Default replacement for missing y
#'   reliabilities. If `NULL`, the mean of the observed reliabilities is used.
#'
#' @return A one-row data frame of Morris-weight summary results.
#'@export
#'

morris_weight_analysis <- function(
  correlations,
  sample_sizes,
  reliability_of_x,
  reliability_of_y,
  significance_levels,
  default_reliability_of_x=1,
  default_reliability_of_y=1
) {
  #INPUTS:
  #Am <- c()                                                       #TODO: For future purposes to included other attenuation factors
  if(is.null(default_reliability_of_x)){
    default_reliability_of_x <- mean(reliability_of_x, na.rm=TRUE)

  }

  if(is.null(default_reliability_of_y)){
    default_reliability_of_y <- mean(reliability_of_y, na.rm=TRUE)
  }
  if(is.nan(default_reliability_of_x)){
    default_reliability_of_x <- 1
  }

  if(is.nan(default_reliability_of_y)){
    default_reliability_of_y <- 1
  }

  reliability_of_x <- replace(reliability_of_x, is.na(reliability_of_x), default_reliability_of_x)
  reliability_of_y <- replace(reliability_of_y, is.na(reliability_of_y), default_reliability_of_y)

  data <- cbind(correlations, sample_sizes, reliability_of_x, reliability_of_y)                                    #Pack data in dataframe, just for checking purposes

  K <- length(sample_sizes)                                                    #Number of studies
  N <- sum(sample_sizes)                                                    #Total sample size

  if(K<2){
    results <- cbind(K, N, round(correlations[1],2), NA, NA)

    colnames(results) <- c("K", "N", "rcmean", "CRlowrc", "CRhighrc")

    results <- as.data.frame(results)

    return(results)
  }

  data <-cbind(correlations, sample_sizes, reliability_of_x, reliability_of_y) 							   		#Pack data in dataframe, just for checking purposes

  K <- length(sample_sizes)                    		 					    #Number of studies
  N <- sum(sample_sizes)                         							#Total sample size

  #Calculating Morris Weights


  # We can't have NA for reliability, so we need to replace it with 1s
  reliability_of_x[is.na(default_reliability_of_x)] <- 1.00000000
  reliability_of_y[is.na(default_reliability_of_y)] <- 1.00000000
  ri <- correlations
  ni <- sample_sizes


  ########################################################################
  #Adding the N - 1 adjustment because it is more accurate
  #Changing credibility intervals to 95%

  ri <- ifelse(ri == 0, 0.0000001, ri) # avoid dividing by zero later on
  M1 <- sum((ni)*ri)/sum(ni) # weighted mean, uncorrected (bare bones)
  # M1 This is what HS produce
  ni <- (ni - 1) # adjusted to be more accurate

  var.i <-((1-M1^2)^2)/(ni) # individual study sampling variance (bare bones) N - 1 otherwise
  ai <- sqrt(reliability_of_x*reliability_of_y) # individual study attenuation for unreliability
  # ai <- 1 If you want to have no reliability corrections
  var.i2 <- var.i/ai^2 # corrected individual study sampling variance for #unreliability

  ri.m <- ri/ai # corrected individual study ES for attenutation due to unreliablity and range restriction
  ri.m <- ifelse(ri.m > 1, 1, ri.m) # avoid over corrections
  ri.m <- ifelse(ri.m < -1, -1, ri.m) # avoid over corrections

  morris.dat <- data.frame(cbind(ri.m,var.i2, ni)) # collect corrected estimates of rxy and error

  result = tryCatch({
    morris1 <- metafor::rma(yi=ri.m,vi=var.i2,data=morris.dat,
                            control=list(maxiter=1000, stepadj=.5)) # run the random-effects meta with REML
    morris1 # print the result
    Morris.M.rho <- morris1$b
    Morris.V.rho <- morris1$tau2 # random-effects variance component
    Morris.SD.rho <- sqrt(morris1$tau2)
    Morris.CR95.L <- (Morris.M.rho - 1.96*Morris.SD.rho) # Lower CR Bound
    Morris.CR95.U <- (Morris.M.rho + 1.96*Morris.SD.rho) # Upper CR Bound
    Morris.CR90.L <- (Morris.M.rho -  1.645*Morris.SD.rho) # Lower CR Bound
    Morris.CR90.U <- (Morris.M.rho +  1.645*Morris.SD.rho) # Upper CR Bound

    round(Morris.CR95.L,2)
    round(Morris.CR95.U,2)
    round(Morris.CR90.L,3)
    round(Morris.CR90.U,3)

    options(max.print=10000)
    inf <- stats::influence(morris1)
    Morris.DFFITS <- inf$inf$dffits
    Morris.Outlier <- inf$inf$inf
    Morris.rcmean  <- morris1$b


    if(Morris.rcmean[1][1]>1){
      Morris.rcmean[1][1] <- 1
    }
    if(-1>Morris.rcmean[1][1]){
      Morris.rcmean[1][1] <- -1
    }

    results <- cbind(K, N, Morris.rcmean[1][1], round(Morris.CR90.L,3)[1][1], round(Morris.CR90.U,3)[1][1])

    colnames(results) <- c("K", "N", "rcmean", "CRlowrc", "CRhighrc")

    results <- as.data.frame(results)

    return(results)
  },error = function(error_condition) {
    results <- cbind(K, N, NA, NA, NA)

    colnames(results) <- c("K", "N", "rcmean", "CRlowrc", "CRhighrc")

    results <- as.data.frame(results)

    return(results)
  })

}

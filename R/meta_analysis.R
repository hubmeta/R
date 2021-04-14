#'Square the file
#'
#'@param correlations: correlations should pass like c(.18, .0, .08, .15, .27, .1, .28, .17, .02, .28)
#'@param sample_sizes: sample sizes should pass like c(426, 328, 122, 284, 472, 154, 372, 674, 110, 116)
#'@param reliability_of_x: reliability of x should pass like c(.85, NA, NA, .86, .80, .79, .91, .85, .92, .85)
#'@param reliability_of_y: reliability of y should pass like c(.63, .63, .62, .39, .24, .85, .89, .48, .68, .84)
#'@param significance_levels: Significance level of 1)Confidence intervals; 2)Credibility intervals and it should pass
#'@param default_reliability_of_x: if you add any number it will be the default value and if you pass NULL it will replace with average of other reliabilities
#'#'@param default_reliability_of_y: if you add any number it will be the default value and if you pass NULL it will replace with average of other reliabilities

#' like c(0.95, 0.80)
#'@return result data frame
#'@export
#'

meta_analysis <- function(
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

  Am <- rep(1, K)                                                 # THIS WILL CHANGE IF Am is entered by user

  correlations <- ifelse(correlations == 0, 0.0000001, correlations)                            #Replace zero values to avoid division by zero
  Af <- reliability_of_x * reliability_of_y * Am                                                #The final attenuation factor
  Af <- ifelse(Af == 0, 0.0000001, Af)                            #Replace zero values to avoid division by zero

  W <- sample_sizes * Af                                                        #Adjusted weights, similar to column V
  Ar <- correlations / sqrt(Af)                                                #Adjusted correlations

  # Calculating un-adjusted for reliability values

  rmean <- sum(sample_sizes * correlations) / N                                             #Mean Effect Size (r ), H-S method unadjusted


  sigmar2 <- sum(sample_sizes * (correlations - rmean)^2) / N                                 #Frequency-weighted average squared error (ơr^2)
  SDr <- sqrt(sigmar2)                                        #Frequency-weighted Standard dev. (SDr )
  SEr <- SDr / sqrt(K)                                            #Standard error of mean r (SEr )

  sigmae2 <- (1 - rmean^2)^2 / (N / K - 1)                            #Sampling Error Variance (estimate)
  sigmarho2 <- ifelse((sigmar2 - sigmae2) < 0, 0, sigmar2 - sigmae2)  # Residual Variance (variance of population) (ơp^2)
  sigmarho <- sqrt(sigmarho2)                                    #Residual S.D. (ơp)
  PercExp <- sigmae2 / sigmar2                                    #Percentage explained (by sampling error)

  siglev <- qnorm(1 - (1 - significance_levels[1]) / 2)                              #Significance level confidence interval
  credsig <- qnorm(1 - (1 - significance_levels[2]) / 2)                             #Significance level credibility interval

  CIlowr <- rmean - siglev * (1 - rmean^2) / sqrt(N - K)                                    #Random-Effects Model, Un-Adjusted , Confidence interval lower bound
  CIhighr <- rmean + siglev * (1 - rmean^2) / sqrt(N - K)                                #Random-Effects Model, Un-Adjusted , Confidence interval upper bound

  CRlowr <- rmean - credsig * sigmarho                                #Random-Effects Model, Un-Adjusted , Credibility interval lower bound
  CRhighr <- rmean + credsig * sigmarho                            ##Random-Effects Model, Un-Adjusted , Credibility interval upper bound

  # Calculating adjusted for reliability values

  rcmean <- sum(W * Ar) / sum(W)                                      #Mean Effect Size (r ), H-S method adjusted

  VarRc <- sum(W * (Ar - rcmean)^2) / sum(W)                            #Var(rc)  Based on Hunter & Schmidt (2004), page 125-126

  AveVe <- sum(W * ((1 - rmean^2)^2 / (sample_sizes - 1) / Af)) / sum(W)            #Ave(Ve)
  VarRho <- ifelse((VarRc - AveVe) < 0, 0, VarRc - AveVe)        #Var (ρ)
  SDrho <- sqrt(VarRho)                                            #SDρ

  SDrc <- sqrt(VarRc)                                            #Adjusted frequency-weighted Standard dev.
  SErc <- SDrc / sqrt(K)                                        #Standard error of mean adjusted r

  CIlowrc <- rcmean - siglev * SErc                                    #Random-Effects Model, Adjusted , Confidence interval lower bound
  CIhighrc <- rcmean + siglev * SErc                                #Random-Effects Model, Adjusted , Confidence interval upper bound

  CRlowrc <- rcmean - credsig * SDrho                                #Random-Effects Model, Adjusted , Credibility interval lower bound
  CRhighrc <- rcmean + credsig * SDrho                                #Random-Effects Model, Adjusted , Credibility interval upper bound

  #Calculating Fisher's Z model estimate, Fixed effects model

  FZ <- 0.5 * log((1 + correlations) / (1 - correlations))                                    #Fisher's Z transform of correlations (effect sizes)
  SEF <- ifelse(sample_sizes < 4, 0.00001, 1 / sqrt(sample_sizes - 3))                #Standard error of Fisher's Z
  AZ <- FZ / sqrt(Af)
  WF <- Af / SEF^2                                                #Weights of fixed-effects Fisher's Z estimates

  AveZ <- sum(WF * AZ) / sum(WF)                                        #Fixed-effects adjusted average weighted effect size, Fisher's Z
  rfisherfixed <- (exp(2 * AveZ) - 1) / (exp(2 * AveZ) + 1)       #Mean Effect Size (r) , adjusted, Fixed-effects Fisher's Z model
  SEZes = sqrt(1 / sum(W))                                            #Standard Error (Z)

  Zc <- sum(FZ / SEF) / sqrt(K)
  FSN <- K * ((Zc / siglev)^2 - 1)                          #Failsafe N Rosenthal(1979)

  CIESlow <- AveZ - siglev * SEZes                                    #Fixed-Effects Model, Adjusted , Confidence interval lower bound Fisher's Z estimate
  CIEShigh <- AveZ + siglev * SEZes                                    #Fixed-Effects Model, Adjusted , Confidence interval upper bound Fisher's Z estimate

  FisherCILow <- (exp(2 * CIESlow) - 1) / (exp(2 * CIESlow) + 1)      #Correlation equivalent of Fixed-Effects Model, Adjusted , Confidence interval lower bound Fisher's Z estimate
  FisherCIHigh <- (exp(2 * CIEShigh) - 1) / (exp(2 * CIEShigh) + 1)   #Correlation equivalent of Fixed-Effects Model, Adjusted , Confidence interval upper bound Fisher's Z estimate

  Q <- sum(W * AZ^2) - sum(W * AZ)^2 / sum(W)                            #Hedges Q test of homogeniety
  dr <- K - 1                                                        #Degrees of Freedom

  I2 <- ifelse(Q <= K - 1, 0, (Q - (K - 1)) / Q)                        #I^2 index (Magnitude of Heterogeneity)
  tau2 <- ifelse(Q <= K - 1, 0, (Q - (K - 1)) / (sum(W) - sum(W^2) / sum(W)))    #Population Variability in Effect Sizes

  #Calculating Fisher's Z model estimate, Random effects model

  FisherW <- 1 / (tau2 + SEF^2 / Af)                                #Weights of random-effects Fisher's Z estimates

  REAveZ <- sum(FisherW * AZ) / sum(FisherW)                        #Random-effects average weighted effect size, Fisher's Z
  rfisherrandom <- (exp(2 * REAveZ) - 1) / (exp(2 * REAveZ) + 1)  #Mean Effect Size (r) , adjusted, random-effects Fisher's Z model
  RESEz = 1 / sqrt(sum(FisherW))                                    #Random-effects standard Error (Z)

  RECIESlow <- REAveZ - siglev * RESEz                                #Random-Effects Model, Adjusted , Confidence interval lower bound Fisher's Z estimate
  RECIEShigh <- REAveZ + siglev * RESEz                            #Random-Effects Model, Adjusted , Confidence interval upper bound Fisher's Z estimate

  REFisherCILow <- (exp(2 * RECIESlow) - 1) / (exp(2 * RECIESlow) + 1)      #Correlation equivalent of random-Effects Model, Adjusted , Confidence interval lower bound Fisher's Z estimate
  REFisherCIHigh <- (exp(2 * RECIEShigh) - 1) / (exp(2 * RECIEShigh) + 1)   #Correlation equivalent of random-Effects Model, Adjusted , Confidence interval lower bound Fisher's Z estimate

  if(rmean>1){
    rmean <- 1
  }
  if(-1>rmean){
    rmean <- -1
  }
  if(rcmean>1){
    rcmean <- 1
  }
  if(-1>rcmean){
    rcmean <- -1
  }

  results <- cbind(K, N, rmean, sigmar2, SDr, SEr, sigmae2, sigmarho2, sigmarho, PercExp, CIlowr, CIhighr,
                   CRlowr, CRhighr, rcmean, VarRc, AveVe, VarRho, SDrho, SDrc, SErc, CIlowrc, CIhighrc,
                   CRlowrc, CRhighrc, AveZ, rfisherfixed, SEZes, CIESlow, CIEShigh, FisherCILow, FisherCIHigh, Q, dr,
                   I2, tau2, REAveZ, rfisherrandom, RESEz, RECIESlow, RECIEShigh, REFisherCILow, REFisherCIHigh, FSN)

  colnames(results) <- c("K", "N", "rmean", "sigmar2", "SDr", "SEr", "sigmae2", "sigmarho2", "sigmarho", "PercExp", "CIlowr", "CIhighr",
                         "CRlowr", "CRhighr", "rcmean", "VarRc", "AveVe", "VarRho", "SDrho", "SDrc", "SErc", "CIlowrc", "CIhighrc",
                         "CRlowrc", "CRhighrc", "AveZ", "rfisherfixed", "SEZes", "CIESlow", "CIEShigh", "FisherCILow", "FisherCIHigh", "Q", "dr",
                         "I2", "tau2", "REAveZ", "rfisherrandom", "RESEz", "RECIESlow", "RECIEShigh", "REFisherCILow", "REFisherCIHigh", "FSN")

  results <- as.data.frame(results)

  return(results)
}


#' Test of RMST for comparing two or more survival curves
#'
#' Performs the test of Restricted Mean Survival Time for two or more survival curves
#' by comparing the difference of areas under survival curves.
#'
#' For each group, the Restricted Mean Survival Time at time tau (RMST(tau)) is computed as the area
#' under the survival curve between time 0 and tau. The test of RMST is a pairwise
#' multiple comparison test. For each pair of groups, it tests whether the difference
#' between the RMST(tau) is zero or not. If the difference is not null, then the
#' survival curves cannot be equal.
#'
#' For exactly two groups, a single test is performed. For more than two survival curves,
#' it compares each survival curve to every other curves and tests the global null
#' hypothesis "all curves are equal" against the hypothesis "the curves are not all equal".
#' 
#' @usage multi_rmst(df, tau, method = p.adjust.methods, nboot = 500)
#'
#' @param df A dataframe with columns :
#'   * `time` : positive numbers, time-to-event;
#'   * `status` : vector of integer, 0 or 1. 0 is (right) censoring, 1 is event;
#'   * `arm` : a factor or object that can be coerced to one. The group the patient 
#'     belongs to. Must have at least two levels.
#' @param tau The truncation time, default is the lowest of the max(time) of the groups;
#' @param nboot Number of bootstrap samples;
#' @param method The correction used for the p-values. Must be in [p.adjust.methods]. 
#'   Default is the Holm correction. Unused if number of groups equals two.
#'
#' @return An object of class `multi_rmst` containing :
#'   * `rmst_mat` RMST estimation for each arm;
#'   * `results` A matrix. Each row represents a comparison of two curves and contains the difference
#'     of RMST, its standard deviation, the p-value and the adjusted p-value;
#'   * `p` The p-value of the global test;
#'   * `nb_tests` The number of performed tests;
#'   * The parameters `tau`, `method` and `nboot`.
#'
#' @references Royston, P., & Parmar, M. K. (2013). Restricted mean survival time:
#'     an alternative to the hazard ratio for the design and analysis of randomized
#'     trials with a time-to-event outcome. BMC medical research methodology, 13, 1-15.
#'
#' @examples
#'   multi_rmst(data_under_PH, tau = 36, nboot = 300)
#'   multi_rmst(data_not_PH, tau = 36, method = "BH", nboot = 300)
#'
#' @export
multi_rmst = function(df, tau, method = p.adjust.methods, nboot = 500){
  method = match.arg(method)
  
  if (!all(c("time", "status", "arm") %in% colnames(df))){
    stop("The dataframe must contain the columns 'time', 'status' and 'arm'.")
  }
  
  if (!all(df$status %in% c(0,1))){stop("'status' must be either 0 or 1.")}
  
  df$arm = as.factor(df$arm)
  lev = levels(df$arm)
  df$arm = as.numeric(df$arm)
  nb_arms = length(lev)
  if (nb_arms < 2){stop("Need at least two groups.")}
  
  if (missing(tau)){tau = min(tapply(X=df$time, INDEX=df$arm, FUN=max))}
  
  rmst_mat = matrix(NA, nrow = nb_arms, ncol = 2)
  label_rmst = rep(NA, nb_arms)
  for (i in 1:nb_arms){
    label_rmst[i] = paste("arm", lev[i])
    df_i = df[df$arm == i,]
    X = boot(df_i, rmst_unique, R = nboot, tau = tau)
    rmst_mat[i, 1] = X$t0
    rmst_mat[i, 2] = sqrt(var(X$t))
  }
  
  nb_tests = nb_arms*(nb_arms - 1) / 2
  label_test = rep(NA, nb_tests)
  results = matrix(NA, nrow = nb_tests, ncol = 4)
  k = 1
  for (i in 1:(nb_arms - 1)){
    for (j in (i+1):nb_arms){
      label_test[k] = paste(lev[i], "VS", lev [j])
      results[k, 1] = rmst_mat[j, 1] - rmst_mat[i, 1]
      results[k, 2] = sqrt(rmst_mat[j, 2]^2 + rmst_mat[i, 2]^2)
      results[k, 3] = 2 * (1 - stats::pnorm(abs(results[k, 1] / results[k, 2])))
      
      k=k + 1
    }
  }
  
  rownames(rmst_mat) = label_rmst
  colnames(rmst_mat) = c("rmst", "sd")
  
  if (nb_tests == 1){
    results = results[, -4, drop=FALSE]
    colnames(results) = c("dRMST","sd","p-value")
    p = results[1,3]
  }
  else {
    p_adjusted = p.adjust(results[,3], method=method)
    results[,4] = p_adjusted
    p = min(p_adjusted)
    rownames(results) = label_test
    colnames(results) = c("dRMST","sd","p","p adjusted")
  }
  
  z = list(results = results, rmst_mat = rmst_mat, tau = tau, p = p,
           nb_tests = nb_tests, method = method, nboot = nboot)
  class(z) = "multi_rmst"
  return(z)
}




#' Print method for the multiple test of RMST
#'
#' @param x An object of class `multi_rmst` as returned by [multi_rmst()];
#' @param ... For compatibility with the `print` method, unused and to be ignored.
#' 
#' @return None
#' 
#' @examples 
#'   x = multi_rmst(data_under_PH, tau = 36, nboot = 300)
#'   print(x)
#'
#' @export
print.multi_rmst = function(x, ...){
  nb_tests=x$nb_tests
  
  cat("(Multiple) test of RMST \n")
  cat("Truncation time :",x$tau," \n")
  if (nb_tests > 1) {cat("Correction :",x$method,"\n\n")}
  cat("RMST estimation for each arm \n")
  print(x$rmst_mat)
  cat("\nPair-wise comparisons \n")
  print(x$results)
  cat("",end="\n")
  if (nb_tests > 1) {cat("p=",x$p,sep="")}
}




rmst_unique = function(df, ind, tau){
  if (missing(ind)){ind = 1:length(df$time)}
  time = df$time[ind]
  status = df$status[ind]
  arm = df$arm[ind] # shuffling the arm variable
  
  KM = survival::survfit(survival::Surv(time, status)~1)
  S = stats::stepfun(KM$time, c(1, KM$surv), right = FALSE)
  time_evt = unique(sort(c(0, KM$time[KM$time < tau], tau)))
  time_int = time_evt[-1]-time_evt[-length(time_evt)]
  surv = S(time_evt[-length(time_evt)])
  rmst = sum(surv * time_int)
  
  return(rmst)
}


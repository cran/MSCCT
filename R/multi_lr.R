#' (Weighted) Log-rank test for comparison of two or more survival curves.
#'
#' Performs a global log-rank test for comparing two or more survival curves.
#'
#' `weights` contains the chosen weights for the test. It must be a vector, a matrix
#' or an object that can be coerced to a matrix, like a data frame (passed as 
#' argument to `as.matrix`). Can be omitted.
#' 
#' If not given (default), then perform either a log-rank test, a Gehan-Wilcoxon
#' test or a Fleming-Harrington test depending on the choice of `test`.
#' 
#' If `weights` is a one-dimension vector, its length must be equal to the
#' number of distinct time of event and throws an error if it is not true. In this
#' case, `multi_lr()` performs a weighted log-rank test with the specified weights.
#' 
#' If `weights` is a matrix (or a two-dimension object), its number of rows
#' must be equal to the number of distinct time of event and throws an error if
#' it is not true. In this case, `multi_lr()` performs as many tests as the
#' number of columns in `weights`. The first test is a weighted log-rank test
#' with weights the first column of `weights`, the second test is a weighted
#' log-rank test with weights the second column of `weights`, and so on.
#' 
#' @usage multi_lr(df, weights, test = c("lr", "gw", "fh"), rho = 1, gamma = 0)
#'
#' @param df A data frame with columns :
#'   * `time` : positive numbers, time-to-event;
#'   * `status` : vector of integer, 0 or 1. 0 is (right) censoring, 1 is event;
#'   * `arm` : a factor or object that can be coerced to one. The group the patient 
#'     belongs to. Must have at least two levels.
#' @param weights An object that can be coerced to a matrix. The weights
#'   used for the tests. Can be omitted (see Details);
#' @param test If `weights` is omitted, specify the test to perform. Possible values are
#'   `lr` for log-rank, `gw` for Gehan-Wilcoxon, and `fh` for Flemming-Harrington;
#' @param rho,gamma The parameters for Flemming-Harrington test. Default is (rho,gamma)=(1,0),
#'   which is also called the Peto-Peto test.
#'
#' @return An object of class `multi_lr` containing:
#'   * `U` : Statistics of tests;
#'   * `p` : The corresponding p-values;
#'   * `degree` : Degrees of freedom of the statistics of tests;
#'   * The argument `test`, changed to "chosen" if weights are given.
#'
#' @export
#'
#' @examples
#'   # Log-rank test
#'   multi_lr(data_not_PH)
#'
#'   # Gehan-Wilcoxon test
#'   multi_lr(data_not_PH, test="gw")
#'
#'   # It is possible to run several tests with different weights at a time
#'   evt_time = unique(data_not_PH$time[data_not_PH$status == 1])
#'   nb_evt_time = length(evt_time)
#'   weights = matrix(runif(nb_evt_time*3), ncol=3)
#'   multi_lr(data_not_PH, weights=weights)
multi_lr = function(df, weights, test = c("lr", "gw", "fh"), rho = 1, gamma = 0){
  test = match.arg(test)
  
  if (!all(c("time", "status", "arm") %in% colnames(df))){
    stop("The dataframe must contain the columns 'time', 'status' and 'arm'.")
  }
  
  if (!all(df$status %in% c(0,1))){stop("'status' must be either 0 or 1.")}
  
  df$arm = as.factor(df$arm)
  nb_arms = nlevels(df$arm)
  if (nb_arms < 2){stop("Needs at least two groups.")}
  
  time = df$time
  status = df$status
  arm = df$arm
  
  evt_time = unique(time[status == 1])
  evt_time_ordered = sort(evt_time)
  
  Y = matrix(0, nrow=nb_arms, ncol=length(evt_time_ordered))
  D = matrix(0, nrow=nb_arms, ncol=length(evt_time_ordered))
  for (i in 1:length(evt_time_ordered)){
    t = evt_time_ordered[i]
    index = time >= t
    M = table(arm[index], time[index] == t)
    D[,i] = M[,"TRUE"]
    Y[,i] = apply(M, 1, sum)
  }
  d_total = apply(D, 2, sum)
  y_total = apply(Y, 2, sum)
  D = D[1:(nb_arms-1), , drop=FALSE]
  Y = Y[1:(nb_arms-1), , drop=FALSE]
  
  
  if (missing(weights)) {
    if (test == "lr"){weights = matrix(1, nrow=length(evt_time_ordered), ncol=1)}
    else if (test == "gw"){weights = matrix(y_total, ncol=1)}
    else if (test == "fh"){
      KM = survival::survfit(Surv(time,status)~1)
      S = stats::stepfun(KM$time, c(1,KM$surv), right=FALSE)
      weights = matrix(S(evt_time_ordered)^rho * (1-S(evt_time_ordered))^gamma, ncol=1)
    }
    else {stop("Incorrect argument for 'test'")}
  } else {
    weights = as.matrix(weights)
    if (nrow(weights) != length(evt_time_ordered)) {
      stop(sprintf("Expected %d weights, got %d", length(evt_time_ordered), nrow(weights)))
    }
    test = "chosen"
  }
  
  E = Y %*% diag(d_total/y_total)
  V = (D-E) %*% weights
  
  sigma = as.list(rep(0,length.out=ncol(weights)))
  coef = t(weights^2) %*% diag(d_total * (y_total-d_total) / (y_total^2 * (y_total-1)))
  for (i in 1:length(evt_time_ordered)){
    Y_i = matrix(Y[,i], ncol=1)
    M_i = y_total[i] * diag(Y[,i], nrow=length(Y[,i])) - Y_i %*% t(Y_i)
    M_ji = lapply(coef[,i], "*", M_i)
    sigma = lapply(seq_along(sigma), function(j){sigma[[j]] + M_ji[[j]]})
  }
  
  U = lapply(seq_along(sigma), function(j){t(V[,j]) %*% solve(sigma[[j]], V[,j])})
  U = as.numeric(U)
  p = 1 - pchisq(U,nb_arms-1)
  
  z = list(U=U, p=p, degree=nb_arms-1, test=test, rho=rho, gamma=gamma)
  class(z) = "multi_lr"
  return(z)
}

#' Print method for the multiple log-rank test
#'
#' @param x An object of class `multi_lr` as returned by [multi_lr()];
#' @param ... For compatibility with the `print` method, unused and to be ignored.
#' 
#' @return None
#' 
#' @examples
#'   x = multi_lr(data_not_PH)
#'   print(x)
#' 
#'
#' @export
print.multi_lr = function(x, ...){
  cat("(Multiple) Weighted log-rank test \n\n")
  cat("Weighting : ")
  if (x$test == "chosen"){cat("Specified weighted \n")}
  else if (x$test == "lr"){cat("Classic log-rank test \n")}
  else if (x$test == "gw"){cat("Gehan-Wilcoxon test \n")}
  else if (x$test == "fh"){
    cat("Flemming-Harrington test \n")
    cat("Parameters : rho =", x$rho, ", gamma = ", x$gamma, "\n")}
  cat("Degrees of freedom :", x$degree, "\n\n")
  M = matrix(c(x$U, x$p), ncol=2, byrow=FALSE)
  labs = rep(NA, length(x$U))
  for (i in 1:length(labs)){labs[i] = paste("Test", i)}
  rownames(M) = labs
  colnames(M) = c("Statistic", "p")
  print(M)
}



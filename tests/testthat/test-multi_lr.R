test_that("The given data frame is correctly shaped", {
  df = data_under_PH
  df["random"] = runif(nrow(df))
  expect_error(multi_lr(df[c("time","status","random")]),
               "The dataframe must contain the columns 'time', 'status' and 'arm'.")
})


test_that("The arguments given are correct", {
  # Stop if status is not 0 or 1
  df = data_under_PH
  df$status[df$status == 0] = -1
  expect_error(multi_lr(df), "'status' must be either 0 or 1.")
  
  
  # Should work if arm and status are factor instead of integers
  df = data_under_PH
  df$status = as.factor(df$status)
  df$arm = as.factor(df$arm)
  expect_no_error(multi_lr(df))
})


test_that("The output is correct", {
  # log-rank with two curves
  index = (data_not_PH$arm == 1) | (data_not_PH$arm == 0)
  df = data_not_PH[index,]
  X = multi_lr(df)
  Y = survival::survdiff(survival::Surv(time,status)~arm, data=df)
  expect_equal(X$U, Y$chisq)
  expect_equal(X$p, Y$pvalue)
  
  # Flemming-Harrington avec rho=c(0,1,2), gamma=0 and three curves
  X = matrix(NA, nrow=3, ncol=2)
  Y = matrix(NA, nrow=3, ncol=2)
  for (rho in 0:2){
    testX = multi_lr(data_not_PH, test="fh", rho=rho, gamma=0)
    testY = survival::survdiff(survival::Surv(time,status)~arm, data=data_not_PH, rho=rho)
    X[rho+1,1] = testX$U
    X[rho+1,2] = testX$p
    Y[rho+1,1] = testY$chisq
    Y[rho+1,2] = testY$pvalue
  }
  expect_equal(X, Y, tolerance=0.01)
})



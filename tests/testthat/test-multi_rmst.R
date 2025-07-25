test_that("The given data frame is correctly shaped", {
  df = data_under_PH
  df["random"] = runif(nrow(df))
  expect_error(multi_rmst(df[c("time","status","random")]),
               "The dataframe must contain the columns 'time', 'status' and 'arm'.")
})

test_that("The arguments given are correct", {
  # Stop if status is not 0 or 1
  df = data_under_PH
  df$status[df$status == 0] = -1
  expect_error(multi_rmst(df), "'status' must be either 0 or 1.")
  
  # Stop if arm is not 0,...,n-1
  df = data_under_PH
  df$arm[df$arm == 1] = 3
  expect_error(multi_rmst(df[1:200,]), "Need at least two groups.")
})


test_that("The output is correct", {
  testthat::skip_if_not_installed("survRM2")
  X = multi_rmst(data_not_PH, tau=30)$results
  X = X[,1:3]
  X = unname(X)
  
  ind=matrix(c(0,1,0,2,1,2),nrow=3,byrow=TRUE)
  Y = matrix(NA, nrow=3, ncol=3)
  for (k in 1:3){
    i=ind[k,1]
    j=ind[k,2]
    index_ij = (data_not_PH$arm == i) | (data_not_PH$arm == j)
    df_ij = data_not_PH[index_ij,]
    df_ij$arm = (df_ij$arm - i) / (j-i)
    
    test = survRM2::rmst2(df_ij$time, df_ij$status, df_ij$arm, tau=30)
    Y[k,1] = test$unadjusted.result[1,1]
    Y[k,2] = sqrt(test$RMST.arm0$rmst[2]^2 + test$RMST.arm1$rmst[2]^2)
    Y[k,3] = test$unadjusted.result[1,4]
  }
  expect_equal(X, Y, tolerance=0.05)
})


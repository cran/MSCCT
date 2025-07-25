test_that("The given data frame is correctly shaped", {
  df = data_under_PH
  df["random"] = runif(nrow(df))
  expect_error(multi_ts(df[c("time","status","random")]),
               "The dataframe must contain the columns 'time', 'status' and 'arm'.")
})


test_that("The arguments given are correct", {
  # Stop if status is not 0 or 1
  df = data_under_PH
  df$status[df$status == 0] = -1
  expect_error(multi_ts(df), "'status' must be either 0 or 1.")
  
  # Stop if arm is not 0,...,n-1 and n >= 2
  df = data_under_PH
  df$arm[df$arm == 1] = 3
  expect_error(multi_ts(data_under_PH[1:200,]), "Need at least two groups.")
})



test_that("The results are correct", {
  testthat::skip_if_not_installed("TSHRC")
  X = multi_ts(data_not_PH, nboot=100, eps=0.25)$results
  X = X[,1:2]
  X = unname(X)
  
  ind=matrix(c(0,1,0,2,1,2),nrow=3,byrow=TRUE)
  Y = matrix(NA, nrow=3, ncol=2)
  for (k in 1:3){
    i=ind[k,1]
    j=ind[k,2]
    index_ij = (data_not_PH$arm == i) | (data_not_PH$arm == j)
    df_ij = data_not_PH[index_ij,]
    df_ij$arm = (df_ij$arm - i) / (j-i)
    
    test = TSHRC::twostage(df_ij$time, df_ij$status, df_ij$arm, nboot=100, eps=0.25)
    Y[k,] = unname(test)[1:2]
  }
  expect_equal(X, Y, tolerance=0.05)
})


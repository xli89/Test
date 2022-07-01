# 
regression.control <- function(fit){
  regression.control <- tryCatch(
    {
      var.cov <- sqrt(diag(solve(fit$hessian)))
      regression.control <- 2
    },error=function(cond) {
      #          print("Regression caused an error")
      #          print("Here's the original error message:")
      #          print(cond)
      # Choose a return value in case of error
      return(regression.control=-1)
    },
    warning=function(cond) {
      #          print("Regression caused a warning:")
      #          print("Here's the original warning message:")
      #          print(cond)
      regression.control=1
      return(regression.control=1)
    }
  )
  return(regression.control)
}
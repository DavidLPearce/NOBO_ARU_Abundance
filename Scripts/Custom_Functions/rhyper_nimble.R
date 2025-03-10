# Define Hypergeometric Sampling Function
rhyper_nimble <- nimbleFunction(
  run = function(n = integer(0), K = double(0), Q = double(0), n_draws = double(0)) {
    returnType(double(0))  # Must match `x` in `dhyper_nimble
    
    if (n != 1) print("rhyper_nimble only allows n = 1; using n = 1.")
    
    x <- 0
    K_remain <- K
    Q_remain <- Q
    
    for (j in 1:n_draws) {
      if (runif(1) < (K_remain / (K_remain + Q_remain))) {
        x <- x + 1
        K_remain <- K_remain - 1
      } else {
        Q_remain <- Q_remain - 1
      }
    }
    return(x)  # Returns a single value instead of a vector
  }
)

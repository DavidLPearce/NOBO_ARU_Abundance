# Define Hypergeometric Probability Function
dhyper_nimble <- nimbleFunction(
  run = function(x = double(0), K = double(0), Q = double(0), n_draws = double(0), log = integer(0, default = 0)) {
    returnType(double(0))
    
    # Ensure valid input values (replace `||` with separate `if` statements)
    if (x < 0) return(-Inf)
    if (x > K) return(-Inf)
    if ((n_draws - x) > Q) return(-Inf)
    
    # Compute hypergeometric probability
    prob <- exp(lgamma(K + 1) - lgamma(x + 1) - lgamma(K - x + 1) +
                  lgamma(Q + 1) - lgamma(n_draws - x + 1) - lgamma(Q - (n_draws - x) + 1) -
                  (lgamma(K + Q + 1) - lgamma(n_draws + 1) - lgamma(K + Q - n_draws + 1)))
    
    # Return log probability if requested
    if (log) return(log(prob))
    return(exp(log_prob))
  }
)

dhyper_nimble <- nimbleFunction(
  run = function(x = double(0), K = double(0), Q = double(0), n_draws = double(0), log = integer(0, default = 0)) {
    returnType(double(0))
    
    # Check for missing values (NA handling)
    if (is.nan(x) | is.nan(K) | is.nan(Q) | is.nan(n_draws)) {
      if (log) return(-Inf)
      return(0)
    }
    
    if (x < 0 | is.na(x)) {
      if (log) return(-Inf)
      return(0)
    }
    if (K < 0 | is.na(K)) {
      if (log) return(-Inf)
      return(0)
    }
    if (Q < 0 | is.na(Q)) {
      if (log) return(-Inf)
      return(0)
    }
    if (n_draws < 0 | is.na(n_draws)) {
      if (log) return(-Inf)
      return(0)
    }
    
    # Ensure valid hypergeometric conditions
    if (x > K) {
      if (log) return(-Inf)
      return(0)
    }
    if ((n_draws - x) > Q) {
      if (log) return(-Inf)
      return(0)
    }
    
    # Compute hypergeometric probability
    log_prob <- lgamma(K + 1) - lgamma(x + 1) - lgamma(K - x + 1) +
      lgamma(Q + 1) - lgamma(n_draws - x + 1) - lgamma(Q - (n_draws - x) + 1) -
      (lgamma(K + Q + 1) - lgamma(n_draws + 1) - lgamma(K + Q - n_draws + 1))
    
    # Return log probability if requested
    if (log) return(log_prob)
    return(exp(log_prob))
  }
)

# Binomial Coefficient Function
choose_nimble <- nimbleFunction(
  run = function(n = double(0), k = double(0)) {  # Removed extra `)`
    returnType(double(0))
    if (k > n) return(0)
    if (k < 0) return(0)
    return(exp(logFactorial(n) - logFactorial(k) - logFactorial(n - k)))
  }
)
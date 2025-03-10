# Log-Factorial Function for numerical stability
logFactorial <- nimbleFunction(
  run = function(n = double(0)) {  # Removed extra `)` here
    returnType(double(0))
    if (n == 0) return(0)
    result <- 0
    for (i in 1:n) {
      result <- result + log(i)
    }
    return(result)
  }
)
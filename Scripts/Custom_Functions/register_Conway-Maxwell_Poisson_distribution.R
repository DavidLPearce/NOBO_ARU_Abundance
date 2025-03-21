# Register Conway-Maxwell Poisson function
registerDistributions(
  list(
    dCMPois = list(   
      BUGSdist = "dCMPois(lambda, nu)",  
      Rdist = "dCMPois(lambda, nu)",  
      types = c('value = double(0)', 'lambda = double(0)', 'nu = double(0)'),
      pqAvail = FALSE
    )
  )
)

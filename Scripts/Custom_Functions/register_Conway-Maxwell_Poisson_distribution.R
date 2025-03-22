# Register Conway-Maxwell Poisson function
registerDistributions(
  list(
    dCMPois_nimble = list(   
      BUGSdist = "dCMPois_nimble(lambda, nu)",  
      Rdist = "dCMPois_nimble(lambda, nu)",  
      types = c('value = double(0)', 'lambda = double(0)', 'nu = double(0)'),
      pqAvail = FALSE
    )
  )
)

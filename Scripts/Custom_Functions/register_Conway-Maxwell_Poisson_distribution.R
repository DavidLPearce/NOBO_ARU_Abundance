# Register Conway-Maxwell Poisson function
registerDistributions(
  list(
    dCOMPois_nimble = list(   
      BUGSdist = "dCOMPois_nimble(lambda, nu)",  
      Rdist = "dCOMPois_nimble(lambda, nu)",  
      types = c('value = double(0)', 'lambda = double(0)', 'nu = double(0)'),
      pqAvail = FALSE
    )
  )
)

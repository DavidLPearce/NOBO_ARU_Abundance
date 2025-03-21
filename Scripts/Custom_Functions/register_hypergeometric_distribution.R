# Register hypergeometric function
registerDistributions(
  list(
    dhyper_nimble = list(   
      BUGSdist = "dhyper_nimble(m, n, k)",  
      Rdist = "dhyper_nimble(m, n, k)",  
      types = c('value = double(0)', 'm = double(0)', 'n_draws = double(0)', 'k = double(0)'),
      pqAvail = FALSE
    )
  )
)

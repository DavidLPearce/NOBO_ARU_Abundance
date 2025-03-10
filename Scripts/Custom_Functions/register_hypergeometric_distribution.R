# Register hypergeometric function
registerDistributions(
  list(
    dhyper_nimble = list(  # Must match the name of your function!
      BUGSdist = "dhyper_nimble(K, Q, n_draws)",  # Adjust parameter names
      Rdist = "dhyper_nimble(K, Q, n_draws)",  
      types = c('value = double(0)', 'K = double(0)', 'Q = double(0)', 'n_draws = double(0)'),
      pqAvail = FALSE
    )
  )
)

# Define Hypergeometric Sampling Function from R
rhyper_nimble <- nimbleRcall(
  # n_blk is n since NIMBLE expects first argument in simulation function to be n
  function(n = double(0), m = double(0), n_draws  = double(0), k = double(0)){},
  Rfun = 'rhyper',
  returnType = double()
)

# Define in global environment
assign('rhyper_nimble', rhyper_nimble, envir = .GlobalEnv)

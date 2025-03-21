# Define Hypergeometric Sampling Function from R
dhyper_nimble <- nimbleRcall(
  function(x = double(0), m = double(0), n = double(0), k = double(0), log = integer(0, default = 0)){},
  Rfun = 'dhyper',
  returnType = double()
)

# Define in global environment
assign('dhyper_nimble', dhyper_nimble, envir = .GlobalEnv)
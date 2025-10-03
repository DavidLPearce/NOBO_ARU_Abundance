# Define Hypergeometric Sampling Function from R
dhyper_nimble <- nimbleRcall(
  # Have to use R notation, but for model
  # x = manually validated calls found to be true (k)
  # m = number of true detections (K)
  # n = (is n in R) number of false positive detections (Q)
  # k = number of calls validated (n)
  function(x = double(0), m = double(0), n = double(0), k = double(0), log = double(0, default = 0)){},
  Rfun = 'dhyper',
  returnType = double()
)

# Define in global environment
assign('dhyper_nimble', dhyper_nimble, envir = .GlobalEnv)
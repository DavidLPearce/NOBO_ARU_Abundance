# Conway-Maxwell Poisson random values function
rCMPois_nimble <- nimbleRcall(
  function(n=double(0), lambda=double(0), nu=double(0)){},
  Rfun = 'rcmp',
  returnType=double()) 

# Define in global environment
assign('rCMPois_nimble', rCMPois_nimble, envir=.GlobalEnv)

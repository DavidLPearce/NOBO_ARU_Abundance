# Conway-Maxwell Poisson random values function
rCMPois <- nimbleRcall(
  function(n=double(0), lambda=double(0), nu=double(0)){},
  Rfun = 'rcmp',
  returnType=double()) 

# Define in global environment
assign('rCMPois', rCMPois, envir=.GlobalEnv)

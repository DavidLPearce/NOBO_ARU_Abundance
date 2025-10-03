# Function is from the COMPoissonReg package
require(COMPoissonReg)

# Conway-Maxwell Poisson random values function
rCOMPois_nimble <- nimbleRcall(
  function(n=double(0), lambda=double(0), nu=double(0)){},
  Rfun = 'rcmp',
  returnType=double()) 

# Define in global environment
assign('rCOMPois_nimble', rCOMPois_nimble, envir=.GlobalEnv)

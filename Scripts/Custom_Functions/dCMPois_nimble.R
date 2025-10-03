# Function is from the COMPoissonReg package
require(COMPoissonReg)

# Conway-Maxwell Poisson density function
dCOMPois_nimble <- nimbleRcall(
  function(x = double(0), lambda = double(0), nu = double(0), log = double(0, default = 0)){},
  Rfun = 'dcmp',
  returnType=double()) 

# Define in global environment
assign('dCOMPois_nimble', dCOMPois_nimble, envir=.GlobalEnv)


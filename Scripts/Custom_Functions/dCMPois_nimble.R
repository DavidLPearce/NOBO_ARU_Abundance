# Conway-Maxwell Poisson density function
dCMPois_nimble <- nimbleRcall(
  function(x = double(0), lambda = double(0), nu = double(0), log = double(0, default = 0)){},
  Rfun = 'dcmp',
  returnType=double()) 

# Define in global environment
assign('dCMPois_nimble', dCMPois_nimble, envir=.GlobalEnv)


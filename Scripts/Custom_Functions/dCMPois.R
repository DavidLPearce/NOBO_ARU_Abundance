# Conway-Maxwell Poisson density function
dCMPois <- nimbleRcall(
  function(x = double(0), lambda = double(0), nu = double(0), log = double(0, default = 0)){},
  Rfun = 'dcmp',
  returnType=double()) 

# Define in global environment
assign('dCMPois', dCMPois, envir=.GlobalEnv)


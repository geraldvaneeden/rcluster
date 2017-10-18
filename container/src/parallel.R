primary <- '172.17.0.4'
machineAddresses <- list(
  list(host=primary,user='johnmount',
       ncore=1),
  list(host='172.17.0.5',user='johnmount',
       ncore=1),
  list(host='172.17.0.6',user='johnmount',
       ncore=1),
  list(host='172.17.0.7',user='johnmount',
       ncore=1)
)

spec <- lapply(machineAddresses,
               function(machine) {
                 rep(list(list(host=machine$host,
                               user=machine$user)),
                     machine$ncore)
               })
spec <- unlist(spec,recursive=FALSE)

parallelCluster <- parallel::makeCluster(type='PSOCK',
                                         master=primary,
                                         spec=spec)
print(parallelCluster)

# calculate HFD based on daily status (wide format)
# 0 - alive and outside hosp
# 1 - alive and in hosp
# 2 - died
get_HFD = function(sim.wide){
  m <- ncol(sim.wide)-2
  HFD <- apply(sim.wide[,paste0("d.",1:m)],1,function(x) sum(x==0))
  return(HFD)
}
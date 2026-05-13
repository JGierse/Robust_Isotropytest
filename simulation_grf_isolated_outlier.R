################################################################################
## Simulation for the GRF with isolated outliers (see section 4.2)            ##
################################################################################

##############################
## Load functions and packages
##############################

library(readxl)
library(RobVario)
library(future) # for parallelisation  
library(future.apply) # for parallelisation 
source("Parameters.R")

########################################################
## Parameter combinations for the scenarios of the paper
########################################################

combs <- read_excel("Parametercombinations.xlsx")
combs.iso <- combs[which(combs$Outlier.type == "isolated"),]

#############
## Simulation
#############

# needed to run the different scenarios on the HPC Cluster
ind <- as.integer(Sys.getenv("PBS_ARRAYID")) 
print(ind)

# load the data for the corresponding scenario
load(paste0("Data/Iso/data.iso.vario", combs.iso$variogram[ind], ".range", combs.iso$range[ind], ".aniso.r", combs.iso$rotation[ind],
            ".s", combs.iso$scale[ind], ".nugget", combs.iso$nugget[ind], ".amount", combs.iso$Amount[ind], ".dist", combs.iso$Dist[ind], ".RData"))


plan(multicore)  # need for the HPC cluster

# save in the needed format
datas <- lapply(1:1000, function(l) cbind(data.iso$data[,l], data.iso$grid))

# define the window sample size
if(combs.iso$Method[ind] == "subsampling"){
  w <- window.subsampling[combs.iso$gridsize[ind],]
}
if(combs.iso$Method[ind] == "blockpermutation"){
  if(combs.iso$gridsize[ind] == 1){
    w <- window.blockpermutation[[combs.iso$gridsize[ind]]] 
  }
  if(combs.iso$gridsize[ind] == 1){
    w <- window.blockpermutation[[combs.iso$gridsize[ind]]][combs.iso$Window[ind],] 
  }
}

res_iso <- future_lapply(1:1000, function(x) isotropy_test(datas[[x]], lagmat = lags[[combs.iso$LAGS[ind]]], A = As[[combs.iso$LAGS[ind]]], 
                                                           method = combs.iso$Method[ind], window.dims = w, var.robust = TRUE, edge.sub = TRUE),
                         future.seed = TRUE)

save(res_non, file = paste0("Results/Iso/", combs.iso$Method[ind],".vario", combs.iso$variogram[ind], ".grid", combs.iso$gridsize[ind]  ,".range", combs.iso$range[ind],
                            ".nugget", combs.iso$nugget[ind],".l", combs.iso$LAGS[ind], ".rotation", combs.iso$rotation[ind], ".scale", combs.iso$scale[ind], 
                            ".amount", combs.iso$Amount[ind], ".dist", combs.iso$Dist[ind], ".RData"))

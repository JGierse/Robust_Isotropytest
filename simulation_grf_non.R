################################################################################
## Simulation for the GRF without outliers (see section 4.1)                  ##
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
combs.non <- combs[which(combs$Outlier.type == "NULL"),]

#############
## Simulation
#############

# needed to run the different scenarios on the HPC Cluster
ind <- as.integer(Sys.getenv("PBS_ARRAYID")) 
print(ind)

# load the data for the corresponding scenario
load(paste0("Data/Non/data.non.vario", combs.non$variogram[ind], ".range", combs.non$range[ind], ".aniso.r", combs.non$rotation[ind],
            ".s", combs.non$scale[ind], ".nugget", combs.non$nugget[ind], ".RData"))

plan(multicore)  # need for the HPC cluster

# save in the needed format
datas <- lapply(1:1000, function(l) cbind(data.non$data[,l], data.non$grid))

# define the window sample size
if(combs.non$Method[ind] == "subsampling"){
  w <- window.subsampling[combs.non$gridsize[ind],]
}
if(combs.non$Method[ind] == "blockpermutation"){
  if(combs.non$gridsize[ind] == 1){
    w <- window.blockpermutation[[combs.non$gridsize[ind]]] 
  }
  if(combs.non$gridsize[ind] == 1){
    w <- window.blockpermutation[[combs.non$gridsize[ind]]][combs.non$Window[ind],] 
  }
}

res_non <- future_lapply(1:1000, function(x) isotropy_test(datas[[x]], lagmat = lags[[combs.non$LAGS[ind]]], A = As[[combs.non$LAGS[ind]]], 
                                                           method = combs.non$Method[ind], window.dims = w, var.robust = TRUE, edge.sub = TRUE),
                         future.seed = TRUE)

save(res_non, file = paste0("Results/Non/", combs.non$Method[ind],".vario", combs.non$variogram[ind], ".grid", combs.non$gridsize[ind]  ,".range", combs.non$range[ind],
                            ".nugget", combs.non$nugget[ind],".l", combs.non$LAGS[ind], 
                            ".rotation", combs.non$rotation[ind], ".scale", combs.non$scale[ind], 
                            ".window", combs.non$Window[ind], ".RData"))

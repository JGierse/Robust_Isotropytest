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
combs.block <- combs[which(combs$Outlier.type == "block"),]

#############
## Simulation
#############

# needed to run the different scenarios on the HPC Cluster
ind <- as.integer(Sys.getenv("PBS_ARRAYID")) 
print(ind)

# load the data for the corresponding scenario
load(paste0("Data/Iso/data.block.", combs.block$Block.type[ind], ".vario", combs.block$variogram[ind], ".range", combs.block$range[ind], ".aniso.r", combs.block$rotation[ind],
            ".s", combs.block$scale[ind], ".nugget", combs.block$nugget[ind], ".amount", combs.block$Amount[ind], ".dist", combs.block$Dist[ind], ".RData"))


plan(multicore)  # need for the HPC cluster

# save in the needed format
datas <- lapply(1:1000, function(l) cbind(data.block$data[,l], data.block$grid))

# define the window sample size
if(combs.block$Method[ind] == "subsampling"){
  w <- window.subsampling[combs.block$gridsize[ind],]
}
if(combs.block$Method[ind] == "blockpermutation"){
  if(combs.block$gridsize[ind] == 1){
    w <- window.blockpermutation[[combs.block$gridsize[ind]]] 
  }
  if(combs.block$gridsize[ind] == 1){
    w <- window.blockpermutation[[combs.block$gridsize[ind]]][combs.block$Window[ind],] 
  }
}

res_block <- future_lapply(1:1000, function(x) isotropy_test(datas[[x]], lagmat = lags[[combs.block$LAGS[ind]]], A = As[[combs.block$LAGS[ind]]], 
                                                           method = combs.block$Method[ind], window.dims = w, var.robust = TRUE, edge.sub = TRUE),
                         future.seed = TRUE)

save(res_block, file = paste0("Results/Block/", combs.block$Method[ind], " block.", combs.block$Block.type[ind], ".vario", combs.block$variogram[ind], ".grid", combs.block$gridsize[ind]  ,".range", combs.block$range[ind],
                            ".nugget", combs.block$nugget[ind],".l", combs.block$LAGS[ind], ".rotation", combs.block$rotation[ind], ".scale", combs.block$scale[ind], 
                            ".window", combs.block$Window[ind], ".amount", combs.block$Amount[ind], ".dist", combs.block$Dist[ind], ".RData"))

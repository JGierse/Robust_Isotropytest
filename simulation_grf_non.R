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

##############
## Evaluation
##############

## calculate the empirical rejection rates
res <- array(dim = c(3,2,2,2,3,3,3,3), dimnames = list("estimator" = c("Matheron", "Genton", "MCD.diff"),
                                                   "method" = c("subsampling", "blockpermutation"),
                                                   "grid" = c("24x24", "40x40"),
                                                   "rotation" = c(0, 0.79),
                                                   "scale" = c(0, 1.41, 2),
                                                   "range"  = c(2,5,8),
                                                   "Lags" = c(1,2,3),
                                                   "windos" = c(1,2,3)))

for(ind in 1:nrow(combs.non)){
  load(paste0("Results/Non/", combs.non$Method[ind],".vario", combs.non$variogram[ind], ".grid", combs.non$gridsize[ind]  ,".range", combs.non$range[ind],
              ".nugget", combs.non$nugget[ind],".l", combs.non$LAGS[ind], 
              ".rotation", combs.non$rotation[ind], ".scale", combs.non$scale[ind], 
              ".window", combs.non$Window[ind], ".RData"))
  
  
  mat <- lapply(res_non, function(x){
    x$Matheron$p.value
  })
  mat <- unlist(mat)
  
  gen <- lapply(res_non, function(x){
    x$Genton$p.value
  })
  gen <- unlist(gen)
  
  mcd <- lapply(res_non, function(x){
    x$MCD$p.value
  })
  mcd <- unlist(mcd)
  
  ifelse(combs.non$Method[ind] == "subsampling", m <- 1, m <- 2)
  if(combs.non$rotation[ind] == 0) r <- 1
  if(combs.non$rotation[ind] == pi/4) r <- 2
  if(combs.non$scale[ind] == 1) s <- 1
  if(combs.non$scale[ind] == sqrt(2)) s <- 2
  if(combs.non$scale[ind] == 2) s <- 3
  if(combs.non$range[ind] == 2) ra <- 1
  if(combs.non$range[ind] == 5) ra <- 2
  if(combs.non$range[ind] == 8) ra <- 3
  
  
  res[1,m,combs.non$gridsize[ind],r,s,ra,combs.non$LAGS[ind],combs.non$Window[ind]] <- mean(mat <= 0.05)
  res[2,m,combs.non$gridsize[ind],r,s,ra,combs.non$LAGS[ind],combs.non$Window[ind]] <- mean(gen <= 0.05)
  res[3,m,combs.non$gridsize[ind],r,s,ra,combs.non$LAGS[ind],combs.non$Window[ind]] <- mean(mcd <= 0.05)
}

save(res, file = "Results/res_non.RData")

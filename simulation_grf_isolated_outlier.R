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


##############
## Evaluation
##############

## calculate the empirical rejection rates
res <- array(dim = c(3,2,2,3,3,3,2,2), dimnames = list("estimator" = c("Matheron", "Genton", "MCD.diff"),
                                                       "method" = c("subsampling", "blockpermutation"),
                                                       "rotation" = c(0, 0.79),
                                                       "scale" = c(0, 1.41, 2),
                                                       "range"  = c(2,5,8),
                                                       "Lags" = c(1,2,3),
                                                       "amount" = c(0.1, 0.2),
                                                       "dist" = c(1,2)))

for(ind in 1:nrow(combs.iso)){
  load(paste0("Results/Iso/", combs.iso$Method[ind],".vario", combs.iso$variogram[ind], ".grid", combs.iso$gridsize[ind]  ,".range", combs.iso$range[ind],
              ".nugget", combs.iso$nugget[ind],".l", combs.iso$LAGS[ind], ".rotation", combs.iso$rotation[ind], ".scale", combs.iso$scale[ind], 
              ".amount", combs.iso$Amount[ind], ".dist", combs.iso$Dist[ind], ".RData"))
  
  
  mat <- lapply(res_iso, function(x){
    x$Matheron$p.value
  })
  mat <- unlist(mat)
  
  gen <- lapply(res_iso, function(x){
    x$Genton$p.value
  })
  gen <- unlist(gen)
  
  mcd <- lapply(res_iso, function(x){
    x$MCD$p.value
  })
  mcd <- unlist(mcd)
  
  ifelse(combs.iso$Method[ind] == "subsampling", m <- 1, m <- 2)
  if(combs.iso$rotation[ind] == 0) r <- 1
  if(combs.iso$rotation[ind] == pi/4) r <- 2
  if(combs.iso$scale[ind] == 1) s <- 1
  if(combs.iso$scale[ind] == sqrt(2)) s <- 2
  if(combs.iso$scale[ind] == 2) s <- 3
  if(combs.iso$range[ind] == 2) ra <- 1
  if(combs.iso$range[ind] == 5) ra <- 2
  if(combs.iso$range[ind] == 8) ra <- 3
  ifelse(combs.iso$Amount[ind] == 0.1, a <- 1, a <- 2)
  
  
  res[1,m,r,s,ra,combs.iso$LAGS[ind],a,combs.iso$Dist[ind]] <- mean(mat <= 0.05)
  res[2,m,r,s,ra,combs.iso$LAGS[ind],a,combs.iso$Dist[ind]] <- mean(gen <= 0.05)
  res[3,m,r,s,ra,combs.iso$LAGS[ind],a,combs.iso$Dist[ind]] <- mean(mcd <= 0.05)
}

save(res, file = "Results/res_iso.RData")


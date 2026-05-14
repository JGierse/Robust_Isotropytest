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
                              ".amount", combs.block$Amount[ind], ".dist", combs.block$Dist[ind], ".RData"))

##############
## Evaluation
##############

## calculate the empirical rejection rates
res <- array(dim = c(3,2,2,3,3,3,2,2,3), dimnames = list("estimator" = c("Matheron", "Genton", "MCD.diff"),
                                                         "method" = c("subsampling", "blockpermutation"),
                                                         "rotation" = c(0, 0.79),
                                                         "scale" = c(0, 1.41, 2),
                                                         "range"  = c(2,5,8),
                                                         "Lags" = c(1,2,3),
                                                         "amount" = c(0.1, 0.2),
                                                         "dist" = c(1,2),
                                                         "blocktype" = c("random", "rectangle", "square")))

for(ind in 1:nrow(combs.block)){
  load(paste0("C:/Users/paulg/Sciebo2/Paper_2/R-Code/Robust-Isotropytest/Paper/Results/Block/", combs.block$Method[ind], ".block.", combs.block$Block.type[ind], ".vario", combs.block$variogram[ind], ".grid", combs.block$gridsize[ind]  ,".range", combs.block$range[ind],
              ".nugget", combs.block$nugget[ind],".l", combs.block$LAGS[ind], ".rotation", combs.block$rotation[ind], ".scale", combs.block$scale[ind], 
              ".amount", combs.block$Amount[ind], ".dist", combs.block$Dist[ind], ".RData"))
  
  
  mat <- lapply(res_block, function(x){
    x$Matheron$p.value
  })
  mat <- unlist(mat)
  
  gen <- lapply(res_block, function(x){
    x$Genton$p.value
  })
  gen <- unlist(gen)
  
  mcd <- lapply(res_block, function(x){
    x$MCD$p.value
  })
  mcd <- unlist(mcd)
  
  ifelse(combs.block$Method[ind] == "subsampling", m <- 1, m <- 2)
  if(combs.block$rotation[ind] == 0) r <- 1
  if(combs.block$rotation[ind] == pi/4) r <- 2
  if(combs.block$scale[ind] == 1) s <- 1
  if(combs.block$scale[ind] == sqrt(2)) s <- 2
  if(combs.block$scale[ind] == 2) s <- 3
  if(combs.block$range[ind] == 2) ra <- 1
  if(combs.block$range[ind] == 5) ra <- 2
  if(combs.block$range[ind] == 8) ra <- 3
  ifelse(combs.block$Amount[ind] == 0.1, a <- 1, a <- 2)
  if(combs.block$Block.type[ind] == "random") b <- 1
  if(combs.block$Block.type[ind] == "rectangle") b <- 2
  if(combs.block$Block.type[ind] == "square") b <- 3
  
  res[1,m,r,s,ra,combs.block$LAGS[ind],a,combs.block$Dist[ind],b] <- mean(mat <= 0.05)
  res[2,m,r,s,ra,combs.block$LAGS[ind],a,combs.block$Dist[ind],b] <- mean(gen <= 0.05)
  res[3,m,r,s,ra,combs.block$LAGS[ind],a,combs.block$Dist[ind],b] <- mean(mcd <= 0.05)
}

save(res, file = "Results/res_block.RData")


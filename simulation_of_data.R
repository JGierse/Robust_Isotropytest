################################################################################
## Simulation of the data, i.e. simulation of the Gaussian random fields      ##
################################################################################


##############################
## Load functions and packages
##############################

library(readxl)
library(RobVario)
source("Parameters.R")

########################################################
## Parameter combinations for the scenarios of the paper
########################################################

combs <- read_excel("Parametercombinations.xlsx")


###########################################
## Gaussian random fields: without outliers
###########################################

combs.non <- combs[which(combs$Outlier.type == "NULL"),]
combs.non <- combs.non[, -c(7,8,9,10,11,12,13,14)] # not necessary for simulation of data
combs.non <- unique(combs.non)
seeds <- rep(0, nrow(combs.non)) # to save the seeds

for(ind in 1:nrow(combs.non)){
  s <- sample(1:10000, 1) # draw a random seed
  seeds[ind] <- s
  
  if(combs.non$rotation[ind] != 0| combs.non$scale[ind] != 1){
    aniso <- c(combs.non$rotation[ind], combs.non$scale[ind])
  } else{
    aniso <- NULL
  }
  
  data.non <- simulate_grf(gridsize = size[combs.non$gridsize[ind],],
                           variogram = combs.non$variogram[ind],
                           param.variogram = c(1, combs.non$range[ind]),
                           nugget = as.numeric(combs.non$nugget[ind]),
                           aniso.param = aniso,
                           seed = s)
  save(data.non, file = paste0("Data/Non/data.non.vario", combs.non$variogram[ind], ".range", combs.non$range[ind], ".aniso.r", combs.non$rotation[ind],
                               ".s", combs.non$scale[ind], ".nugget", combs.non$nugget[ind], ".RData"))
}
save(seeds, file = "Data/Non/seeds.RData")


################################################################
## Gaussian random fields: for simulation with isolated outliers
################################################################

combs.iso <- combs[which(combs$Outlier.type == "isolated"),]
combs.iso <- combs.iso[, -c(7,8,9,12,13,14)] # not necessary for simulation of data
combs.iso <- unique(combs.iso)
seeds <- rep(0, nrow(combs.iso))

for(ind in 1:nrow(combs.iso)){
  s <- sample(1:10000, 1) # draw a random seed
  seeds[ind] <- s
  
  if(combs.iso$rotation[ind] != 0| combs.iso$scale[ind] != 1){
    aniso <- c(combs.iso$rotation[ind], combs.iso$scale[ind])
  } else{
    aniso <- NULL
  }
  
  data.iso <- sim.data(gridsize = size[combs.iso$gridsize[ind],],
                       variogram = combs.iso$variogram[ind],
                       param.variogram = c(1, combs.iso$range[ind]),
                       nugget = as.numeric(combs.iso$nugget[ind]),
                       aniso.param = aniso,
                       out.type = "isolated",
                       amount = as.numeric(combs.iso$Amount[ind]),
                       param.outlier = dists[combs.iso$Dist[ind],])
  save(data.iso, file = paste0("Data/Iso/data.iso.vario", combs.iso$variogram[ind], ".range", combs.iso$range[ind], ".aniso.r", combs.iso$rotation[ind],
                               ".s", combs.iso$scale[ind], ".nugget", combs.iso$nugget[ind], 
                               ".amount", combs.iso$Amount[ind], ".dist", combs.iso$Dist[ind], ".RData"))
}

save(seeds, file = "Data/Iso/seeds.RData")



##############################################################
## Gaussian random fiels: for simulation with an outlier block
##############################################################

combs.block <- combs[which(combs$Outlier.type == "block"),]
combs.block <- combs.block[, -c(7,8,9,12,14)] # not necessary for simulation of data
combs.block <- unique(combs.block)
seeds <- rep(0, nrow(combs.block))


for(ind in 1:nrow(combs.block)){
  s <- sample(1:10000, 1) # draw a random seed
  seeds[ind] <- s
  
  if(combs.block$rotation[ind] != 0| combs.block$scale[ind] != 1){
    aniso <- c(combs.block$rotation[ind], combs.block$scale[ind])
  } else{
    aniso <- NULL
  }
  
  data.block <- sim.data(gridsize = size[combs.block$gridsize[ind],],
                         variogram = combs.block$variogram[ind],
                         param.variogram = c(1, combs.block$range[ind]),
                         nugget = as.numeric(combs.block$nugget[ind]),
                         aniso.param = aniso,
                         out.type = "block",
                         block.type = combs.block$Block.type[ind],
                         amount = as.numeric(combs.block$Amount[ind]),
                         param.outlier = dists[combs.block$Dist[ind],])
  
  save(data.block, file = paste0("Data/Block/data.block.", combs.block$Block.type[ind], ".vario", combs.block$variogram[ind], ".range", combs.block$range[ind], ".aniso.r", combs.block$rotation[ind],
                                 ".s", combs.block$scale[ind], ".nugget", combs.block$nugget[ind], 
                                 ".amount", combs.block$Amount[ind], ".dist", combs.block$Dist[ind], ".RData"))
}

save(seeds, file = "Data/Iso/seeds.RData")



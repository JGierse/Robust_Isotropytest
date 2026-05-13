################################################################################
## Runtimes of the different tests                                            ##
################################################################################

################
## Load packages
################

library(microbenchmark)
library(RobVario)


#############
## Parameters
#############

# gridsize 
grids <- rbind(c("nx" = 24, "ny" = 24),
               c("nx" = 40, "ny" = 40))


# lagmat 
lags <- list(rbind(c(1,0), c(0,1)),
             rbind(c(1,0), c(0,1), c(1,1), c(1, -1)),
             rbind(c(1,0), c(0,1), c(1,1), c(1,-1), c(2,1), c(-1, 2), c(1, 2), c(-2, 1)))

# A
As <- list(rbind(c(1, -1)),
           rbind(c(1, -1, 0, 0), c(0, 0, 1, -1)),
           rbind(c(1, -1, 0, 0, 0, 0, 0, 0), c(0, 0, 1, -1, 0, 0, 0, 0), c(0, 0, 0, 0, 1, -1, 0, 0), c(0, 0, 0, 0, 0, 0, 1, -1)))


# estimator 
estimator <- c("Matheron", "Genton", "MCD")


# window.dims
window.subsampling <- rbind(c(5,5),
                            c(6,6),
                            c(8,8))

window.blockpermutation <- list(matrix(c(6,6), nrow = 1),
                                rbind(c(5,5), c(8,8), c(10,10)),
                                rbind(c(6,6), c(10,10), c(12,12)))


# method
Method <- c("Subsampling", "Blockpermutation")


# grid
para.grid <- expand.grid("size" = 1:3, "lag" = 1:3, "estimator" = 1:3, "method" = 1:2)


#############
## Simulation
#############

# needed to run the different scenarios on the HPC Cluster
ind <- as.integer(Sys.getenv("PBS_ARRAYID")) 
print(ind)

# load Data 
if(para.grid[ind, "size"] == 1){
  load(paste0("Data/Non/data.non.variospherical.range5.aniso.r0.s0.nugget0.RData"))
  
} else{
  load(paste0("Data/Non/data.non.grid2variospherical.range5.aniso.r0.s0.nugget0.RData"))
}

datas <- lapply(1:1000, function(l) cbind(data.non$data[,l], data.non$grid))


if(para.grid[ind, "method"] == 1){ # Subsampling
  RunTime <- microbenchmark(isotropy_subsampling(data = datas[[1]], lagmat = lags[[para.grid[ind, "lag"]]], A = As[[para.grid[ind, "lag"]]],
                                                 estimator = estimator[para.grid[ind, "estimator"]],
                                                 window.dims = window.subsampling[para.grid[ind, "size"] ,], edge = TRUE))
  save(RunTime, file = paste0("/work/smjnkoen/Robust_Isotropytest/Laufzeit/grid", para.grid[ind, "size"], ".lag", para.grid[ind, "lag"], ".estimator", para.grid[ind, "estimator"],
                              ".methodSubsampling.RData"))
  
} else{ # blockpermutation
  RunTime1 <- microbenchmark(isotropy_blockpermutation(data = datas[[1]], lagmat = lags[[para.grid[ind, "lag"]]], A = As[[para.grid[ind, "lag"]]],
                                                       estimator = estimator[para.grid[ind, "estimator"]], 
                                                       window.dims = window.blockpermutation[[para.grid[ind, "size"]]][1,], corr.block = FALSE))
  
  if(para.grid[ind, "size"] != 1){
    RunTime2 <- microbenchmark(isotropy_blockpermutation(data = datas[[1]], lagmat = lags[[para.grid[ind, "lag"]]], A = As[[para.grid[ind, "lag"]]],
                                                         estimator = estimator[para.grid[ind, "estimator"]], 
                                                         window.dims = window.blockpermutation[[para.grid[ind, "size"]]][2,], corr.block = FALSE))
    
    RunTime3 <- microbenchmark(isotropy_blockpermutation(data = datas[[1]], lagmat = lags[[para.grid[ind, "lag"]]], A = As[[para.grid[ind, "lag"]]],
                                                         estimator = estimator[para.grid[ind, "estimator"]], 
                                                         window.dims = window.blockpermutation[[para.grid[ind, "size"]]][3,], corr.block = FALSE))
    
    RunTime <- list(RunTime1, RunTime2, RunTime3)
  } else{
    RunTime <- RunTime1
  }
  
  save(RunTime, file = paste0("Runtimes/grid", para.grid[ind, "size"], ".lag", para.grid[ind, "lag"], ".estimator", para.grid[ind, "estimator"],
                              ".methodBlockpermutation.RData"))
}

##############
## Evaluation
##############

## 1. 24 times 24

para.grid1 <- para.grid[which(para.grid$size == 1),]

## Subampling
para.grid1.sub <- para.grid1[which(para.grid1$method == 1),]
nrow(para.grid1.sub)
#[1] 9

# Matheron
# lag 1
# Unit: seconds
# expr
# isotropy_subsampling(data = datas[[1]], lagmat = lags[[para.grid[ind,      "lag"]]], A = As[[para.grid[ind, "lag"]]], estimator = estimator[para.grid[ind,      "estimator"]], window.dims = window.subsampling[para.grid[ind,      "size"], ], edge = TRUE)
# min       lq     mean   median      uq      max neval
# 1.321058 1.361069 1.393351 1.378453 1.41112 1.762024   100

# lag 2
# Unit: seconds
# expr
# isotropy_subsampling(data = datas[[1]], lagmat = lags[[para.grid[ind,      "lag"]]], A = As[[para.grid[ind, "lag"]]], estimator = estimator[para.grid[ind,      "estimator"]], window.dims = window.subsampling[para.grid[ind,      "size"], ], edge = TRUE)
# min       lq     mean   median       uq      max neval
# 1.373082 1.401565 1.428602 1.415313 1.438604 1.874489   100

# lag 3
# Unit: seconds
# expr
# isotropy_subsampling(data = datas[[1]], lagmat = lags[[para.grid[ind,      "lag"]]], A = As[[para.grid[ind, "lag"]]], estimator = estimator[para.grid[ind,      "estimator"]], window.dims = window.subsampling[para.grid[ind,      "size"], ], edge = TRUE)
# min       lq    mean   median      uq      max neval
# 1.604075 1.628492 1.65965 1.651753 1.67727 2.019275   100


# Genton 
# lag 1
# Unit: seconds
# expr
# isotropy_subsampling(data = datas[[1]], lagmat = lags[[para.grid[ind,      "lag"]]], A = As[[para.grid[ind, "lag"]]], estimator = estimator[para.grid[ind,      "estimator"]], window.dims = window.subsampling[para.grid[ind,      "size"], ], edge = TRUE)
# min       lq     mean   median       uq      max neval
# 1.371988 1.407551 1.432468 1.422169 1.439726 1.803238   100

# lag 2
# Unit: seconds
# expr
# isotropy_subsampling(data = datas[[1]], lagmat = lags[[para.grid[ind,      "lag"]]], A = As[[para.grid[ind, "lag"]]], estimator = estimator[para.grid[ind,      "estimator"]], window.dims = window.subsampling[para.grid[ind,      "size"], ], edge = TRUE)
# min     lq     mean   median       uq      max neval
# 1.516844 1.5385 1.559581 1.551214 1.575822 1.929336   100

# lag3
# Unit: seconds
# expr
# isotropy_subsampling(data = datas[[1]], lagmat = lags[[para.grid[ind,      "lag"]]], A = As[[para.grid[ind, "lag"]]], estimator = estimator[para.grid[ind,      "estimator"]], window.dims = window.subsampling[para.grid[ind,      "size"], ], edge = TRUE)
# min       lq    mean   median       uq     max neval
# 1.684165 1.703821 1.73087 1.719759 1.738561 2.09545   100

# MCD
# lag 1
# Unit: seconds
# expr
# isotropy_subsampling(data = datas[[1]], lagmat = lags[[para.grid[ind,      "lag"]]], A = As[[para.grid[ind, "lag"]]], estimator = estimator[para.grid[ind,      "estimator"]], window.dims = window.subsampling[para.grid[ind,      "size"], ], edge = TRUE)
# min       lq     mean   median       uq      max neval
# 2.374824 2.409855 2.420646 2.413966 2.420549 2.827614   100

# lag 2
# Unit: seconds
# expr
# isotropy_subsampling(data = datas[[1]], lagmat = lags[[para.grid[ind,      "lag"]]], A = As[[para.grid[ind, "lag"]]], estimator = estimator[para.grid[ind,      "estimator"]], window.dims = window.subsampling[para.grid[ind,      "size"], ], edge = TRUE)
# min       lq     mean   median       uq      max neval
# 3.038524 3.070549 3.096837 3.083112 3.109052 3.505705   100

# lag 3
# Unit: seconds
# expr
# isotropy_subsampling(data = datas[[1]], lagmat = lags[[para.grid[ind,      "lag"]]], A = As[[para.grid[ind, "lag"]]], estimator = estimator[para.grid[ind,      "estimator"]], window.dims = window.subsampling[para.grid[ind,      "size"], ], edge = TRUE)
# min       lq     mean   median       uq      max neval
# 3.880916 3.926061 3.955887 3.951318 3.974012 4.299655   100


## Blockpermutation
para.grid1.block <- para.grid1[which(para.grid1$method == 2),]
nrow(para.grid1.block)
#[1] 9

# Matheron
# lag 1
# Unit: seconds
# expr
# isotropy_blockpermutation(data = datas[[1]], lagmat = lags[[para.grid[ind,      "lag"]]], A = As[[para.grid[ind, "lag"]]], estimator = estimator[para.grid[ind,      "estimator"]], window.dims = window.blockpermutation[[para.grid[ind,      "size"]]][1, ], corr.block = FALSE)
# min       lq     mean   median       uq      max neval
# 20.70209 20.80373 20.85759 20.84988 20.89941 21.16118   100

# lag 2
# Unit: seconds
# expr
# isotropy_blockpermutation(data = datas[[1]], lagmat = lags[[para.grid[ind,      "lag"]]], A = As[[para.grid[ind, "lag"]]], estimator = estimator[para.grid[ind,      "estimator"]], window.dims = window.blockpermutation[[para.grid[ind,      "size"]]][1, ], corr.block = FALSE)
# min      lq     mean   median       uq      max neval
# 34.778 34.9632 35.03039 35.01847 35.09853 35.61206   100

# lag 3
# Unit: seconds
# expr
# isotropy_blockpermutation(data = datas[[1]], lagmat = lags[[para.grid[ind,      "lag"]]], A = As[[para.grid[ind, "lag"]]], estimator = estimator[para.grid[ind,      "estimator"]], window.dims = window.blockpermutation[[para.grid[ind,      "size"]]][1, ], corr.block = FALSE)
# min       lq     mean  median       uq      max neval
# 59.41301 59.71785 59.82419 59.8192 59.93488 60.20918   100

# Genton
# lag 1
# Unit: seconds
# expr
# isotropy_blockpermutation(data = datas[[1]], lagmat = lags[[para.grid[ind,      "lag"]]], A = As[[para.grid[ind, "lag"]]], estimator = estimator[para.grid[ind,      "estimator"]], window.dims = window.blockpermutation[[para.grid[ind,      "size"]]][1, ], corr.block = FALSE)
# min       lq     mean   median       uq     max neval
# 21.61986 21.79088 21.84116 21.83406 21.89206 22.2421   100

# lag 2
# Unit: seconds
# expr
# isotropy_blockpermutation(data = datas[[1]], lagmat = lags[[para.grid[ind,      "lag"]]], A = As[[para.grid[ind, "lag"]]], estimator = estimator[para.grid[ind,      "estimator"]], window.dims = window.blockpermutation[[para.grid[ind,      "size"]]][1, ], corr.block = FALSE)
# min       lq     mean   median       uq      max neval
# 37.69514 37.85784 38.27917 38.00813 38.99116 39.31519   100

# lag 3
# Unit: seconds
# expr
# isotropy_blockpermutation(data = datas[[1]], lagmat = lags[[para.grid[ind,      "lag"]]], A = As[[para.grid[ind, "lag"]]], estimator = estimator[para.grid[ind,      "estimator"]], window.dims = window.blockpermutation[[para.grid[ind,      "size"]]][1, ], corr.block = FALSE)
# min       lq     mean   median       uq      max neval
# 69.03808 69.40069 69.49113 69.47558 69.62106 69.85337   100

# MCD
# lag 1
# Unit: seconds
# expr
# isotropy_blockpermutation(data = datas[[1]], lagmat = lags[[para.grid[ind,      "lag"]]], A = As[[para.grid[ind, "lag"]]], estimator = estimator[para.grid[ind,      "estimator"]], window.dims = window.blockpermutation[[para.grid[ind,      "size"]]][1, ], corr.block = FALSE)
# min       lq     mean   median       uq      max neval
# 24.03562 24.15783 24.21518 24.21842 24.26015 24.66413   100

# lag 2
# Unit: seconds
# expr
# isotropy_blockpermutation(data = datas[[1]], lagmat = lags[[para.grid[ind,      "lag"]]], A = As[[para.grid[ind, "lag"]]], estimator = estimator[para.grid[ind,      "estimator"]], window.dims = window.blockpermutation[[para.grid[ind,      "size"]]][1, ], corr.block = FALSE)
# min      lq     mean   median       uq      max neval
# 41.20129 41.4704 42.94472 41.60266 45.27984 45.71867   100

# lag 3
# Unit: seconds
# expr
# isotropy_blockpermutation(data = datas[[1]], lagmat = lags[[para.grid[ind,      "lag"]]], A = As[[para.grid[ind, "lag"]]], estimator = estimator[para.grid[ind,      "estimator"]], window.dims = window.blockpermutation[[para.grid[ind,      "size"]]][1, ], corr.block = FALSE)
# min       lq     mean   median       uq      max neval
# 72.66369 72.93366 73.05981 73.04748 73.17898 73.68235   100


## 2. 40 times 40
para.grid2 <- para.grid[which(para.grid$size == 2),]

## Subampling
para.grid2.sub <- para.grid2[which(para.grid1$method == 1),]
nrow(para.grid2.sub)
#[1] 9

# Matheron
# lag 1
# Unit: seconds
# expr
# isotropy_subsampling(data = datas[[1]], lagmat = lags[[para.grid[ind,      "lag"]]], A = As[[para.grid[ind, "lag"]]], estimator = estimator[para.grid[ind,      "estimator"]], window.dims = window.subsampling[para.grid[ind,      "size"], ], edge = TRUE)
# min       lq    mean  median       uq      max neval
# 26.56672 27.17356 27.5832 27.5411 27.85346 37.94405   100

# lag 2
# Unit: seconds
# expr
# isotropy_subsampling(data = datas[[1]], lagmat = lags[[para.grid[ind,      "lag"]]], A = As[[para.grid[ind, "lag"]]], estimator = estimator[para.grid[ind,      "estimator"]], window.dims = window.subsampling[para.grid[ind,      "size"], ], edge = TRUE)
# min       lq     mean   median       uq      max neval
# 25.05635 27.47657 28.46627 27.88119 29.17565 38.65957   100

# lag 3
# Unit: seconds
# expr
# isotropy_subsampling(data = datas[[1]], lagmat = lags[[para.grid[ind,      "lag"]]], A = As[[para.grid[ind, "lag"]]], estimator = estimator[para.grid[ind,      "estimator"]], window.dims = window.subsampling[para.grid[ind,      "size"], ], edge = TRUE)
# min       lq     mean   median       uq      max neval
# 28.2343 29.76543 31.10639 30.23837 30.57446 43.23387   100

# Genton
# lag 1
# Unit: seconds
# expr
# isotropy_subsampling(data = datas[[1]], lagmat = lags[[para.grid[ind,      "lag"]]], A = As[[para.grid[ind, "lag"]]], estimator = estimator[para.grid[ind,      "estimator"]], window.dims = window.subsampling[para.grid[ind,      "size"], ], edge = TRUE)
# min       lq     mean   median       uq      max neval
# 26.96153 27.03306 27.40924 27.44829 27.58377 37.72213   100

# lag 2
# Unit: seconds
# expr
# isotropy_subsampling(data = datas[[1]], lagmat = lags[[para.grid[ind,      "lag"]]], A = As[[para.grid[ind, "lag"]]], estimator = estimator[para.grid[ind,      "estimator"]], window.dims = window.subsampling[para.grid[ind,      "size"], ], edge = TRUE)
# min       lq     mean   median       uq      max neval
# 23.82389 27.13256 27.33647 27.78609 28.00978 38.43767   100

# lag 3
# Unit: seconds
# expr
# isotropy_subsampling(data = datas[[1]], lagmat = lags[[para.grid[ind,      "lag"]]], A = As[[para.grid[ind, "lag"]]], estimator = estimator[para.grid[ind,      "estimator"]], window.dims = window.subsampling[para.grid[ind,      "size"], ], edge = TRUE)
# min       lq    mean   median       uq      max neval
# 27.27809 28.03305 28.6016 28.20661 28.66673 39.31418   100

# MCD
# lag 1
# Unit: seconds
# expr
# isotropy_subsampling(data = datas[[1]], lagmat = lags[[para.grid[ind,      "lag"]]], A = As[[para.grid[ind, "lag"]]], estimator = estimator[para.grid[ind,      "estimator"]], window.dims = window.subsampling[para.grid[ind,      "size"], ], edge = TRUE)
# min      lq     mean   median       uq      max neval
# 29.96669 30.7239 30.99365 30.99873 31.01718 41.04103   100

# lag 2
# Unit: seconds
# expr
# isotropy_subsampling(data = datas[[1]], lagmat = lags[[para.grid[ind,      "lag"]]], A = As[[para.grid[ind, "lag"]]], estimator = estimator[para.grid[ind,      "estimator"]], window.dims = window.subsampling[para.grid[ind,      "size"], ], edge = TRUE)
# min       lq    mean   median       uq      max neval
# 39.84789 47.52947 47.4575 47.63005 47.75526 57.43872   100

# lag 3
# Unit: seconds
# expr
# isotropy_subsampling(data = datas[[1]], lagmat = lags[[para.grid[ind,      "lag"]]], A = As[[para.grid[ind, "lag"]]], estimator = estimator[para.grid[ind,      "estimator"]], window.dims = window.subsampling[para.grid[ind,      "size"], ], edge = TRUE)
# min       lq     mean   median       uq      max neval
# 42.24968 48.77988 48.90358 48.95966 49.14741 65.52118   100

## Blockpermutation
para.grid2.block <- para.grid2[which(para.grid2$method == 2),]
nrow(para.grid2.block)
#[1] 9

# Matheron
# lag 1
# [[1]]
# Unit: seconds
# expr
# isotropy_blockpermutation(data = datas[[1]], lagmat = lags[[para.grid[ind,      "lag"]]], A = As[[para.grid[ind, "lag"]]], estimator = estimator[para.grid[ind,      "estimator"]], window.dims = window.blockpermutation[[para.grid[ind,      "size"]]][1, ], corr.block = FALSE)
# min       lq     mean   median       uq    max neval
# 94.6923 96.20465 97.90572 97.16752 99.65665 106.65   100
# 
# [[2]]
# Unit: seconds
# expr
# isotropy_blockpermutation(data = datas[[1]], lagmat = lags[[para.grid[ind,      "lag"]]], A = As[[para.grid[ind, "lag"]]], estimator = estimator[para.grid[ind,      "estimator"]], window.dims = window.blockpermutation[[para.grid[ind,      "size"]]][2, ], corr.block = FALSE)
# min       lq     mean   median      uq      max neval
# 87.72233 90.46506 91.70159 90.98669 91.8544 97.51199   100
# 
# [[3]]
# Unit: seconds
# expr
# isotropy_blockpermutation(data = datas[[1]], lagmat = lags[[para.grid[ind,      "lag"]]], A = As[[para.grid[ind, "lag"]]], estimator = estimator[para.grid[ind,      "estimator"]], window.dims = window.blockpermutation[[para.grid[ind,      "size"]]][3, ], corr.block = FALSE)
# min       lq     mean  median       uq      max neval
# 84.86892 88.72021 89.69959 89.6716 90.62287 94.83175   100

# lag 2
# [[1]]
# Unit: seconds
# expr
# isotropy_blockpermutation(data = datas[[1]], lagmat = lags[[para.grid[ind,      "lag"]]], A = As[[para.grid[ind, "lag"]]], estimator = estimator[para.grid[ind,      "estimator"]], window.dims = window.blockpermutation[[para.grid[ind,      "size"]]][1, ], corr.block = FALSE)
# min       lq     mean  median       uq     max neval
# 138.1982 139.8396 141.6395 141.014 142.6756 158.844   100
# 
# [[2]]
# Unit: seconds
# expr
# isotropy_blockpermutation(data = datas[[1]], lagmat = lags[[para.grid[ind,      "lag"]]], A = As[[para.grid[ind, "lag"]]], estimator = estimator[para.grid[ind,      "estimator"]], window.dims = window.blockpermutation[[para.grid[ind,      "size"]]][2, ], corr.block = FALSE)
# min       lq     mean   median       uq      max neval
# 139.3084 141.2412 143.8838 143.0652 145.5419 155.5316   100
# 
# [[3]]
# Unit: seconds
# expr
# isotropy_blockpermutation(data = datas[[1]], lagmat = lags[[para.grid[ind,      "lag"]]], A = As[[para.grid[ind, "lag"]]], estimator = estimator[para.grid[ind,      "estimator"]], window.dims = window.blockpermutation[[para.grid[ind,      "size"]]][3, ], corr.block = FALSE)
# min       lq     mean  median       uq      max neval
# 140.9242 143.4593 146.1665 145.262 148.1749 158.0984   100

# lag 3
# [[1]]
# Unit: seconds
# expr
# isotropy_blockpermutation(data = datas[[1]], lagmat = lags[[para.grid[ind,      "lag"]]], A = As[[para.grid[ind, "lag"]]], estimator = estimator[para.grid[ind,      "estimator"]], window.dims = window.blockpermutation[[para.grid[ind,      "size"]]][1, ], corr.block = FALSE)
# min       lq     mean   median       uq     max neval
# 204.5343 207.5055 215.5999 209.3438 222.4422 239.321   100
# 
# [[2]]
# Unit: seconds
# expr
# isotropy_blockpermutation(data = datas[[1]], lagmat = lags[[para.grid[ind,      "lag"]]], A = As[[para.grid[ind, "lag"]]], estimator = estimator[para.grid[ind,      "estimator"]], window.dims = window.blockpermutation[[para.grid[ind,      "size"]]][2, ], corr.block = FALSE)
# min      lq     mean   median       uq      max neval
# 231.3092 257.369 265.4439 265.2403 278.3804 286.9242   100
# 
# [[3]]
# Unit: seconds
# expr
# isotropy_blockpermutation(data = datas[[1]], lagmat = lags[[para.grid[ind,      "lag"]]], A = As[[para.grid[ind, "lag"]]], estimator = estimator[para.grid[ind,      "estimator"]], window.dims = window.blockpermutation[[para.grid[ind,      "size"]]][3, ], corr.block = FALSE)
# min       lq    mean   median       uq      max neval
# 244.3054 258.4493 260.657 260.4009 261.9757 286.7263   100

# Genton
# lag 1
# [[1]]
# Unit: seconds
# expr
# isotropy_blockpermutation(data = datas[[1]], lagmat = lags[[para.grid[ind,      "lag"]]], A = As[[para.grid[ind, "lag"]]], estimator = estimator[para.grid[ind,      "estimator"]], window.dims = window.blockpermutation[[para.grid[ind,      "size"]]][1, ], corr.block = FALSE)
# min       lq     mean   median       uq      max neval
# 100.8779 102.4171 105.0594 104.1956 107.2411 117.3879   100
# 
# [[2]]
# Unit: seconds
# expr
# isotropy_blockpermutation(data = datas[[1]], lagmat = lags[[para.grid[ind,      "lag"]]], A = As[[para.grid[ind, "lag"]]], estimator = estimator[para.grid[ind,      "estimator"]], window.dims = window.blockpermutation[[para.grid[ind,      "size"]]][2, ], corr.block = FALSE)
# min       lq     mean   median      uq      max neval
# 95.11661 96.76792 99.89749 99.72636 102.361 108.3459   100
# 
# [[3]]
# Unit: seconds
# expr
# isotropy_blockpermutation(data = datas[[1]], lagmat = lags[[para.grid[ind,      "lag"]]], A = As[[para.grid[ind, "lag"]]], estimator = estimator[para.grid[ind,      "estimator"]], window.dims = window.blockpermutation[[para.grid[ind,      "size"]]][3, ], corr.block = FALSE)
# min       lq     mean   median       uq      max neval
# 93.5827 96.63911 98.10373 97.63685 99.22819 104.7439   100

# lag 2
# [[1]]
# Unit: seconds
# expr
# isotropy_blockpermutation(data = datas[[1]], lagmat = lags[[para.grid[ind,      "lag"]]], A = As[[para.grid[ind, "lag"]]], estimator = estimator[para.grid[ind,      "estimator"]], window.dims = window.blockpermutation[[para.grid[ind,      "size"]]][1, ], corr.block = FALSE)
# min       lq     mean  median       uq      max neval
# 129.7486 131.7493 132.7475 132.379 133.6898 148.9582   100
# 
# [[2]]
# Unit: seconds
# expr
# isotropy_blockpermutation(data = datas[[1]], lagmat = lags[[para.grid[ind,      "lag"]]], A = As[[para.grid[ind, "lag"]]], estimator = estimator[para.grid[ind,      "estimator"]], window.dims = window.blockpermutation[[para.grid[ind,      "size"]]][2, ], corr.block = FALSE)
# min       lq     mean   median       uq      max neval
# 132.0182 133.5638 134.2483 134.3171 134.9301 137.1762   100
# 
# [[3]]
# Unit: seconds
# expr
# isotropy_blockpermutation(data = datas[[1]], lagmat = lags[[para.grid[ind,      "lag"]]], A = As[[para.grid[ind, "lag"]]], estimator = estimator[para.grid[ind,      "estimator"]], window.dims = window.blockpermutation[[para.grid[ind,      "size"]]][3, ], corr.block = FALSE)
# min       lq     mean median       uq      max neval
# 133.2796 134.6045 135.4132 135.24 135.9808 139.3237   100

# lag 3
# [[1]]
# Unit: seconds
# expr
# isotropy_blockpermutation(data = datas[[1]], lagmat = lags[[para.grid[ind,      "lag"]]], A = As[[para.grid[ind, "lag"]]], estimator = estimator[para.grid[ind,      "estimator"]], window.dims = window.blockpermutation[[para.grid[ind,      "size"]]][1, ], corr.block = FALSE)
# min       lq    mean   median       uq      max neval
# 211.8114 213.8628 215.188 214.4145 215.3162 238.0787   100
# 
# [[2]]
# Unit: seconds
# expr
# isotropy_blockpermutation(data = datas[[1]], lagmat = lags[[para.grid[ind,      "lag"]]], A = As[[para.grid[ind, "lag"]]], estimator = estimator[para.grid[ind,      "estimator"]], window.dims = window.blockpermutation[[para.grid[ind,      "size"]]][2, ], corr.block = FALSE)
# min       lq     mean   median       uq      max neval
# 241.672 245.2374 256.6011 254.1359 265.5966 283.3743   100
# 
# [[3]]
# Unit: seconds
# expr
# isotropy_blockpermutation(data = datas[[1]], lagmat = lags[[para.grid[ind,      "lag"]]], A = As[[para.grid[ind, "lag"]]], estimator = estimator[para.grid[ind,      "estimator"]], window.dims = window.blockpermutation[[para.grid[ind,      "size"]]][3, ], corr.block = FALSE)
# min       lq     mean   median       uq      max neval
# 252.3798 256.4098 263.0681 261.6867 269.9338 274.7048   100


# MCD
# lag 1
# [[1]]
# Unit: seconds
# expr
# isotropy_blockpermutation(data = datas[[1]], lagmat = lags[[para.grid[ind,      "lag"]]], A = As[[para.grid[ind, "lag"]]], estimator = estimator[para.grid[ind,      "estimator"]], window.dims = window.blockpermutation[[para.grid[ind,      "size"]]][1, ], corr.block = FALSE)
# min       lq     mean   median       uq      max neval
# 89.57649 98.34298 98.98201 99.23888 100.4572 108.2438   100
# 
# [[2]]
# Unit: seconds
# expr
# isotropy_blockpermutation(data = datas[[1]], lagmat = lags[[para.grid[ind,      "lag"]]], A = As[[para.grid[ind, "lag"]]], estimator = estimator[para.grid[ind,      "estimator"]], window.dims = window.blockpermutation[[para.grid[ind,      "size"]]][2, ], corr.block = FALSE)
# min       lq     mean   median       uq      max neval
# 84.26201 91.30857 91.86091 91.88899 92.40804 99.04361   100
# 
# [[3]]
# Unit: seconds
# expr
# isotropy_blockpermutation(data = datas[[1]], lagmat = lags[[para.grid[ind,      "lag"]]], A = As[[para.grid[ind, "lag"]]], estimator = estimator[para.grid[ind,      "estimator"]], window.dims = window.blockpermutation[[para.grid[ind,      "size"]]][3, ], corr.block = FALSE)
# min       lq     mean   median       uq     max neval
# 86.62673 90.14665 91.43077 90.65766 91.79623 98.1642   100

# lag 2
# [[1]]
# Unit: seconds
# expr
# isotropy_blockpermutation(data = datas[[1]], lagmat = lags[[para.grid[ind,      "lag"]]], A = As[[para.grid[ind, "lag"]]], estimator = estimator[para.grid[ind,      "estimator"]], window.dims = window.blockpermutation[[para.grid[ind,      "size"]]][1, ], corr.block = FALSE)
# min      lq     mean   median       uq      max neval
# 134.201 135.588 136.8727 136.4634 137.3999 152.1335   100
# 
# [[2]]
# Unit: seconds
# expr
# isotropy_blockpermutation(data = datas[[1]], lagmat = lags[[para.grid[ind,      "lag"]]], A = As[[para.grid[ind, "lag"]]], estimator = estimator[para.grid[ind,      "estimator"]], window.dims = window.blockpermutation[[para.grid[ind,      "size"]]][2, ], corr.block = FALSE)
# min      lq    mean   median       uq      max neval
# 135.9072 137.414 137.974 137.9002 138.3417 142.4844   100
# 
# [[3]]
# Unit: seconds
# expr
# isotropy_blockpermutation(data = datas[[1]], lagmat = lags[[para.grid[ind,      "lag"]]], A = As[[para.grid[ind, "lag"]]], estimator = estimator[para.grid[ind,      "estimator"]], window.dims = window.blockpermutation[[para.grid[ind,      "size"]]][3, ], corr.block = FALSE)
# min       lq     mean   median       uq      max neval
# 137.1764 138.8455 139.5467 139.6601 140.1817 144.3218   100

# lag 3
# [[1]]
# Unit: seconds
# expr
# isotropy_blockpermutation(data = datas[[1]], lagmat = lags[[para.grid[ind,      "lag"]]], A = As[[para.grid[ind, "lag"]]], estimator = estimator[para.grid[ind,      "estimator"]], window.dims = window.blockpermutation[[para.grid[ind,      "size"]]][1, ], corr.block = FALSE)
# min       lq     mean   median       uq      max neval
# 198.5297 202.2416 205.0531 204.8939 207.9596 222.2736   100
# 
# [[2]]
# Unit: seconds
# expr
# isotropy_blockpermutation(data = datas[[1]], lagmat = lags[[para.grid[ind,      "lag"]]], A = As[[para.grid[ind, "lag"]]], estimator = estimator[para.grid[ind,      "estimator"]], window.dims = window.blockpermutation[[para.grid[ind,      "size"]]][2, ], corr.block = FALSE)
# min       lq     mean  median      uq      max neval
# 222.4803 223.8139 224.8792 224.787 225.812 227.6893   100
# 
# [[3]]
# Unit: seconds
# expr
# isotropy_blockpermutation(data = datas[[1]], lagmat = lags[[para.grid[ind,      "lag"]]], A = As[[para.grid[ind, "lag"]]], estimator = estimator[para.grid[ind,      "estimator"]], window.dims = window.blockpermutation[[para.grid[ind,      "size"]]][3, ], corr.block = FALSE)
# min       lq     mean   median       uq      max neval
# 231.8339 236.1655 239.7716 240.8093 242.9954 245.9798   100




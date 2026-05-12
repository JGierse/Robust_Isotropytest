################################################################################
## Parameter values for the simulations presented in the paper                ##
################################################################################


## grid size
size <- rbind(c("nx" = 24, "ny" = 24),
              c("nx" = 40, "ny" = 40))

## (an)isotropy parameters
aniso <- cbind(c(0,2), c(0, sqrt(2)), c(0, 2), c(pi/4, sqrt(2)), c(pi/4, 2))


## variogrammodel
model <- "spherical"
nugget <- 0
sill <- 1
range <- c(2, 5, 8)


## lag set: relavant for Test/Simulation
lags <- list(rbind(c(1,0), c(0,1)),
             rbind(c(1,0), c(0,1), c(1,1), c(1, -1)),
             rbind(c(1,0), c(0,1), c(1,1), c(1,-1), c(2,1), c(-1, 2), c(1, 2), c(-2, 1)))

## contrast matrices
As <- list(rbind(c(1, -1)),
           rbind(c(1, -1, 0, 0), c(0, 0, 1, -1)),
           rbind(c(1, -1, 0, 0, 0, 0, 0, 0), c(0, 0, 1, -1, 0, 0, 0, 0), c(0, 0, 0, 0, 1, -1, 0, 0), c(0, 0, 0, 0, 0, 0, 1, -1)))

## variogram estimators
estimator <- c("Matheron", "Genton", "MCD")

## Test used 
Method <- c("Subsampling", "Blockpermutation")


## Subsampling: subsample size
window.subsampling <- rbind(c(5,5),
                            c(6,6))

## Blockpermutation: block size
window.blockpermutation <- list(c(6,6),
                                rbind(c(5,5), c(8,8), c(10,10)))


## Parameters outliers
dists <- rbind(c(0,5),
               c(5,1))

amount <- c(0.1, 0.2)



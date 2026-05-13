################################################################################
## Application of the test to the Satelite Data (Section 5)                   ##
################################################################################

##############################
## Load functions and packages
##############################

# Necessary packages for preparing the satellite data
library(gdalcubes)
library(magrittr)
library(xts)
library(magick)

# Necessary for the variogram estimation
library(tidyverse)
library(RobVario)


# graphs
library(ggplot2)
library(tidyr)
library(dplyr)
library(patchwork)

####################################
## Preparation of the satellite data
####################################

## Load the satellite data
## Data can be downloaded from https://earthexplorer.usgs.gov/
## the used data here have the following dimensions in space and time:
# A data cube view object
# 
# Dimensions:
#   low              high count pixel_size
# t        2013-04-12        2019-11-30  2424        P1D
# y -1096876.57818711 -993876.578187112   206        500
# x -7370181.54841633 -7235181.54841633   270        500
# 
# SRS: "EPSG:3857"
# Temporal aggregation method: "first"
# Spatial resampling method: "near"
# IMAGE_DIR = "C:/Users/jkoenig/sciebo2/L8_cropped"  # Path of the data 

# load the structure of the data
col <- create_image_collection(list.files(IMAGE_DIR, recursive = TRUE, pattern=".tif", full.names  = TRUE), "L8_SR")

## Selection of pixels based on quality band
# only pixel classified as clear
L8.clear_mask = image_mask("PIXEL_QA", values=c(322, 386, 834, 898, 1346, 324, 388, 836, 900, 1348), invert = TRUE)

# Define what the data set should look like
v = cube_view(srs="EPSG:3857", extent=col, dx=500, dy=500, dt="P1D")

# Transform data into the required format and select the required spectral bands
comp_mask <- raster_cube(col, v, mask = L8.clear_mask) %>% select_bands((names(raster_cube(col,v))[2:8]))
comp <- raster_cube(col, v) %>% select_bands((names(raster_cube(col,v))[2:8]))
v.date <-cube_view(view = v, extent=list(left = -7331535-400*30, right = -7331535+400*30, 
                                         bottom = -999642-400*30, top = -999642 + 150*30, t0 = "2016-10-20", t1 = "2016-10-20"), dx=30, dy=30, dt="P1D") 
comp.date <- raster_cube(col, v.date) %>% select_bands(c("B02", "B03", "B04"))


# Function for calculating the NDVI for a subregion
NDVI_sub <- function(pixel_x, pixel_y, time, mask = NULL){
  delta = 30
  v.sub = cube_view(view = v, extent=list(left=pixel_x-(20 * delta), right=pixel_x+(40 * delta),
                                          bottom=pixel_y-(20 * delta),top=pixel_y+(40 * delta),
                                          t0 = time, t1 = time), dx=delta, dy=delta, dt="P1D")
  raster_cube(col, v.sub, mask=mask) %>%
    select_bands(c("B04", "B05")) %>%
    apply_pixel("(B05-B04)/(B05+B04)", "NDVI")  -> tseries_plot
  
  raster_cube(col, v.sub, mask=mask) %>%
    select_bands(c("B04", "B05")) %>%
    apply_pixel("(B05-B04)/(B05+B04)", "NDVI") %>%
    as_array() %>% drop -> tseries
  
  return(list(tseries = tseries, tseries_plot = tseries_plot))
}


#############
## 60x60 grid 
#############

# all data in the subregion
delta = 30
v.sub1 = cube_view(view = v, extent=list(left=-7331535-20*delta, right=-7331535+40*delta, 
                                         bottom=-999642-20*delta,top=-999642+40*delta,
                                         t0 = "2016-10-20", t1 = "2016-10-20"), dx=delta, dy=delta, dt="P1D") 
sub1 <- raster_cube(col, v.sub1) %>% select_bands(c("B02", "B03", "B04"))

# only clear data in the subregion
sub.clear1 <- raster_cube(col, v.sub1, mask = L8.clear_mask) %>% select_bands(c("B02", "B03", "B04"))

# calculate the pixlewise NDVI for the subregion
reg1 <-  NDVI_sub(-7331535, -999642, "2016-10-20")

# extract the NDVI Index
reg1.data = reg1$tseries

# rename the y-cordinates
rownames(reg1.data) <- sapply(60:1, function(x) paste0("y", x))

# add the coordinates
reg1.data <- cbind(sapply(60:1, function(x) paste0("y", x)), reg1.data)

# rename the x-coordinates
colnames(reg1.data) <- c("y", sapply(1:60, function(x) paste0("x", x)))

# data set with a row for each data point
reg1.data.long <- as_tibble(reg1.data) %>% pivot_longer(cols = starts_with("x"), values_to = "data") %>%  
  mutate(data = as.numeric(data), ID = row_number()) %>% rename(x = name) %>% 
  mutate(x = as.numeric(sub("x", "", x)), y = as.numeric(sub("y", "", y))) 

# estimate the standard deviation of the data with an robust estimator
sd.all <- mad(as.vector(reg1$tseries))

# standardize the data with the standard deviation fo the complete data 
reg1.data.long <- reg1.data.long %>% mutate(data.sd = data/sd.all)

# Graphic of the NDVI of the subregion (Figure 9)
pdf("Graphs/NDVI.pdf")
plot(reg1$tseries_plot, key.pos = 4, breaks = seq(0.6, 0.9, 0.02), nbreaks = 16, col = rev(hcl.colors(15, palette = "Greens 2")))
dev.off()



#############
## 30x30 grid: Divide into 4 subgrids
#############

#############
#### 1. grid: bottom left

# all data bottom left
v.sub1.1 = cube_view(view = v, extent=list(left=-7331535-20*delta, right=-7331535+10*delta, 
                                           bottom=-999642-20*delta,top=-999642+10*delta,
                                           t0 = "2016-10-20", t1 = "2016-10-20"), dx=delta, dy=delta, dt="P1D") 
sub1.1 <- raster_cube(col, v.sub1.1) %>% select_bands(c("B02", "B03", "B04"))

# Function for calculating the NDVI in the bottom left
NDVI_sub.1 <- function(pixel_x, pixel_y, time, mask = NULL){
  delta = 30
  v.sub = cube_view(view = v, extent=list(left=pixel_x-(20 * delta), right=pixel_x+(10 * delta),
                                          bottom=pixel_y-(20 * delta),top=pixel_y+(10 * delta),
                                          t0 = time, t1 = time), dx=delta, dy=delta, dt="P1D")
  raster_cube(col, v.sub, mask=mask) %>%
    select_bands(c("B04", "B05")) %>%
    apply_pixel("(B05-B04)/(B05+B04)", "NDVI")  -> tseries_plot
  
  raster_cube(col, v.sub, mask=mask) %>%
    select_bands(c("B04", "B05")) %>%
    apply_pixel("(B05-B04)/(B05+B04)", "NDVI") %>%
    as_array() %>% drop -> tseries
  
  return(list(tseries = tseries, tseries_plot = tseries_plot))
}

reg1.1 <-  NDVI_sub.1(-7331535, -999642, "2016-10-20")

# extract the NDVI Index
reg1.1.data = reg1.1$tseries

# rename the y-cordinates
rownames(reg1.1.data) <- sapply(30:1, function(x) paste0("y", x))

# add the coordinates
reg1.1.data <- cbind(sapply(30:1, function(x) paste0("y", x)), reg1.1.data)

# rename the x-coordinates
colnames(reg1.1.data) <- c("y", sapply(1:30, function(x) paste0("x", x)))

# data set with a row for each data point
reg1.1.data.long <- as_tibble(reg1.1.data) %>% pivot_longer(cols = starts_with("x"), values_to = "data") %>%  
  mutate(data = as.numeric(data), ID = row_number()) %>% rename(x = name) %>% 
  mutate(x = as.numeric(sub("x", "", x)), y = as.numeric(sub("y", "", y))) 

# estimate the standard deviation of the data with an robust estimator
sd.all.1 <- mad(as.vector(reg1.1$tseries))

# standardaize the data with the standard deviation fo the complete data 
reg1.1.data.long <- reg1.1.data.long %>% mutate(data.sd = data/sd.all)

pdf("Graphs/NDVI_UL.pdf")
plot(reg1.1$tseries_plot, key.pos = 4, breaks = seq(0.6, 0.9, 0.02), nbreaks = 16, col = rev(hcl.colors(15, palette = "Greens 2")))
dev.off()

#############
#### 2. grid: top left

# all data top left
v.sub1.2 = cube_view(view = v, extent=list(left=-7331535-20*delta, right=-7331535+10*delta, 
                                           bottom=-999642+10*delta,top=-999642+40*delta,
                                           t0 = "2016-10-20", t1 = "2016-10-20"), dx=delta, dy=delta, dt="P1D") 
sub1.2 <- raster_cube(col, v.sub1.2) %>% select_bands(c("B02", "B03", "B04"))

# Function for calculating the NDVI in the top left
NDVI_sub.2 <- function(pixel_x, pixel_y, time, mask = NULL){
  delta = 30
  v.sub = cube_view(view = v, extent=list(left=pixel_x-(20 * delta), right=pixel_x+(10 * delta),
                                          bottom=pixel_y+(10 * delta),top=pixel_y+(40 * delta),
                                          t0 = time, t1 = time), dx=delta, dy=delta, dt="P1D")
  raster_cube(col, v.sub, mask=mask) %>%
    select_bands(c("B04", "B05")) %>%
    apply_pixel("(B05-B04)/(B05+B04)", "NDVI")  -> tseries_plot
  
  raster_cube(col, v.sub, mask=mask) %>%
    select_bands(c("B04", "B05")) %>%
    apply_pixel("(B05-B04)/(B05+B04)", "NDVI") %>%
    as_array() %>% drop -> tseries
  
  return(list(tseries = tseries, tseries_plot = tseries_plot))
}

reg1.2 <-  NDVI_sub.2(-7331535, -999642, "2016-10-20")

# extract the NDVI Index
reg1.2.data = reg1.2$tseries

# rename the y-cordinates
rownames(reg1.2.data) <- sapply(30:1, function(x) paste0("y", x))

# add the coordinates
reg1.2.data <- cbind(sapply(30:1, function(x) paste0("y", x)), reg1.2.data)

# rename the x-coordinates
colnames(reg1.2.data) <- c("y", sapply(1:30, function(x) paste0("x", x)))

# data set with a row for each data point
reg1.2.data.long <- as_tibble(reg1.2.data) %>% pivot_longer(cols = starts_with("x"), values_to = "data") %>%  
  mutate(data = as.numeric(data), ID = row_number()) %>% rename(x = name) %>% 
  mutate(x = as.numeric(sub("x", "", x)), y = as.numeric(sub("y", "", y))) 

# estimate the standard deviation of the data with an robust estimator
sd.all.2 <- mad(as.vector(reg1.2$tseries))

# standardaize the data with the standard deviation fo the complete data 
reg1.2.data.long <- reg1.2.data.long %>% mutate(data.sd = data/sd.all)

pdf("Graphs/NDVI_OL.pdf")
plot(reg1.2$tseries_plot, key.pos = 4, breaks = seq(0.6, 0.9, 0.02), nbreaks = 16, col = rev(hcl.colors(15, palette = "Greens 2")))
dev.off()

#############
#### 3. grid: bottom right 

# all data bottom right
v.sub1.3 = cube_view(view = v, extent=list(left=-7331535+10*delta, right=-7331535+40*delta, 
                                           bottom=-999642-20*delta,top=-999642+10*delta,
                                           t0 = "2016-10-20", t1 = "2016-10-20"), dx=delta, dy=delta, dt="P1D") 
sub1.3 <- raster_cube(col, v.sub1.3) %>% select_bands(c("B02", "B03", "B04"))

# Function for calculating the NDVI in the bottom right
NDVI_sub.3 <- function(pixel_x, pixel_y, time, mask = NULL){
  delta = 30
  v.sub = cube_view(view = v, extent=list(left=pixel_x+(10 * delta), right=pixel_x+(40 * delta),
                                          bottom=pixel_y-(20 * delta),top=pixel_y+(10 * delta),
                                          t0 = time, t1 = time), dx=delta, dy=delta, dt="P1D")
  raster_cube(col, v.sub, mask=mask) %>%
    select_bands(c("B04", "B05")) %>%
    apply_pixel("(B05-B04)/(B05+B04)", "NDVI")  -> tseries_plot
  
  raster_cube(col, v.sub, mask=mask) %>%
    select_bands(c("B04", "B05")) %>%
    apply_pixel("(B05-B04)/(B05+B04)", "NDVI") %>%
    as_array() %>% drop -> tseries
  
  return(list(tseries = tseries, tseries_plot = tseries_plot))
}

reg1.3 <-  NDVI_sub.3(-7331535, -999642, "2016-10-20")

# extract the NDVI Index
reg1.3.data = reg1.3$tseries

# rename the y-cordinates
rownames(reg1.3.data) <- sapply(30:1, function(x) paste0("y", x))

# add the coordinates
reg1.3.data <- cbind(sapply(30:1, function(x) paste0("y", x)), reg1.3.data)

# rename the x-coordinates
colnames(reg1.3.data) <- c("y", sapply(1:30, function(x) paste0("x", x)))

# data set with a row for each data point
reg1.3.data.long <- as_tibble(reg1.3.data) %>% pivot_longer(cols = starts_with("x"), values_to = "data") %>%  
  mutate(data = as.numeric(data), ID = row_number()) %>% rename(x = name) %>% 
  mutate(x = as.numeric(sub("x", "", x)), y = as.numeric(sub("y", "", y))) 

# estimate the standard deviation of the data with an robust estimator
sd.all.3 <- mad(as.vector(reg1.3$tseries))

# standardaize the data with the standard deviation fo the complete data 
reg1.3.data.long <- reg1.3.data.long %>% mutate(data.sd = data/sd.all)

pdf("Graphs/NDVI_UR.pdf")
plot(reg1.3$tseries_plot, key.pos = 4, breaks = seq(0.6, 0.9, 0.02), nbreaks = 16, col = rev(hcl.colors(15, palette = "Greens 2")))
dev.off()

#############
#### 4. grid: top right 

# all data top right
v.sub1.4 = cube_view(view = v, extent=list(left=-7331535+10*delta, right=-7331535+40*delta, 
                                           bottom=-999642+10*delta,top=-999642+40*delta,
                                           t0 = "2016-10-20", t1 = "2016-10-20"), dx=delta, dy=delta, dt="P1D") 
sub1.4 <- raster_cube(col, v.sub1.4) %>% select_bands(c("B02", "B03", "B04"))

# Function for calculating the NDVI in the top right
NDVI_sub.4 <- function(pixel_x, pixel_y, time, mask = NULL){
  delta = 30
  v.sub = cube_view(view = v, extent=list(left=pixel_x+(10 * delta), right=pixel_x+(40 * delta),
                                          bottom=pixel_y+(10 * delta),top=pixel_y+(40 * delta),
                                          t0 = time, t1 = time), dx=delta, dy=delta, dt="P1D")
  raster_cube(col, v.sub, mask=mask) %>%
    select_bands(c("B04", "B05")) %>%
    apply_pixel("(B05-B04)/(B05+B04)", "NDVI")  -> tseries_plot
  
  raster_cube(col, v.sub, mask=mask) %>%
    select_bands(c("B04", "B05")) %>%
    apply_pixel("(B05-B04)/(B05+B04)", "NDVI") %>%
    as_array() %>% drop -> tseries
  
  return(list(tseries = tseries, tseries_plot = tseries_plot))
}

reg1.4 <-  NDVI_sub.4(-7331535, -999642, "2016-10-20")

# extract the NDVI Index
reg1.4.data = reg1.4$tseries

# rename the y-cordinates
rownames(reg1.4.data) <- sapply(30:1, function(x) paste0("y", x))

# add the coordinates
reg1.4.data <- cbind(sapply(30:1, function(x) paste0("y", x)), reg1.4.data)

# rename the x-coordinates
colnames(reg1.4.data) <- c("y", sapply(1:30, function(x) paste0("x", x)))

# data set with a row for each data point
reg1.4.data.long <- as_tibble(reg1.4.data) %>% pivot_longer(cols = starts_with("x"), values_to = "data") %>%  
  mutate(data = as.numeric(data), ID = row_number()) %>% rename(x = name) %>% 
  mutate(x = as.numeric(sub("x", "", x)), y = as.numeric(sub("y", "", y))) 

# estimate the standard deviation of the data with an robust estimator
sd.all.4 <- mad(as.vector(reg1.4$tseries))

# standardaize the data with the standard deviation fo the complete data 
reg1.4.data.long <- reg1.4.data.long %>% mutate(data.sd = data/sd.all)

pdf("Graphs/NDVI_OR.pdf")
plot(reg1.4$tseries_plot, key.pos = 4, breaks = seq(0.6, 0.9, 0.02), nbreaks = 16, col = rev(hcl.colors(15, palette = "Greens 2")))
dev.off()

#######################
### Variogramestimation: for the lags of the test
#######################

data <- as.data.frame(cbind(reg1.data.long$data.sd, reg1.data.long$x, reg1.data.long$y))

varog <- variogram_est_general(data, h = rbind(c(1,0), c(0,1), c(1,1), c(1,-1), c(2,1), c(-1, 2), c(1, 2), c(-2, 1)),
                               estimator = c("Matheron", "Genton", "MCD.diff"))

cbind(varog$`0,1`, varog$`1,0`)
# lag.x lag.y    n variogram estimator lag.x lag.y    n variogram estimator
# 1     0     1 3540 2.6743634  matheron     1     0 3540 3.1834009  matheron
# 2     0     1 3540 0.5605953    genton     1     0 3540 0.6730721    genton
# 3     0     1 3540 0.4547845  mcd.diff     1     0 3540 0.5555432  mcd.diff

cbind(varog$`1,-1`, varog$`1,1`)
# lag.x lag.y    n variogram estimator lag.x lag.y    n variogram estimator
# 1     1    -1 3481 4.4829499  matheron     1     1 3481 3.0606587  matheron
# 2     1    -1 3481 0.9321184    genton     1     1 3481 0.9007944    genton
# 3     1    -1 3481 0.7446536  mcd.diff     1     1 3481 0.7594240  mcd.diff

## Frage sind die hier richtig?!
cbind(varog$`0.5,1`, varog$`-1,0.5`)
# lag.x lag.y    n variogram estimator lag.x lag.y    n variogram estimator
# 1     1     2 3422 4.2051347  matheron    -2     1 3422  5.673798  matheron
# 2     1     2 3422 1.1944042    genton    -2     1 3422  1.342198    genton
# 3     1     2 3422 0.9773029  mcd.diff    -2     1 3422  1.036118  mcd.diff

cbind(varog$`1,0.5`, varog$`-0.5,1`)
# lag.x lag.y    n variogram estimator lag.x lag.y    n variogram estimator
# 1     2     1 3422  4.597875  matheron    -1     2 3422 5.2284805  matheron
# 2     2     1 3422  1.306475    genton    -1     2 3422 1.1774949    genton
# 3     2     1 3422  1.029415  mcd.diff    -1     2 3422 0.9147454  mcd.diff


####################
### Test of isotropy
####################

####################
### Blockpermutation

## complete data set
test.block1 <- isotropy_blockpermutation(data, lagmat =  rbind(c(1,0), c(0,1)), A =  matrix(c(1, -1), nrow = 1))

test.block2 <- isotropy_blockpermutation(data)

test.block3 <- isotropy_blockpermutation(data, lagmat =  rbind(c(1,0), c(0,1), c(1,1), c(1,-1), c(2,1), c(-1, 2), c(1, 2), c(-2, 1)),
                                         A =  rbind(c(1, -1, 0, 0, 0, 0, 0, 0), c(0, 0, 1, -1, 0, 0, 0, 0), c(0, 0, 0, 0, 1, -1, 0, 0), c(0, 0, 0, 0, 0, 0, 1, -1)))

## 1. grid: bottom left 
data.1 <- as.data.frame(cbind(reg1.1.data.long$data.sd, reg1.1.data.long$x, reg1.1.data.long$y))

test.block1.1 <- isotropy_blockpermutation(data.1, lagmat =  rbind(c(1,0), c(0,1)), A =  matrix(c(1, -1), nrow = 1))

test.block2.1 <- isotropy_blockpermutation(data.1)

test.block3.1 <- isotropy_blockpermutation(data.1, lagmat =  rbind(c(1,0), c(0,1), c(1,1), c(1,-1), c(2,1), c(-1, 2), c(1, 2), c(-2, 1)),
                                           A =  rbind(c(1, -1, 0, 0, 0, 0, 0, 0), c(0, 0, 1, -1, 0, 0, 0, 0), c(0, 0, 0, 0, 1, -1, 0, 0), c(0, 0, 0, 0, 0, 0, 1, -1)))

## 2 grid: top left
data.2 <- as.data.frame(cbind(reg1.2.data.long$data.sd, reg1.2.data.long$x, reg1.2.data.long$y))

test.block1.2 <- isotropy_blockpermutation(data.2, lagmat =  rbind(c(1,0), c(0,1)), A =  matrix(c(1, -1), nrow = 1))

test.block2.2 <- isotropy_blockpermutation(data.2)

test.block3.2 <- isotropy_blockpermutation(data.2, lagmat =  rbind(c(1,0), c(0,1), c(1,1), c(1,-1), c(2,1), c(-1, 2), c(1, 2), c(-2, 1)),
                                           A =  rbind(c(1, -1, 0, 0, 0, 0, 0, 0), c(0, 0, 1, -1, 0, 0, 0, 0), c(0, 0, 0, 0, 1, -1, 0, 0), c(0, 0, 0, 0, 0, 0, 1, -1)))

## 3. grid: bottom right
data.3 <- as.data.frame(cbind(reg1.3.data.long$data.sd, reg1.3.data.long$x, reg1.3.data.long$y))

test.block1.3 <- isotropy_blockpermutation(data.3, lagmat =  rbind(c(1,0), c(0,1)), A =  matrix(c(1, -1), nrow = 1))

test.block2.3 <- isotropy_blockpermutation(data.3)

test.block3.3 <- isotropy_blockpermutation(data.3, lagmat =  rbind(c(1,0), c(0,1), c(1,1), c(1,-1), c(2,1), c(-1, 2), c(1, 2), c(-2, 1)),
                                           A =  rbind(c(1, -1, 0, 0, 0, 0, 0, 0), c(0, 0, 1, -1, 0, 0, 0, 0), c(0, 0, 0, 0, 1, -1, 0, 0), c(0, 0, 0, 0, 0, 0, 1, -1)))


## 4. grid: top right
data.4 <- as.data.frame(cbind(reg1.4.data.long$data.sd, reg1.4.data.long$x, reg1.4.data.long$y))

test.block1.4 <- isotropy_blockpermutation(data.4, lagmat =  rbind(c(1,0), c(0,1)), A =  matrix(c(1, -1), nrow = 1))

test.block2.4 <- isotropy_blockpermutation(data.4)

test.block3.4 <- isotropy_blockpermutation(data.4, lagmat =  rbind(c(1,0), c(0,1), c(1,1), c(1,-1), c(2,1), c(-1, 2), c(1, 2), c(-2, 1)),
                                           A =  rbind(c(1, -1, 0, 0, 0, 0, 0, 0), c(0, 0, 1, -1, 0, 0, 0, 0), c(0, 0, 0, 0, 1, -1, 0, 0), c(0, 0, 0, 0, 0, 0, 1, -1)))


###############
### Subsampling

## complete data set
test.sub1 <- isotropy_subsampling(data, lagmat =  rbind(c(1,0), c(0,1)), A =  matrix(c(1, -1), nrow = 1), window.dims = c(8,8))

test.sub2 <- isotropy_subsampling(data, window.dims = c(8,8))


test.sub3 <- isotropy_subsampling(data, lagmat =  rbind(c(1,0), c(0,1), c(1,1), c(1,-1), c(2,1), c(-1, 2), c(1, 2), c(-2, 1)),
                                  A =  rbind(c(1, -1, 0, 0, 0, 0, 0, 0), c(0, 0, 1, -1, 0, 0, 0, 0), c(0, 0, 0, 0, 1, -1, 0, 0), c(0, 0, 0, 0, 0, 0, 1, -1)),
                                  window.dims = c(8,8))

## 1. grid: bottom left 
test.sub1.1 <- isotropy_subsampling(data.1, lagmat =  rbind(c(1,0), c(0,1)), A =  matrix(c(1, -1), nrow = 1))

test.sub2.1 <- isotropy_subsampling(data.1)

test.sub3.1 <- isotropy_subsampling(data.1, lagmat =  rbind(c(1,0), c(0,1), c(1,1), c(1,-1), c(2,1), c(-1, 2), c(1, 2), c(-2, 1)),
                                    A =  rbind(c(1, -1, 0, 0, 0, 0, 0, 0), c(0, 0, 1, -1, 0, 0, 0, 0), c(0, 0, 0, 0, 1, -1, 0, 0), c(0, 0, 0, 0, 0, 0, 1, -1)))

## 2 grid: top left
test.sub1.2 <- isotropy_subsampling(data.2, lagmat =  rbind(c(1,0), c(0,1)), A =  matrix(c(1, -1), nrow = 1))

test.sub2.2 <- isotropy_subsampling(data.2)

test.sub3.2 <- isotropy_subsampling(data.2, lagmat =  rbind(c(1,0), c(0,1), c(1,1), c(1,-1), c(2,1), c(-1, 2), c(1, 2), c(-2, 1)),
                                    A =  rbind(c(1, -1, 0, 0, 0, 0, 0, 0), c(0, 0, 1, -1, 0, 0, 0, 0), c(0, 0, 0, 0, 1, -1, 0, 0), c(0, 0, 0, 0, 0, 0, 1, -1)))

## 3. grid: bottom right
test.sub1.3 <- isotropy_subsampling(data.3, lagmat =  rbind(c(1,0), c(0,1)), A =  matrix(c(1, -1), nrow = 1))

test.sub2.3 <- isotropy_subsampling(data.3)

test.sub3.3 <- isotropy_subsampling(data.3, lagmat =  rbind(c(1,0), c(0,1), c(1,1), c(1,-1), c(2,1), c(-1, 2), c(1, 2), c(-2, 1)),
                                    A =  rbind(c(1, -1, 0, 0, 0, 0, 0, 0), c(0, 0, 1, -1, 0, 0, 0, 0), c(0, 0, 0, 0, 1, -1, 0, 0), c(0, 0, 0, 0, 0, 0, 1, -1)))

## 4. grid: top right
test.sub1.4 <- isotropy_subsampling(data.4, lagmat =  rbind(c(1,0), c(0,1)), A =  matrix(c(1, -1), nrow = 1))

test.sub2.4 <- isotropy_subsampling(data.4)

test.sub3.4 <- isotropy_subsampling(data.4, lagmat =  rbind(c(1,0), c(0,1), c(1,1), c(1,-1), c(2,1), c(-1, 2), c(1, 2), c(-2, 1)),
                                    A =  rbind(c(1, -1, 0, 0, 0, 0, 0, 0), c(0, 0, 1, -1, 0, 0, 0, 0), c(0, 0, 0, 0, 1, -1, 0, 0), c(0, 0, 0, 0, 0, 0, 1, -1)))

##########
### Graphs
##########

data <- as.data.frame(cbind(reg1.data.long$data.sd, reg1.data.long$x, reg1.data.long$y))

####################
# variogram estimation
varog <- variogram_est(data, hmax = c(4,4,3,3), estimator = c("Matheron", "Genton", "MCD.diff"))


## lags
# S-N
lag.vec.SN <- cbind(rep(0, 4), 1:4)
lags.SN <- apply(lag.vec.SN, 1, function(x) sqrt(t(x) %*% x))

# E-W
lag.vec.EW <- cbind(1:4, rep(0, 4))
lags.EW <- apply(lag.vec.EW, 1, function(x) sqrt(t(x) %*% x))

# SE-NW
lag.vec.SENW <- cbind(1:3, -(1:3))
lags.SENW <- apply(lag.vec.SENW, 1, function(x) sqrt(t(x) %*% x))

# SW-NE
lag.vec.SWNE <- cbind(1:3, 1:3)
lags.SWNE <- apply(lag.vec.SWNE, 1, function(x) sqrt(t(x) %*% x))

lags <- list(lags.SN, lags.EW, lags.SWNE, lags.SENW)   

varo.long <- data.frame(lag = NA, direction = NA, estimator = NA, vario = NA)
dic <- c("S-N", "E-W", "SW-NE", "SE-NW")
est <- c("matheron", "genton", "mcd.diff")

# prepare for ggplot
for(d in 1:4){
  for(e in 1:3){
    ls <- lags[[d]]
    varo.long <- rbind(varo.long, data.frame(lag = ls, direction = dic[d], estimator = est[e], vario = varog[[d]][which(varog[[d]]$estimator == est[e]), "variogram"])) 
  }
}

# colours
cbbPalette2 <- c("#E69F00", "#56B4E9",  "#D55E00","#0072B2")

S_N <- filter(varo.long, direction == "S-N")
S_Ng <- ggplot(S_N, mapping = aes(x = lag, y = vario, col = estimator, shape = estimator)) +
  geom_line(linewidth = 0.25) + geom_point(size= 2) + ylab(expression(2 * gamma(h))) + xlab("||h||") + scale_colour_manual(values=cbbPalette2) +
  theme_minimal(base_size = 10) + scale_shape_manual(values = c(3, 17, 8)) + ylim(c(0,7)) + labs(title = "S-N")


E_W <- filter(varo.long, direction == "E-W")
E_Wg <- ggplot(E_W, mapping = aes(x = lag, y = vario, col = estimator, shape = estimator)) +
  geom_line(linewidth = 0.25) + geom_point(size= 2) + ylab(expression(2 * gamma(h))) + xlab("||h||") + scale_colour_manual(values=cbbPalette2) +
  theme_minimal(base_size = 10)  + scale_shape_manual(values = c(3, 17, 8)) + ylim(c(0,7)) + labs(title = "E-W")

SW_NE <- filter(varo.long, direction == "SW-NE")
SW_NEg <- ggplot(SW_NE, mapping = aes(x = lag, y = vario, col = estimator, shape = estimator)) +
  geom_line(linewidth = 0.25) + geom_point(size= 2) + ylab(expression(2 * gamma(h))) + xlab("||h||") + scale_colour_manual(values=cbbPalette2) +
  theme_minimal(base_size = 10)  + scale_shape_manual(values = c(3, 17, 8)) + ylim(c(0,7)) + labs(title = "SW-NE")

SE_NW <- filter(varo.long, direction == "SE-NW" )
SE_NWg <- ggplot(SE_NW, mapping = aes(x = lag, y = vario, col = estimator, shape = estimator)) +
  geom_line(linewidth = 0.25) + geom_point(size= 2) + ylab(expression(2 * gamma(h))) + xlab("||h||") + scale_colour_manual(values=cbbPalette2) +
  theme_minimal(base_size = 10)  + scale_shape_manual(values = c(3, 17, 8)) + ylim(c(0,7)) + labs(title = "SE-NW")

(S_Ng | E_Wg) / (SW_NEg | SE_NWg)  + plot_layout(guides = "collect") & theme(legend.position = 'bottom')
ggsave("Graphs/Anwendung_MCD.pdf", width = 16.5, units = "cm") # Figure 11


######################
### Amount of outliers
######################

# complete data set
# calculate the pixlewise NDVI for the subregion with only clear data
reg1.clear <- NDVI_sub(-7331535, -999642, "2016-10-20", L8.clear_mask)

# extract the NDVI Index
reg1.clear.data = reg1.clear$tseries 

# rename the y-cordinates
rownames(reg1.clear.data) <- sapply(60:1, function(x) paste0("y", x)) 

# add the coordinates
reg1.clear.data <- cbind(sapply(60:1, function(x) paste0("y", x)), reg1.clear.data)

# rename the x-coordinates
colnames(reg1.clear.data) <- c("y", sapply(1:60, function(x) paste0("x", x)))

# data set with a row for each data point
reg1.clear.data.long <- as_tibble(reg1.clear.data) %>% pivot_longer(cols = starts_with("x"), values_to = "data") %>%  
  mutate(data = as.numeric(data), ID = row_number()) %>% rename(x = name) %>% 
  mutate(x = as.numeric(sub("x", "", x)), y = as.numeric(sub("y", "", y))) 

# standardaize the data with the standard deviation fo the complete data 
reg1.clear.data.long <- reg1.clear.data.long %>% mutate(data.sd = data/sd.all)

# amount of missing values
sum(is.na(reg1.clear.data.long$data))/(60*60)
#[1] 0.1730556

#############
# 30x30 grids

## 1. grid: bottom left
reg1.1.clear <-  NDVI_sub.1(-7331535, -999642, "2016-10-20", L8.clear_mask)

# extract the NDVI Index
reg1.1.clear.data = reg1.1.clear$tseries

sum(is.na(reg1.1.clear.data))/(30*30)
# 0

## 2 grid: top left
reg1.2.clear <-  NDVI_sub.2(-7331535, -999642, "2016-10-20", L8.clear_mask)

# extract the NDVI Index
reg1.2.clear.data = reg1.2.clear$tseries

sum(is.na(reg1.2.clear.data))/(30*30)
# [1] 0.1888889

## 3. grid: bottom right
reg1.3.clear <-  NDVI_sub.3(-7331535, -999642, "2016-10-20", L8.clear_mask)

# extract the NDVI Index
reg1.3.clear.data = reg1.3.clear$tseries

sum(is.na(reg1.3.clear.data))/(30*30)
# [1] 0.1222222

## 4. grid: top right
reg1.4.clear <-  NDVI_sub.4(-7331535, -999642, "2016-10-20", L8.clear_mask)

# extract the NDVI Index
reg1.4.clear.data = reg1.4.clear$tseries

sum(is.na(reg1.4.clear.data))/(30*30)
# [1] 0.3811111


################################################################################
#### Graphics and tables for the paper                                      ####
################################################################################

################
## Load packages
################

library(readxl)
library(ggplot2)
library(tidyr)
library(dplyr)
library(xtable)
library(showtext)
library(geoR)
library(cowplot)

##############
#### Chapter 4
##############

## colors
cbbPalette <- c( "#56B4E9", "#D55E00")


#####################
## 4.1: Gaussian Data
#####################


##########
## Table 3: Percentatges of rejection

load("Results/res_non.RData")

# grid 24x24
res_24 <- res[,,1,,,,,1]
xtable(round(res_24[,1,1,1,1,], 2)*100, digits = 0)

##########
## Table 4: Percentatges of rejection size corrected

load("Results/res.size.RData")

# grid 24x24
res.size_24 <- res[,,1,,,,,1]
xtable(round(res.size_24[,1,1,1,1,], 2)*100, digits = 0)

###################################
## 4.1: Conatmination with outliers
###################################

############
## Figure 2: Isolated outliers

load("Results/res_iso.RData")

# range: 5
res.ia <- res[,,,,2,,,] 
# list("estimator" = c("Matheron", "Genton", "MCD.diff"),
#      "method" = c("subsampling", "block permutation"),
#      "variogram" = c("spherical", "exponential"),
#      "nugget" = c(0, 0.15),
#      "rotation" = c(0, 0.79, 1.18),
#      "scale" = c(0, 1.41, 2),
#      "range"  = c(2,5,8),
#      "Lags" = c(1,2,3),
#      "amount" = c(0.1, 0.2),
#      "dist" = 1:3)
estimator <- c("Matheron", "Genton", "MCD")
method <- c("Subsampling", "Blockpermutation")
rotation <- c("\u03B8 = 0", "\u03B8 = \u03C0/4")
Scale <- c("b = 1", "b = 1.41", "b = 2")
amount <- c(0.1, 0.2)
dist <- c("N(0, 5)", "N(5, 1)")

res.long.IA <- data.frame(Isotropy = NA, rotation= NA, scale = NA, Lag = NA, Amount = NA, Dist = NA, estimator = NA, Method = NA, rejection = NA)

for(i in 1:5){
  if(i == 1){rot <- 1; s <- 1; iso <- "\u03B8 = 0, b = 1"}
  if(i == 2){rot <- 1; s <- 2; iso <- "\u03B8 = 0, b = 1.41"}
  if(i == 3){rot <- 1; s <- 3; iso <- "\u03B8 = 0, b = 2"}
  if(i == 4){rot <- 2; s <- 2; iso <- "\u03B8 = \u03C0/4, b = 1.41"}
  if(i == 5){rot <- 2; s <- 3; iso <- "\u03B8 = \u03C0/4, b = 2"}
  for(a in 1:2){
    for(d in 1:2){
      for(l in 1:3){
        for(es in 1:3){
          for(m in 1:2){
            res.long.IA <- rbind(res.long.IA, data.frame(Isotropy = iso, rotation  = rotation[rot], scale = Scale[s], Amount = amount[a], Dist = dist[d], Lag = l, estimator = estimator[es], Method = method[m],  rejection = res.ia[es,m,rot,s, l, a, d]))
          }
        }
      }
    }
  }
}
res.long.IA$fill_col <- ifelse(res.long.IA$rotation == "\u03B8 = 0", res.long.IA$Method, "open")

# Legende 
legend_plot <- ggplot() + xlim(0, 1) + ylim(0, 1) + annotate("text", x = 0, y = 1, label = "Method", fontface = "bold", hjust = 0) +
  annotate("point", x = 0.1, y = 0.9, shape = 16, color = "#D55E00", size = 2) +
  annotate("text", x = 0.2, y = 0.9, label = "Block permutation", hjust = 0) +
  annotate("point", x = 0.1, y = 0.8, shape = 16, color = "#56B4E9", size = 2) +
  annotate("text", x = 0.2, y = 0.8, label = "Subsampling", hjust = 0) +
  annotate("text", x = 0, y = 0.7, label = "Rotation", fontface = "bold", hjust = 0) +
  annotate("point", x = 0.1, y = 0.6, shape = 16, size = 2) + 
  annotate("text", x = 0.2, y = 0.6, label = expression(theta == 0), hjust = 0) + 
  annotate("point", x = 0.1, y = 0.5, shape = 21, size = 2) +
  annotate("text", x = 0.2, y = 0.5, label = expression(theta == pi/4), hjust = 0) + 
  annotate("text", x = 0, y = 0.4, label = "Scale", fontface = "bold", hjust = 0) +
  annotate("point", x = 0.1, y = 0.3, shape = 21, size = 2) + 
  annotate("text", x = 0.2, y = 0.3, label = expression(b == 1), hjust = 0) + 
  annotate("point", x = 0.1, y = 0.2, shape = 22, size = 2) +
  annotate("text", x = 0.2, y = 0.2, label = expression(b == sqrt(2)), hjust = 0) + 
  annotate("point", x = 0.1, y = 0.1, shape = 24, size = 2) +
  annotate("text", x = 0.2, y = 0.1, label = expression(b == 2), hjust = 0) + 
  theme_void()

leer <- ggplot() + xlim(0, 1) + ylim(0, 1) +   theme_void()

# Amount 0.1, Dist 2
L1.D2.A1 <-  filter(res.long.IA, Amount == 0.1, Dist == "N(5, 1)", Lag == 1)
PL1.D2.A1 <- ggplot(data = L1.D2.A1, mapping = aes(x = estimator, y = rejection, col = Method, shape = scale, fill = fill_col)) + geom_point(size = 2, alpha = 1, stroke = 1.2) + 
  ylab("Proportion of Rejections") + labs(title = expression(epsilon == 0.1 ~ ", N(5, 1)")) + theme_minimal(base_size = 10) + geom_hline(yintercept = 0.05) + ylim(0, 1) + 
  scale_shape_manual(values = c(21, 22, 24))  + scale_fill_manual(name = "rotation", values = c("Subsampling" = "#56B4E9", "Blockpermutation" = "#D55E00", "open" = "white")) +
  scale_color_manual(values = c("Subsampling" = "#56B4E9", "Blockpermutation" = "#D55E00")) + theme(legend.position = "none")

L2.D2.A1 <-  filter(res.long.IA, Amount == 0.1, Dist == "N(5, 1)", Lag == 2)
PL2.D2.A1 <- ggplot(data = L2.D2.A1, mapping = aes(x = estimator, y = rejection, col = Method, shape = scale, fill = fill_col)) + geom_point(size = 2, alpha = 1, stroke = 1.2) + 
  ylab("Proportion of Rejections")  + theme_minimal(base_size = 10) + geom_hline(yintercept = 0.05) + ylim(0, 1) + 
  scale_shape_manual(values = c(21, 22, 24))  + scale_fill_manual(name = "rotation", values = c("Subsampling" = "#56B4E9", "Blockpermutation" = "#D55E00", "open" = "white")) +
  scale_color_manual(values = c("Subsampling" = "#56B4E9", "Blockpermutation" = "#D55E00")) + theme(legend.position = "none")

# Amount 0.2, Dist 2
L1.D2.A2 <-  filter(res.long.IA, Amount == 0.2, Dist == "N(5, 1)", Lag == 1)
PL1.D2.A2 <- ggplot(data = L1.D2.A2, mapping = aes(x = estimator, y = rejection, col = Method, shape = scale, fill = fill_col)) + geom_point(size = 2, alpha = 1, stroke = 1.2) + 
  ylab("Proportion of Rejections") + labs(title = expression(epsilon == 0.2 ~ ", N(5, 1)")) + theme_minimal(base_size = 10) + geom_hline(yintercept = 0.05) + ylim(0, 1) + 
  scale_shape_manual(values = c(21, 22, 24))  + scale_fill_manual(name = "rotation", values = c("Subsampling" = "#56B4E9", "Blockpermutation" = "#D55E00", "open" = "white")) +
  scale_color_manual(values = c("Subsampling" = "#56B4E9", "Blockpermutation" = "#D55E00")) + theme(legend.position = "none")


L2.D2.A2 <-  filter(res.long.IA, Amount == 0.2, Dist == "N(5, 1)", Lag == 2)
PL2.D2.A2 <- ggplot(data = L2.D2.A2, mapping = aes(x = estimator, y = rejection, col = Method, shape = scale, fill = fill_col)) + geom_point(size = 2, alpha = 1, stroke = 1.2) + 
  ylab("Proportion of Rejections")  + theme_minimal(base_size = 10) + geom_hline(yintercept = 0.05) + ylim(0, 1) + 
  scale_shape_manual(values = c(21, 22, 24))  + scale_fill_manual(name = "rotation", values = c("Subsampling" = "#56B4E9", "Blockpermutation" = "#D55E00", "open" = "white")) +
  scale_color_manual(values = c("Subsampling" = "#56B4E9", "Blockpermutation" = "#D55E00")) + theme(legend.position = "none")

# Amount 0.1, Dist 1
L1.D1.A1 <-  filter(res.long.IA, Amount == 0.1, Dist == "N(0, 5)", Lag == 1)
PL1.D1.A1 <- ggplot(data = L1.D1.A1, mapping = aes(x = estimator, y = rejection, col = Method, shape = scale, fill = fill_col)) + geom_point(size = 2, alpha = 1, stroke = 1.2) + 
  ylab("Proportion of Rejections") + labs(title = expression(epsilon == 0.1 ~ ", N(0, 5)")) + theme_minimal(base_size = 10) + geom_hline(yintercept = 0.05) + ylim(0, 1) + 
  scale_shape_manual(values = c(21, 22, 24))  + scale_fill_manual(name = "rotation", values = c("Subsampling" = "#56B4E9", "Blockpermutation" = "#D55E00", "open" = "white")) +
  scale_color_manual(values = c("Subsampling" = "#56B4E9", "Blockpermutation" = "#D55E00")) + theme(legend.position = "none")


L2.D1.A1 <-  filter(res.long.IA, Amount == 0.1, Dist == "N(0, 5)", Lag == 2)
PL2.D1.A1 <- ggplot(data = L2.D1.A1, mapping = aes(x = estimator, y = rejection, col = Method, shape = scale, fill = fill_col)) + geom_point(size = 2, alpha = 1, stroke = 1.2) + 
  ylab("Proportion of Rejections") + theme_minimal(base_size = 10) + geom_hline(yintercept = 0.05) + ylim(0, 1) + 
  scale_shape_manual(values = c(21, 22, 24))  + scale_fill_manual(name = "rotation", values = c("Subsampling" = "#56B4E9", "Blockpermutation" = "#D55E00", "open" = "white")) +
  scale_color_manual(values = c("Subsampling" = "#56B4E9", "Blockpermutation" = "#D55E00")) + theme(legend.position = "none")

plot_grid(PL1.D2.A1, PL2.D2.A1, leer, PL1.D2.A2, PL2.D2.A2, legend_plot, PL1.D1.A1, PL2.D1.A1, leer, ncol = 3, rel_widths = c(2,2,1))
ggsave("Graphs/IA.pdf", width = 22, height = 22, units = "cm")


############
## Figure 3: Random outlier block

load("Results/res_block.RData")

# range: 5, random block
res_BA_rand <-  res[,,,,2,,,,1] 

# list("estimator" = c("Matheron", "Genton", "MCD.diff"),
#      "method" = c("subsampling", "blockpermutation"),
#      "nugget" = c(0, 0.15),
#      "rotation" = c(0, 0.79),
#      "scale" = c(0, 1.41),
#      "range"  = c(2,5,8),
#      "Lags" = c(1,2,3),
#      "amount" = c(0.1, 0.2),
#      "dist" = 1:2)
estimator <- c("Matheron", "Genton", "MCD")
method <- c("Subsampling", "Blockpermutation")
rotation <- c("\u03B8 = 0", "\u03B8 = \u03C0/4")
Scale <- c("b = 1", "b = 1.41", "b = 2")
amount <- c(0.1, 0.2)
dist <- c("N(0, 5)", "N(5, 1)")

res.long.BA.rand <- data.frame(Isotropy = NA, rotation= NA, scale = NA, Lag = NA, Amount = NA, Dist = NA, estimator = NA, Method = NA, rejection = NA)

for(i in 1:5){
  if(i == 1){rot <- 1; s <- 1; iso <- "\u03B8 = 0, b = 1"}
  if(i == 2){rot <- 1; s <- 2; iso <- "\u03B8 = 0, b = 1.41"}
  if(i == 3){rot <- 1; s <- 3; iso <- "\u03B8 = 0, b = 2"}
  if(i == 4){rot <- 2; s <- 2; iso <- "\u03B8 = \u03C0/4, b = 1.41"}
  if(i == 5){rot <- 2; s <- 3; iso <- "\u03B8 = \u03C0/4, b = 2"}  
  for(a in 1:2){
    for(d in 1:2){
      for(l in 1:3){
        for(es in 1:3){
          for(m in 1:2){
            res.long.BA.rand <- rbind(res.long.BA.rand, data.frame(Isotropy = iso, rotation  = rotation[rot], scale = Scale[s], Amount = amount[a], Dist = dist[d], Lag = l, estimator = estimator[es], Method = method[m],  rejection = res_BA_rand[es,m,rot,s, l, a, d]))
          }
        }
      }
    }
  }
}

res.long.BA.rand$fill_col <- ifelse(res.long.BA.rand$rotation == "\u03B8 = 0", res.long.BA.rand$Method, "open")

# Amount 0.1, Dist 2
L1.D2.A1 <-  filter(res.long.BA.rand, Amount == 0.1, Dist == "N(5, 1)", Lag == 1)
PL1.D2.A1 <- ggplot(data = L1.D2.A1, mapping = aes(x = estimator, y = rejection, col = Method, shape = scale, fill = fill_col)) + geom_point(size = 2, alpha = 1, stroke = 1.2) + 
  ylab("Proportion of Rejections") + labs(title = expression(epsilon == 0.1 ~ ", N(5, 1)")) + theme_minimal(base_size = 10) + geom_hline(yintercept = 0.05) + ylim(0, 1) + 
  scale_shape_manual(values = c(21, 22, 24))  + scale_fill_manual(name = "rotation", values = c("Subsampling" = "#56B4E9", "Blockpermutation" = "#D55E00", "open" = "white")) +
  scale_color_manual(values = c("Subsampling" = "#56B4E9", "Blockpermutation" = "#D55E00")) + theme(legend.position = "none")

L2.D2.A1 <-  filter(res.long.BA.rand, Amount == 0.1, Dist == "N(5, 1)", Lag == 2)
PL2.D2.A1 <- ggplot(data = L2.D2.A1, mapping = aes(x = estimator, y = rejection, col = Method, shape = scale, fill = fill_col)) + geom_point(size = 2, alpha = 1, stroke = 1.2) + 
  ylab("Proportion of Rejections")  + theme_minimal(base_size = 10) + geom_hline(yintercept = 0.05) + ylim(0, 1) + 
  scale_shape_manual(values = c(21, 22, 24))  + scale_fill_manual(name = "rotation", values = c("Subsampling" = "#56B4E9", "Blockpermutation" = "#D55E00", "open" = "white")) +
  scale_color_manual(values = c("Subsampling" = "#56B4E9", "Blockpermutation" = "#D55E00")) + theme(legend.position = "none")

# Amount 0.2, Dist 2
L1.D2.A2 <-  filter(res.long.BA.rand, Amount == 0.2, Dist == "N(5, 1)", Lag == 1)
PL1.D2.A2 <- ggplot(data = L1.D2.A2, mapping = aes(x = estimator, y = rejection, col = Method, shape = scale, fill = fill_col)) + geom_point(size = 2, alpha = 1, stroke = 1.2) + 
  ylab("Proportion of Rejections") + labs(title = expression(epsilon == 0.2 ~ ", N(5, 1)")) + theme_minimal(base_size = 10) + geom_hline(yintercept = 0.05) + ylim(0, 1) + 
  scale_shape_manual(values = c(21, 22, 24))  + scale_fill_manual(name = "rotation", values = c("Subsampling" = "#56B4E9", "Blockpermutation" = "#D55E00", "open" = "white")) +
  scale_color_manual(values = c("Subsampling" = "#56B4E9", "Blockpermutation" = "#D55E00")) + theme(legend.position = "none")


L2.D2.A2 <-  filter(res.long.BA.rand, Amount == 0.2, Dist == "N(5, 1)", Lag == 2)
PL2.D2.A2 <- ggplot(data = L2.D2.A2, mapping = aes(x = estimator, y = rejection, col = Method, shape = scale, fill = fill_col)) + geom_point(size = 2, alpha = 1, stroke = 1.2) + 
  ylab("Proportion of Rejections")  + theme_minimal(base_size = 10) + geom_hline(yintercept = 0.05) + ylim(0, 1) + 
  scale_shape_manual(values = c(21, 22, 24))  + scale_fill_manual(name = "rotation", values = c("Subsampling" = "#56B4E9", "Blockpermutation" = "#D55E00", "open" = "white")) +
  scale_color_manual(values = c("Subsampling" = "#56B4E9", "Blockpermutation" = "#D55E00")) + theme(legend.position = "none")

# Amount 0.1, Dist 1
L1.D1.A1 <-  filter(res.long.BA.rand, Amount == 0.1, Dist == "N(0, 5)", Lag == 1)
PL1.D1.A1 <- ggplot(data = L1.D1.A1, mapping = aes(x = estimator, y = rejection, col = Method, shape = scale, fill = fill_col)) + geom_point(size = 2, alpha = 1, stroke = 1.2) + 
  ylab("Proportion of Rejections") + labs(title = expression(epsilon == 0.1 ~ ", N(0, 5)")) + theme_minimal(base_size = 10) + geom_hline(yintercept = 0.05) + ylim(0, 1) + 
  scale_shape_manual(values = c(21, 22, 24))  + scale_fill_manual(name = "rotation", values = c("Subsampling" = "#56B4E9", "Blockpermutation" = "#D55E00", "open" = "white")) +
  scale_color_manual(values = c("Subsampling" = "#56B4E9", "Blockpermutation" = "#D55E00")) + theme(legend.position = "none")


L2.D1.A1 <-  filter(res.long.BA.rand, Amount == 0.1, Dist == "N(0, 5)", Lag == 2)
PL2.D1.A1 <- ggplot(data = L2.D1.A1, mapping = aes(x = estimator, y = rejection, col = Method, shape = scale, fill = fill_col)) + geom_point(size = 2, alpha = 1, stroke = 1.2) + 
  ylab("Proportion of Rejections") + theme_minimal(base_size = 10) + geom_hline(yintercept = 0.05) + ylim(0, 1) + 
  scale_shape_manual(values = c(21, 22, 24))  + scale_fill_manual(name = "rotation", values = c("Subsampling" = "#56B4E9", "Blockpermutation" = "#D55E00", "open" = "white")) +
  scale_color_manual(values = c("Subsampling" = "#56B4E9", "Blockpermutation" = "#D55E00")) + theme(legend.position = "none")

plot_grid(PL1.D2.A1, PL2.D2.A1, leer, PL1.D2.A2, PL2.D2.A2, legend_plot, PL1.D1.A1, PL2.D1.A1, leer, ncol = 3, rel_widths = c(2,2,1))
ggsave("Graphs/BA.pdf", width = 22, height = 22, units = "cm")


############
## Figure 4: Rectangle outlier block

load("Results/res_block.RData")

# range: 5, random block, scale 1 & 1.41, amount 0.2, dist 2
res_BA_rec <-  res[,,,1:2,2,,2,2,2] 


estimator <- c("Matheron", "Genton", "MCD")
method <- c("Subsampling", "Blockpermutation")
rotation <- c("\u03B8 = 0", "\u03B8 = \u03C0/4")
Scale <- c("b = 1", "b = 1.41")
amount <- c(0.1, 0.2)
dist <- c("N(0, 5)", "N(5, 1)")

res.long.BA.rec <- data.frame(Isotropy = NA, rotation= NA, scale = NA, Lag = NA, Amount = NA, Dist = NA, estimator = NA, Method = NA, rejection = NA)

for(i in 1:3){
  if(i == 1){rot <- 1; s <- 1; iso <- "\u03B8 = 0, b = 1"}
  if(i == 2){rot <- 1; s <- 2; iso <- "\u03B8 = 0, b = 1.41"}
  if(i == 3){rot <- 2; s <- 2; iso <- "\u03B8 = \u03C0/4, b = 1.41"}
  for(l in 1:3){
    for(es in 1:3){
      for(m in 1:2){
        res.long.BA.rec <- rbind(res.long.BA.rec, data.frame(Isotropy = iso, rotation  = rotation[rot], scale = Scale[s], Amount = amount[2], Dist = dist[2], Lag = l, estimator = estimator[es], Method = method[m],  rejection = res_BA_rec[es,m,rot,s, l]))
      }
    }
  }
}

res.long.BA.rec$fill_col <- ifelse(res.long.BA.rec$rotation == "\u03B8 = 0", res.long.BA.rec$Method, "open")

# Legende 
legend_plot <- ggplot() + xlim(0, 1) + ylim(0, 1) + annotate("text", x = 0, y = 1, label = "Method", fontface = "bold", hjust = 0) +
  annotate("point", x = 0.1, y = 0.9, shape = 16, color = "#D55E00", size = 2) +
  annotate("text", x = 0.2, y = 0.9, label = "Block permutation", hjust = 0) +
  annotate("point", x = 0.1, y = 0.8, shape = 16, color = "#56B4E9", size = 2) +
  annotate("text", x = 0.2, y = 0.8, label = "Subsampling", hjust = 0) +
  annotate("text", x = 0, y = 0.7, label = "Rotation", fontface = "bold", hjust = 0) +
  annotate("point", x = 0.1, y = 0.6, shape = 16, size = 2) + 
  annotate("text", x = 0.2, y = 0.6, label = expression(theta == 0), hjust = 0) + 
  annotate("point", x = 0.1, y = 0.5, shape = 21, size = 2) +
  annotate("text", x = 0.2, y = 0.5, label = expression(theta == pi/4), hjust = 0) + 
  annotate("text", x = 0, y = 0.4, label = "Scale", fontface = "bold", hjust = 0) +
  annotate("point", x = 0.1, y = 0.3, shape = 21, size = 2) + 
  annotate("text", x = 0.2, y = 0.3, label = expression(b == 1), hjust = 0) + 
  annotate("point", x = 0.1, y = 0.2, shape = 22, size = 2) +
  annotate("text", x = 0.2, y = 0.2, label = expression(b == sqrt(2)), hjust = 0) + 
  theme_void()


# Amount 0.2, Dist 2
L1.D2.A2 <-  filter(res.long.BA.rec, Amount == 0.2, Dist == "N(5, 1)", Lag == 1)
PL1.D2.A2 <- ggplot(data = L1.D2.A2, mapping = aes(x = estimator, y = rejection, col = Method, shape = scale, fill = fill_col)) + geom_point(size = 2, alpha = 1, stroke = 1.2) + 
  ylab("Proportion of Rejections") + labs(title = expression(epsilon == 0.2 ~ ", N(5, 1)")) + theme_minimal(base_size = 10) + geom_hline(yintercept = 0.05) + ylim(0, 1) + 
  scale_shape_manual(values = c(21, 22))  + scale_fill_manual(name = "rotation", values = c("Subsampling" = "#56B4E9", "Blockpermutation" = "#D55E00", "open" = "white")) +
  scale_color_manual(values = c("Subsampling" = "#56B4E9", "Blockpermutation" = "#D55E00")) + theme(legend.position = "none")


L2.D2.A2 <-  filter(res.long.BA.rec, Amount == 0.2, Dist == "N(5, 1)", Lag == 2)
PL2.D2.A2 <- ggplot(data = L2.D2.A2, mapping = aes(x = estimator, y = rejection, col = Method, shape = scale, fill = fill_col)) + geom_point(size = 2, alpha = 1, stroke = 1.2) + 
  ylab("Proportion of Rejections")  + theme_minimal(base_size = 10) + geom_hline(yintercept = 0.05) + ylim(0, 1) + 
  scale_shape_manual(values = c(21, 22))  + scale_fill_manual(name = "rotation", values = c("Subsampling" = "#56B4E9", "Blockpermutation" = "#D55E00", "open" = "white")) +
  scale_color_manual(values = c("Subsampling" = "#56B4E9", "Blockpermutation" = "#D55E00")) + theme(legend.position = "none")

plot_grid(PL1.D2.A2, PL2.D2.A2, legend_plot, ncol = 3, rel_widths = c(2,2,1))
ggsave("Graphs/BA_rec.pdf", width = 22, height = 8, units = "cm")



##############
#### Appendix
##############

############
## figure 1: anisotropie & true variogram

# spherical variogram
sp.var <- function(h, a = 5, c = 0.4, n = 0.1){
  if(h < a){res <- n + c*(((3*h)/(2*a)) - (1/2) * (h/a)^3)}
  if(h >= a){res <- n + c}
  if(h == 0){res <- 0}
  return(res)
}

x.lag <- seq(0,6, 0.1)
y.lag <- sapply(1:length(x.lag), function(x) sp.var(x.lag[x], a = 5, c = 1, n = 0))

true.long <- data.frame("h" = x.lag, 
                        "vario" = y.lag)

plot.sph <- ggplot(true.long, aes(x = h, y = vario)) + geom_line(linewidth = 0.25) +
  xlab("||h||") + ylab(expression(gamma(h))) + 
  theme_minimal(base_size = 10) + xlim(c(0,6)) + ylim(c(0,1))

# visualisation anisotropie
theta <- seq(0, 2*pi, length.out = 1000)

#1. theta = 0, b = 1
vecs <- cbind(1*cos(theta), 1*sin(theta))
vecs.trans <- apply(vecs, 1, function(x){coords.aniso(matrix(x, ncol = 2), c(0, 1), reverse = TRUE)})
vecs.trans <- data.frame("x" = vecs.trans[1,], "y" = vecs.trans[2,])

Iso.1 <- ggplot(vecs.trans, aes(x, y)) + geom_path(color = "black") +
  coord_fixed() + theme_minimal(base_size = 10) + labs(title = expression(theta == 0 ~ "," ~ b == 1)) +
  ylab("y lag") + xlab("x lag") + ylim(c(-2,2)) + xlim(c(-2,2)) + geom_hline(yintercept = 0,  linewidth = 1) +  
  geom_vline(xintercept = 0, linewidth = 1)  

#2. theta = 0, b = sqrt(2)
vecs <- cbind(1*cos(theta), 1*sin(theta))
vecs.trans <- apply(vecs, 1, function(x){coords.aniso(matrix(x, ncol = 2), c(0, sqrt(2)), reverse = TRUE)})
vecs.trans <- data.frame("x" = vecs.trans[1,], "y" = vecs.trans[2,])

Iso.2 <- ggplot(vecs.trans, aes(x, y)) + geom_path(color = "black") +
  coord_fixed() + theme_minimal(base_size = 10) + labs(title = expression(theta == 0 ~ "," ~ b == sqrt(2))) +
  ylab("y lag") + xlab("x lag") + ylim(c(-2,2)) + xlim(c(-2,2)) + geom_hline(yintercept = 0,  linewidth = 1) +  
  geom_vline(xintercept = 0, linewidth = 1)  

#3. theta = 0, b = 2
vecs <- cbind(1*cos(theta), 1*sin(theta))
vecs.trans <- apply(vecs, 1, function(x){coords.aniso(matrix(x, ncol = 2), c(0, 2), reverse = TRUE)})
vecs.trans <- data.frame("x" = vecs.trans[1,], "y" = vecs.trans[2,])

Iso.3 <- ggplot(vecs.trans, aes(x, y)) + geom_path(color = "black") +
  coord_fixed() + theme_minimal(base_size = 10) + labs(title = expression(theta == 0 ~ "," ~ b == 2)) +
  ylab("y lag") + xlab("x lag")  + ylim(c(-2,2)) + xlim(c(-2,2)) + geom_hline(yintercept = 0,  linewidth = 1) +  
  geom_vline(xintercept = 0, linewidth = 1)  

#4. theta = pi/4, b = sqrt(2)
vecs <- cbind(1*cos(theta), 1*sin(theta))
vecs.trans <- apply(vecs, 1, function(x){coords.aniso(matrix(x, ncol = 2), aniso.pars = c(pi/4, sqrt(2)), reverse = TRUE)})
vecs.trans <- data.frame("x" = vecs.trans[1,], "y" = vecs.trans[2,]) 

Iso.4 <- ggplot(vecs.trans, aes(x, y)) + geom_path(color = "black") +
  coord_fixed() + theme_minimal(base_size = 10) + labs(title = expression(theta == pi/4 ~ "," ~ b == sqrt(2))) +
  ylab("y lag") + xlab("x lag") + ylim(c(-2,2)) + xlim(c(-2,2)) + geom_hline(yintercept = 0,  linewidth = 1) +  
  geom_vline(xintercept = 0, linewidth = 1)  


#5. theta = pi*0.375, b = sqrt(2)
vecs <- cbind(1*cos(theta), 1*sin(theta))
vecs.trans <- apply(vecs, 1, function(x){coords.aniso(matrix(x, ncol = 2), c(pi*0.375, sqrt(2)), reverse = TRUE)})
vecs.trans <- data.frame("x" = vecs.trans[1,], "y" = vecs.trans[2,])

Iso.5 <-  ggplot(vecs.trans, aes(x, y)) + geom_path(color = "black") +
  coord_fixed() + theme_minimal(base_size = 10) + labs(title = expression(theta == 0.375*pi ~ "," ~ b == sqrt(2))) +
  ylab("y lag") + xlab("x lag") + ylim(c(-2,2)) + xlim(c(-2,2)) + geom_hline(yintercept = 0,  linewidth = 1) +  
  geom_vline(xintercept = 0, linewidth = 1)  

library(patchwork)
(plot.sph | Iso.1 | Iso.2)/(Iso.3 | Iso.4 | Iso.5) + plot_layout(guides = "collect") & theme(legend.position = 'right')
ggsave("Graphs/Anisotropy.pdf", width = 22, height = 11, units = "cm")


############
## figure 2: isolated outliers, labmda 3

legend_plot <- ggplot() + xlim(0, 1) + ylim(0, 1) + annotate("text", x = 0, y = 1, label = "Method", fontface = "bold", hjust = 0) +
  annotate("point", x = 0.1, y = 0.9, shape = 16, color = "#D55E00", size = 2) +
  annotate("text", x = 0.2, y = 0.9, label = "Block permutation", hjust = 0) +
  annotate("point", x = 0.1, y = 0.8, shape = 16, color = "#56B4E9", size = 2) +
  annotate("text", x = 0.2, y = 0.8, label = "Subsampling", hjust = 0) +
  annotate("text", x = 0, y = 0.7, label = "Rotation", fontface = "bold", hjust = 0) +
  annotate("point", x = 0.1, y = 0.6, shape = 16, size = 2) + 
  annotate("text", x = 0.2, y = 0.6, label = expression(theta == 0), hjust = 0) + 
  annotate("point", x = 0.1, y = 0.5, shape = 21, size = 2) +
  annotate("text", x = 0.2, y = 0.5, label = expression(theta == pi/4), hjust = 0) + 
  annotate("text", x = 0, y = 0.4, label = "Scale", fontface = "bold", hjust = 0) +
  annotate("point", x = 0.1, y = 0.3, shape = 21, size = 2) + 
  annotate("text", x = 0.2, y = 0.3, label = expression(b == 1), hjust = 0) + 
  annotate("point", x = 0.1, y = 0.2, shape = 22, size = 2) +
  annotate("text", x = 0.2, y = 0.2, label = expression(b == sqrt(2)), hjust = 0) + 
  annotate("point", x = 0.1, y = 0.1, shape = 24, size = 2) +
  annotate("text", x = 0.2, y = 0.1, label = expression(b == 2), hjust = 0) + 
  theme_void()

leer <- ggplot() + xlim(0, 1) + ylim(0, 1) +   theme_void()

L3.D2.A1 <-  filter(res.long.IA, Amount == 0.1, Dist == "N(5, 1)", Lag == 3)
PL3.D2.A1 <- ggplot(data = L3.D2.A1, mapping = aes(x = estimator, y = rejection, col = Method, shape = scale, fill = fill_col)) + geom_point(size = 2, alpha = 1, stroke = 1.2) + 
  ylab("Proportion of Rejections") + labs(title = expression(epsilon == 0.1 ~ ", N(5, 1)")) + theme_minimal(base_size = 10) + geom_hline(yintercept = 0.05) + ylim(0, 1) + 
  scale_shape_manual(values = c(21, 22, 24))  + scale_fill_manual(name = "rotation", values = c("Subsampling" = "#56B4E9", "Blockpermutation" = "#D55E00", "open" = "white")) +
  scale_color_manual(values = c("Subsampling" = "#56B4E9", "Blockpermutation" = "#D55E00")) + theme(legend.position = "none")


L3.D2.A2 <-  filter(res.long.IA, Amount == 0.2, Dist == "N(5, 1)", Lag == 3)
PL3.D2.A2 <- ggplot(data = L3.D2.A2, mapping = aes(x = estimator, y = rejection, col = Method, shape = scale, fill = fill_col)) + geom_point(size = 2, alpha = 1, stroke = 1.2) + 
  ylab("Proportion of Rejections") + labs(title = expression(epsilon == 0.2 ~ ", N(5, 1)")) + theme_minimal(base_size = 10) + geom_hline(yintercept = 0.05) + ylim(0, 1) + 
  scale_shape_manual(values = c(21, 22, 24))  + scale_fill_manual(name = "rotation", values = c("Subsampling" = "#56B4E9", "Blockpermutation" = "#D55E00", "open" = "white")) +
  scale_color_manual(values = c("Subsampling" = "#56B4E9", "Blockpermutation" = "#D55E00")) + theme(legend.position = "none")


L3.D1.A1 <-  filter(res.long.IA, Amount == 0.1, Dist == "N(0, 5)", Lag == 3)
PL3.D1.A1 <-  ggplot(data = L3.D1.A1, mapping = aes(x = estimator, y = rejection, col = Method, shape = scale, fill = fill_col)) + geom_point(size = 2, alpha = 1, stroke = 1.2) + 
  ylab("Proportion of Rejections") + labs(title = expression(epsilon == 0.1 ~ ", N(0, 5)")) + theme_minimal(base_size = 10) + geom_hline(yintercept = 0.05) + ylim(0, 1) + 
  scale_shape_manual(values = c(21, 22, 24))  + scale_fill_manual(name = "rotation", values = c("Subsampling" = "#56B4E9", "Blockpermutation" = "#D55E00", "open" = "white")) +
  scale_color_manual(values = c("Subsampling" = "#56B4E9", "Blockpermutation" = "#D55E00")) + theme(legend.position = "none")


plot_grid(PL3.D2.A1, leer, PL3.D2.A2, legend_plot, PL3.D1.A1, leer, ncol = 2, rel_widths = c(2,1))
ggsave("Graphs/IA2.pdf", width = 22, height = 22, units = "cm")


############
## figure 3: random block outliers, labmda 3

L3.D2.A1 <-  filter(res.long.BA.rand, Amount == 0.1, Dist == "N(5, 1)", Lag == 3)
PL3.D2.A1 <- ggplot(data = L3.D2.A1, mapping = aes(x = estimator, y = rejection, col = Method, shape = scale, fill = fill_col)) + geom_point(size = 2, alpha = 1, stroke = 1.2) + 
  ylab("Proportion of Rejections") + labs(title = expression(epsilon == 0.1 ~ ", N(5, 1)")) + theme_minimal(base_size = 10) + geom_hline(yintercept = 0.05) + ylim(0, 1) + 
  scale_shape_manual(values = c(21, 22, 24))  + scale_fill_manual(name = "rotation", values = c("Subsampling" = "#56B4E9", "Blockpermutation" = "#D55E00", "open" = "white")) +
  scale_color_manual(values = c("Subsampling" = "#56B4E9", "Blockpermutation" = "#D55E00")) + theme(legend.position = "none")

L3.D2.A2 <-  filter(res.long.BA.rand, Amount == 0.2, Dist == "N(5, 1)", Lag == 3)
PL3.D2.A2 <- ggplot(data = L3.D2.A2, mapping = aes(x = estimator, y = rejection, col = Method, shape = scale, fill = fill_col)) + geom_point(size = 2, alpha = 1, stroke = 1.2) + 
  ylab("Proportion of Rejections") + labs(title = expression(epsilon == 0.2 ~ ", N(5, 1)")) + theme_minimal(base_size = 10) + geom_hline(yintercept = 0.05) + ylim(0, 1) + 
  scale_shape_manual(values = c(21, 22, 24))  + scale_fill_manual(name = "rotation", values = c("Subsampling" = "#56B4E9", "Blockpermutation" = "#D55E00", "open" = "white")) +
  scale_color_manual(values = c("Subsampling" = "#56B4E9", "Blockpermutation" = "#D55E00")) + theme(legend.position = "none")

L3.D1.A1 <-  filter(res.long.BA.rand, Amount == 0.1, Dist == "N(0, 5)", Lag == 3)
PL3.D1.A1 <-  ggplot(data = L3.D1.A1, mapping = aes(x = estimator, y = rejection, col = Method, shape = scale, fill = fill_col)) + geom_point(size = 2, alpha = 1, stroke = 1.2) + 
  ylab("Proportion of Rejections") + labs(title = expression(epsilon == 0.1 ~ ", N(0, 5)")) + theme_minimal(base_size = 10) + geom_hline(yintercept = 0.05) + ylim(0, 1) + 
  scale_shape_manual(values = c(21, 22, 24))  + scale_fill_manual(name = "rotation", values = c("Subsampling" = "#56B4E9", "Blockpermutation" = "#D55E00", "open" = "white")) +
  scale_color_manual(values = c("Subsampling" = "#56B4E9", "Blockpermutation" = "#D55E00")) + theme(legend.position = "none")

plot_grid(PL3.D2.A1, leer, PL3.D2.A2, legend_plot, PL3.D1.A1, leer, ncol = 2, rel_widths = c(2,1))
ggsave("Graphs/BA2.pdf", width = 22, height = 22, units = "cm")


############
## figure 4: quadratic block outliers

# Legende 
legend_plot <- ggplot() + xlim(0, 1) + ylim(0, 1) + annotate("text", x = 0, y = 1, label = "Method", fontface = "bold", hjust = 0) +
  annotate("point", x = 0.1, y = 0.9, shape = 16, color = "#D55E00", size = 2) +
  annotate("text", x = 0.2, y = 0.9, label = "Block permutation", hjust = 0) +
  annotate("point", x = 0.1, y = 0.8, shape = 16, color = "#56B4E9", size = 2) +
  annotate("text", x = 0.2, y = 0.8, label = "Subsampling", hjust = 0) +
  annotate("text", x = 0, y = 0.7, label = "Rotation", fontface = "bold", hjust = 0) +
  annotate("point", x = 0.1, y = 0.6, shape = 16, size = 2) + 
  annotate("text", x = 0.2, y = 0.6, label = expression(theta == 0), hjust = 0) + 
  annotate("point", x = 0.1, y = 0.5, shape = 21, size = 2) +
  annotate("text", x = 0.2, y = 0.5, label = expression(theta == pi/4), hjust = 0) + 
  annotate("text", x = 0, y = 0.4, label = "Scale", fontface = "bold", hjust = 0) +
  annotate("point", x = 0.1, y = 0.3, shape = 21, size = 2) + 
  annotate("text", x = 0.2, y = 0.3, label = expression(b == 1), hjust = 0) + 
  annotate("point", x = 0.1, y = 0.2, shape = 22, size = 2) +
  annotate("text", x = 0.2, y = 0.2, label = expression(b == sqrt(2)), hjust = 0) + 
  theme_void()

load("Results/res_block.RData")

# scale 1.41, amount 2 , dist 2, range 5
res_BA_sq <- res[,,,1:2,2,,2,2,3]

estimator <- c("Matheron", "Genton", "MCD")
method <- c("Subsampling", "Block permutation")
rotation <- c("\u03B8 = 0", "\u03B8 = \u03C0/4")
Scale <- c("b = 1", "b = 1.41", "b = 2")
amount <- c(0.1, 0.2)
dist <- c("N(0, 5)", "N(5, 1)")

res.long.BA.sq <- data.frame(Isotropy = NA, rotation= NA, scale = NA, Lag = NA, Amount = NA, Dist = NA, estimator = NA, Method = NA, rejection = NA)

for(i in 1:3){
  if(i == 1){rot <- 1; s <- 1; iso <- "\u03B8 = 0, b = 1"}
  if(i == 2){rot <- 1; s <- 2; iso <- "\u03B8 = 0, b = 1.41"}
  if(i == 3){rot <- 2; s <- 2; iso <- "\u03B8 = \u03C0/4, b = 1.41"}
  for(l in 1:3){
    for(es in 1:3){
      for(m in 1:2){
        res.long.BA.sq <- rbind(res.long.BA.sq, data.frame(Isotropy = iso, rotation  = rotation[rot], scale = Scale[s], Amount = amount[2], Dist = dist[2], Lag = l, estimator = estimator[es], Method = method[m],  rejection = res_BA_sq[es,m,rot,s, l]))
      }
    }
  }
}

res.long.BA.sq$fill_col <- ifelse(res.long.BA.sq$rotation == "\u03B8 = 0", res.long.BA.sq$Method, "open")


# Amount 0.2, Dist 2
L1.D2.A2 <-  filter(res.long.BA.sq, Amount == 0.2, Dist == "N(5, 1)", Lag == 1)
PL1.D2.A2 <- ggplot(data = L1.D2.A2, mapping = aes(x = estimator, y = rejection, col = Method, shape = scale, fill = fill_col)) + geom_point(size = 2, alpha = 1, stroke = 1.2) + 
  ylab("Proportion of Rejections") + labs(title = expression(epsilon == 0.2 ~ ", N(5, 1)")) + theme_minimal(base_size = 10) + geom_hline(yintercept = 0.05) + ylim(0, 0.6) + 
  scale_shape_manual(values = c(21, 22))  + scale_fill_manual(name = "rotation", values = c("Subsampling" = "#56B4E9", "Blockpermutation" = "#D55E00", "open" = "white")) +
  scale_color_manual(values = c("Subsampling" = "#56B4E9", "Block permutation" = "#D55E00")) + theme(legend.position = "none")


L2.D2.A2 <-  filter(res.long.BA.sq, Amount == 0.2, Dist == "N(5, 1)", Lag == 2)
PL2.D2.A2 <- ggplot(data = L2.D2.A2, mapping = aes(x = estimator, y = rejection, col = Method, shape = scale, fill = fill_col)) + geom_point(size = 2, alpha = 1, stroke = 1.2) + 
  ylab("Proportion of Rejections")  + theme_minimal(base_size = 10) + geom_hline(yintercept = 0.05) + ylim(0, 0.6) + 
  scale_shape_manual(values = c(21, 22))  + scale_fill_manual(name = "rotation", values = c("Subsampling" = "#56B4E9", "Blockpermutation" = "#D55E00", "open" = "white")) +
  scale_color_manual(values = c("Subsampling" = "#56B4E9", "Block permutation" = "#D55E00")) + theme(legend.position = "none")

L3.D2.A2 <-  filter(res.long.BA.sq, Amount == 0.2, Dist == "N(5, 1)", Lag == 3)
PL3.D2.A2 <- ggplot(data = L3.D2.A2, mapping = aes(x = estimator, y = rejection, col = Method, shape = scale, fill = fill_col)) + geom_point(size = 2, alpha = 1, stroke = 1.2) + 
  ylab("Proportion of Rejections") + labs(title = expression(epsilon == 0.2 ~ ", N(5, 1)")) + theme_minimal(base_size = 10) + geom_hline(yintercept = 0.05) + ylim(0, 0.6) + 
  scale_shape_manual(values = c(21, 22))  + scale_fill_manual(name = "rotation", values = c("Subsampling" = "#56B4E9", "Blockpermutation" = "#D55E00", "open" = "white")) +
  scale_color_manual(values = c("Subsampling" = "#56B4E9", "Block permutation" = "#D55E00")) + theme(legend.position = "none") 

plot_grid(PL1.D2.A2, PL2.D2.A2, PL3.D2.A2, legend_plot, ncol = 4, rel_widths = c(2,2,2,1))
ggsave("Graphs/BA_quad.pdf", width = 22, height = 8, units = "cm")




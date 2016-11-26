# =========================================================
# Title: Soil faunal diversity in the Tibetan grasslands
#
# Author details: Author: Xin Jing 
# Contact details: jingxin0123@gmail.com
#
# First created: July 19, 2015
# Last modified: Nov. 26, 2016
#
# Script info: This script performs statistical analyses on soil faunal community data.  
# Data info: Data consists of soil faunal community and environmental data.
# Data was collected on the Tibetan Plateau in the growing seasaon of 2011.
#========================================================== 
# Removing objects from the envrionment
rm(list = ls())

#==========================================================
# Loading library
library(plyr) # Calls: dlply
library(dplyr) # Calls: select, mutate
library(reshape2) # Calls: melt
library(ggplot2) # Calls: ggplot
library(randomForest) # Calls: randomForest, varImpPlot, importance etc.
library(rfUtilities) # Calls: rf.significance
library(rfPermute) # Calls: rfPermute, rp.importance
library(SpatialPack)# modified.ttest

#==========================================================
# load functions

# Scatterplot matrices
# put (absolute) correlations on the upper panels,
# with size proportional to the correlations.
panel.cor <- function(x, y, digits = 2, prefix = "", cex.cor, ...)
{
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(0, 1, 0, 1))
  r <- cor(x, y)
  # r <- abs(cor(x, y))
  txt <- format(c(r, 0.123456789), digits = digits)[1]
  txt <- paste0(prefix, txt)
  if(missing(cex.cor)) cex.cor <- 1.5/strwidth(txt)
  text(0.5, 0.5, txt, cex = cex.cor * abs(r))
}

## Normality test for univariate variable
# Shapiro-Wilk normality test
normTest <- function (xvar){
  xvar=xvar + 1
  if(shapiro.test(xvar)$p.value > 0.05){
    return("normal distribution")
  }
  else {lxvar=log(xvar)
  if(shapiro.test(lxvar)$p.value > 0.05){
    return("log-transformation")
  }
  else{sxvar=sqrt(xvar)
  if(shapiro.test(sxvar)$p.value > 0.05){
    return("sqrt-transformation")
  }
  else{
    return("non-parametric test")}
  }
  }
}

#==========================================================
# Loading files
# fauna data
fauna.dat <- read.table("soil_fauna_data_final.txt", header=TRUE)
# climate data
clim.dat <- read.table("site_clim_worldClim.txt", header = TRUE)

# Data transformation
cols <- paste("bio", c(1:2, 5:11), sep = "")
clim.dat[cols] <- clim.dat[cols]/10

# Merging data
all.fauna.dat <- merge(fauna.dat, clim.dat, by = "Site")

# calculate nematode and arthropd biomass (g biomass C m-2)
all.fauna.dat <- mutate(all.fauna.dat, nema.biom = 0.1*(Aphelenchida + Dorylaimida + 
                                           Enoplida + Rhabditida +
                                           Tylenchida)/1000/1000*0.5)
all.fauna.dat <- mutate(all.fauna.dat, arth.biom = (Araneae*570 + Coleoptera*700 + 
                                     Diptera*450 + Hemiptera*340 + 
                                     Hymenoptera*500 + Isopoda*110 +
                                     Lepidoptera*1980)/1000/1000*0.5)

###########################################################
# Table 1
# BIO12 = Annual Precipitation
# BIO15 = Precipitation Seasonality (Coefficient of Variation)
# BIO1 = Annual Mean Temperature
# BIO4 = Temperature Seasonality (standard deviation *100)

df.table1 <- all.fauna.dat %>%
  select(Long, Lat, Ele.m, pH, SOC, SCNratio, bio12, bio15, bio1, bio4,
         plant.species.richness) 
names(df.table1)
# mean, min and max
summary(df.table1)[c(4, 1, 6), ]

# standard error
se.df <- function(x) sd(x)/length(x) 
apply(df.table1, 2, se.df)

############################################################
# Figure S2
df.fig.s2 <- all.fauna.dat %>%
  select(pH, SOC, bio12, bio15, bio1, bio4, 
         plant.species.richness, Nematoda, nema.biom, 
         Arthropoda, arth.biom)
names(df.fig.s2) <- c("pH", "SOC", "MAP", "SEP", "MAT", "SET",
                      "plantSR", "nema.ab", "nema.bio", "arth.ab", 
                      "arth.bio")
# pdf("./outputs/supplementary_Fig2.pdf")
pairs(df.fig.s2, lower.panel=panel.smooth, 
      upper.panel=panel.cor, pch=19, cex=1.5,
      col=rgb(0,0,100,50, maxColorValue=255))
# dev.off()

###########################################################
# Data for randomForest analysis
# select soil, climate and pland SR
sub.fauna.dat <- all.fauna.dat %>%
  select(pH, SOC, bio12, bio4, plant.species.richness, 
         Nematoda, Arthropoda, Long, Lat)
names(sub.fauna.dat) <- c("pH", "SOC", "MAP", "SET",
                      "PlantSR", "nema.abund", "arth.abund", 
                      "Long", "Lat")

# Check normality for each varibale 
nsub.fauna.dat <- length(names(sub.fauna.dat))
for (i in 1:nsub.fauna.dat){
  out <- normTest(sub.fauna.dat[,i])
  cat(names(sub.fauna.dat[i]), out, "\n", sep="+++")
}

# data transformation
sub.fauna.dat$SOC <- log(sub.fauna.dat$SOC)
sub.fauna.dat$SET <- log(sub.fauna.dat$SET)
sub.fauna.dat$nema.abund <- sqrt(sub.fauna.dat$nema.abund)
sub.fauna.dat$arth.abund <- sqrt(sub.fauna.dat$arth.abund)
sub.fauna.dat$PlantSR <- sqrt(sub.fauna.dat$PlantSR)

#############################
# randomForest

# Fig 3 and supplementary figs
# by nematodes
set.seed(123)
fit1 <- randomForest(nema.abund ~ SOC + MAP + SET + PlantSR + pH, ntree = 500, 
                    keep.forest = TRUE, importance = TRUE, data = sub.fauna.dat)
print(fit1) # view results 
# varImpPlot(fit1, main = "Variable importance on nematode abundance") # view results
imp1 <- importance(fit1) # importance of each predictor
impvar1 <- rownames(imp1)[order(imp1[, 1], decreasing = TRUE)]
plot(fit1, main = "Error Rates-Nematode Random Forest")

# perform significance test for the random forest model
set.seed(125)
rf.perm1 <- rf.significance(fit1, sub.fauna.dat[c("pH", "SOC", "MAP", "SET", "PlantSR")], 
                           nperm = 9999, ntree = 500)
rf.perm1
# Number of permutations:  9999 
# p-value:  0 
# Model signifiant at p = 0 
# Model R-square:  0.3541961 
# Random R-square:  -0.1634808 
# Random R-square variance:  0.01361567

# estimate significance of importance metrics
# https://rdrr.io/cran/rfPermute/man/rfPermute.html
set.seed(127)
fit1.rp <- rfPermute(nema.abund ~ pH + SOC + MAP + SET + PlantSR, 
          ntree = 500, nrep = 9999, # num.cores = 3,
          data = sub.fauna.dat)
fit1.rp
# Type of random forest: regression
# Number of trees: 500
# No. of variables tried at each split: 1
# 
# Mean of squared residuals: 569.3448
# % Var explained: 36.91

# plot the scaled importance distributions and highlight significant predictors
plot(rp.importance(fit1.rp, scale = TRUE))
imp.rp1 <- rp.importance(fit1.rp, scale = TRUE)
df.imp1 <- as.data.frame(imp.rp1[, 1:2])
df.imp1 <- mutate(df.imp1, var.nam = row.names(df.imp1))
df.imp1 <- mutate(df.imp1, sig.sign = ifelse (df.imp1[2] < 0.05, "sig", "neu"))
row.names(df.imp1) <- row.names(imp.rp1)
colnames(df.imp1) <- c("IncMSE", "IncMSE.pval", "var.nam", "sig.sign")

# by arthropods
set.seed(124)
fit2 <- randomForest(arth.abund ~ pH + SOC + MAP + SET + PlantSR, ntree = 500, 
                    keep.forest = TRUE, importance = TRUE, data = sub.fauna.dat)
print(fit2) # view results 
# varImpPlot(fit2, main = "Variable importance on arthropod abundance") # view results
imp2 <- importance(fit2) # importance of each predictor
impvar2 <- rownames(imp2)[order(imp2[, 1], decreasing = TRUE)]
plot(fit2, main = "Error Rates-Arthropod Random Forest")

# perform significance test for the random forest model
set.seed(126)
rf.perm2 <- rf.significance(fit2, sub.fauna.dat[c("pH", "SOC", "MAP", "SET", "PlantSR")], 
                            nperm = 9999, ntree = 500)
rf.perm2
# Number of permutations:  9999 
# p-value:  0.009 
# Model signifiant at p = 0.009 
# Model R-square:  0.1177793 
# Random R-square:  -0.1612064 
# Random R-square variance:  0.01327428 

# estimate significance of importance metrics
set.seed(128)
fit2.rp <- rfPermute(arth.abund ~ pH + SOC + MAP + SET + PlantSR, 
                     ntree = 500, nrep = 9999, #num.cores = 3,
                     data = sub.fauna.dat)
fit2.rp
# Type of random forest: regression
# Number of trees: 500
# No. of variables tried at each split: 1
# 
# Mean of squared residuals: 177.2887
# % Var explained: 11.58

# plot the scaled importance distributions and highlight significant predictors
plot(rp.importance(fit2.rp, scale = TRUE))
imp.rp2 <- rp.importance(fit2.rp, scale = TRUE)
df.imp2 <- imp.rp2[, 1:2]
df.imp2 <- as.data.frame(imp.rp2[, 1:2])
df.imp2 <- mutate(df.imp2, var.nam = row.names(df.imp2))
df.imp2 <- mutate(df.imp2, sig.sign = ifelse (df.imp2[2] < 0.05, "sig", "neu"))
row.names(df.imp2) <- row.names(imp.rp2)
colnames(df.imp2) <- c("IncMSE", "IncMSE.pval", "var.nam", "sig.sign")

# plot relative importance 
df.imp1$var.nam <- factor(df.imp1$var.nam, 
                          levels = df.imp1$var.nam[order(df.imp1$IncMSE, decreasing = TRUE)])
p1 <- ggplot(data = df.imp1, aes(x = var.nam, y = IncMSE, color = sig.sign)) + 
  geom_segment(data = df.imp1, 
               aes(x = var.nam, y = 0, xend = var.nam, 
                   yend = IncMSE), size = 10, lineend = "butt") +
  ylim(0, 16) +
  scale_color_manual(values = c("blue", "red")) +
  xlab("") + ylab("Increase in MSE (%)") +
  annotate("text", x = Inf, y = 15, label = "(a)",
           vjust = 1.5, size = 6) +
  scale_x_discrete(limit = c("PlantSR", "SOC", "SET", "MAP", "pH"),
    labels = c("Plant species richness",
               "Soil organic carbon",
               "Temperature seasonality",
               "Mean annual precipitation", "Soil pH")) +
  coord_flip() +
  theme_bw(base_size = 16) +
  theme(legend.position = "none")

# plot relative importance 
df.imp2$var.nam <- factor(df.imp2$var.nam, 
                          levels = df.imp2$var.nam[order(df.imp2$IncMSE, decreasing = TRUE)])
p2 <- ggplot(data = df.imp2, aes(x = var.nam, y = IncMSE, color = sig.sign)) + 
  geom_segment(data = df.imp2, 
               aes(x = var.nam, y = 0, xend = var.nam, 
                   yend = IncMSE), size = 10, lineend = "butt") +
  ylim(0, 16) +
  scale_color_manual(values = c("blue", "red")) +
  xlab("") + ylab("Increase in MSE (%)") +
  annotate("text", x = Inf, y = 15, label = "(b)",
           vjust = 1.5, size = 6) +
  scale_x_discrete(limit = c("pH", "SOC", "SET", "MAP", "PlantSR"),
                   labels = c("Soil pH",
                              "Soil organic carbon",
                              "Temperature seasonality",
                              "Mean annual precipitation", 
                              "Plant species richness")) +
  coord_flip() +
  theme_bw(base_size = 16) +
  theme(legend.position = "none")

# save figs
# pdf("./outputs/fig_3_importance_IncMSE.pdf")
gridExtra::grid.arrange(gridExtra::arrangeGrob(p1, p2, ncol = 1))
# dev.off()

#############################
# partial response plot
# http://luthuli.cs.uiuc.edu/~daf/courses/Optimization/Papers/2699986.pdf

# pdf("./outputs/rf_partial_plot.pdf", width = 12.5, height = 6.5)
op <- par(mfrow = c(2, 5), mar = c(4.5, 4.5, 4, 1.5))
for (i in seq_along(impvar1)) {
  partialPlot(fit1, sub.fauna.dat, impvar1[i], xlab = impvar1[i],
              ylab = expression(Nematode~abundance~
                                (number~of~individuals~m^{-2}~of~soil)),
              main = paste("Partial Dependence on", impvar1[i]),
              ylim = c(35, 60))
}

for (i in seq_along(impvar2)) {
  partialPlot(fit2, sub.fauna.dat, impvar2[i], xlab = impvar2[i],
              ylab = expression(Arthropod~abundance~
                                  (number~of~individuals~m^{-2}~of~soil)),
              main = paste("Partial Dependence on", impvar2[i]),
              ylim = c(20, 40))
}
par(op)
# dev.off()

#############################
# Error rates-random forest
# pdf("./outputs/supplementary_Fig3.pdf", width = 9, height = 5)
op <- par(mfrow = c(1, 2))
plot(fit1, cex = 2.5, lwd = 1.2,
     main = "Error Rates-Nematode")
text(x = 400, y = 1550, "(a)", cex = 1.5)
plot(fit2, cex = 2.5, lwd = 1.2, 
     main = "Error Rates-Arthropod")
text(x = 400, y = 300, "(b)", cex = 1.5)
par(op)
# dev.off()

#############################
# random forest for each taxa of nematodes and arthropods
# select nematode and arthropod community
df.sp <- all.fauna.dat %>%
  select(Site, Aphelenchida, Dorylaimida,
         Enoplida, Rhabditida,
         Tylenchida, Araneae,
         Coleoptera, Diptera,
         Hemiptera, Hymenoptera,
         Isopoda, Lepidoptera)
# data transformation
df.sp[-1] <- sqrt(df.sp[-1])

# select environmental variables
df.env <- all.fauna.dat %>%
  select(Site, pH, SOC, bio12, bio4, plant.species.richness)
names(df.env) <- c("Site", "pH", "SOC", "MAP", "SET", "PlantSR")
df.env$SOC <- log(df.env$SOC)
df.env$SET <- log(df.env$SET)
df.env$PlantSR <- sqrt(df.env$PlantSR)

# merge the two data sets
df.all <- merge(df.env, df.sp, by = "Site")

# random forest Table 4
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~
j <- 1
k <- 100
m <- 1000
for(i in 7:18) {
  print("#---------------------------------------------------------#")
  print(names(df.all[i]))
  # perform randomForest
  set.seed(j)
  fit <- randomForest(df.all[, i] ~ SOC + MAP + SET + PlantSR + pH, ntree = 500, 
                      keep.forest = TRUE, importance = TRUE, data = df.all)
  j <- j + 1
  # perform significance test for the random forest model
  set.seed(k)
  rf.fit <- rf.significance(fit, df.all[c("pH", "SOC", "MAP", "SET", "PlantSR")], 
                  nperm = 9999, ntree = 500)
  print(rf.fit)
  k <- k + 1
  # estimate significance of importance metrics
  set.seed(m)
  rp.fit <- rfPermute(df.all[, i] ~ pH + SOC + MAP + SET + PlantSR, 
                      ntree = 500, nrep = 9999, # num.cores = 3,
                      data = df.all)
  print(rp.importance(rp.fit, scale = TRUE))
  m <- m + 1
}

###########################################################
# fig. 1 bivariate plot
# reshape the data frame
df1 <- melt(sub.fauna.dat, id.vars = c("Long", "Lat", "nema.abund", "arth.abund"))
names(df1) <- c("Long", "Lat", "nema.abund", "arth.abund", "EnvVars", "EnvValues")
levels(df1$EnvVars) <- c("Soil pH", "Soil organic carbon", "Mean annual precipitation",
                                  "Temperature seasonality", "Plant species richness")
df2 <- melt(df1, id.vars = c("Long", "Lat", "EnvVars", "EnvValues"))
levels(df2$variable) <- c("Nematode abundance", "Arthropod abundance")

ggplot(data = df2, aes(x = EnvValues, y = value)) +
  geom_point(size = 3.5, shape = 1, colour = "blue") +
  geom_smooth(method = "lm", colour = "red", se = TRUE) +
  facet_grid(variable ~ EnvVars, scales = "free") +
  xlab("Biotic and abiotic variables") +
  ylab("Soil faunal abundance") +
  theme_bw(base_size = 16) +
  theme(panel.grid = element_blank(),
        legend.key = element_rect(colour = NA),
        strip.background = element_blank())
# ggsave("./outputs/fig_1.pdf", width = 13.34, height = 6.88)

#############################
# Modified t test "SpatialPack"
## Dutilleul, P. (1993), Modifying the t test for assessing the correlation between two spatial processes. Biometrics 49, 305–314.
dlply(df2, .(variable, EnvVars),
      function (x) modified.ttest(x$EnvValues, x$value, x[, 1:2])[c(1, 2, 4)])

###########################################################
# Fig. 2

# select environmental variables
df.fig.2.env <- all.fauna.dat %>%
  select(Site, pH, SOC, bio12, bio4, plant.species.richness)

# data transformation
cols <- c("SOC", "bio4")
df.fig.2.env[cols] <- log(df.fig.2.env[cols])
df.fig.2.env$plant.species.richness <- sqrt(df.fig.2.env$plant.species.richness)
# scaling and centering for ggplot
df.fig.2.env[-1] <- scale(df.fig.2.env[-1], center = TRUE, scale = TRUE)

# select nematode and arthropod community
df.fig.2.sp <- all.fauna.dat %>%
  select(Site, Aphelenchida, Dorylaimida,
         Enoplida, Rhabditida,
         Tylenchida, Araneae,
         Coleoptera, Diptera,
         Hemiptera, Hymenoptera,
         Isopoda, Lepidoptera)

# merge envrionment and community data
df.fig.2 <- merge(df.fig.2.env, df.fig.2.sp, by = "Site")

# reshaping data
df.fig.2.long <- melt(df.fig.2, id.vars = c("Site", "pH", "SOC", "bio12", "bio4", "plant.species.richness"))

# add a new colum
df.fig.2.long$taxa.raw <- ifelse(df.fig.2.long$variable %in% c("Tylenchida", "Aphelenchida",
                                                                 "Rhabditida", "Enoplida", "Dorylaimida"),"Nematode", "Arthropod")
df.fig.2.long$taxa.raw <- factor(df.fig.2.long$taxa.raw, ordered = TRUE,
                                  levels = c("Nematode", "Arthropod"))
names(df.fig.2.long) <- c("Site", "pH", "SOC", "bio12", "bio4", "plant.species.richness",
                           "taxa.fine", "abundance", "taxa.raw")
# reshaping data
df1 <- melt(df.fig.2.long, id.vars = c("Site", "taxa.fine", "taxa.raw", "abundance"))
levels(df1$variable) <- c("Soil pH", "Soil organic carbon", "Mean annual precipitation",
                          "Temperature seasonality", "Plant species richness")
df1$abundance <- sqrt(df1$abundance)

# ggplot
ggplot(df1, aes(x = abundance, y = value, colour = taxa.fine, shape = taxa.raw)) +
  geom_point(size = 2.5) +
  geom_smooth(method = "lm", se = FALSE) +
  facet_wrap(taxa.raw ~ variable, ncol = 5) +
  coord_flip() +
  xlab(expression(Abundance~(number~of~individuals~m^{-2}~of~soil))) +
  ylab("") +
  theme_bw(base_size = 15) +
  theme(panel.grid = element_blank(),
        legend.key = element_rect(colour = NA),
        strip.background=element_blank())
# ggsave("./outputs/Fig2.pdf", width = 13.348, height = 6.88)

#############################
# table 3

# linear mixed models
# select environmental variables
df.table3.env <- all.fauna.dat %>%
  select(Site, Long, Lat, pH, SOC, bio12, bio4, plant.species.richness)

# data transformation
cols <- c("SOC", "bio4")
df.table3.env[cols] <- log(df.table3.env[cols])
df.table3.env$plant.species.richness <- sqrt(df.table3.env$plant.species.richness)

# select nematode and arthropod community
df.table3.sp <- all.fauna.dat %>%
  select(Site, Aphelenchida, Dorylaimida,
         Enoplida, Rhabditida,
         Tylenchida, Araneae,
         Coleoptera, Diptera,
         Hemiptera, Hymenoptera,
         Isopoda, Lepidoptera)

# merge envrionment and community data
df.table3 <- merge(df.table3.env, df.table3.sp, by = "Site")

# reshaping data
df.table3.long <- melt(df.table3, id.vars = c("Site", "Long", "Lat", "pH", "SOC", "bio12", "bio4", "plant.species.richness"))

# add a new colum
df.table3.long$taxa.raw <- ifelse(df.table3.long$variable %in% c("Tylenchida", "Aphelenchida",
                                                                 "Rhabditida", "Enoplida", "Dorylaimida"),"Nematode", "Arthropod")
df.table3.long$taxa.raw <- factor(df.table3.long$taxa.raw, ordered = TRUE,
                                  levels = c("Nematode", "Arthropod"))
names(df.table3.long) <- c("Site", "Long", "Lat", "pH", "SOC", "MAP", "SET", "PlantSR",
                           "taxa.fine", "abundance", "taxa.raw")
df1 <- df.table3.long
df1$abundance <- sqrt(df1$abundance)

## Modified t test "SpatialPack"
## Dutilleul, P. (1993), Modifying the t test for assessing the correlation between two spatial processes. Biometrics 49, 305–314.
dlply(df1, .(taxa.raw, taxa.fine),
      function (x) modified.ttest(x$pH, x$abundance, x[, 2:3])[c(1, 2, 4)])
dlply(df1, .(taxa.raw, taxa.fine),
      function (x) modified.ttest(x$SOC, x$abundance, x[, 2:3])[c(1, 2, 4)])
dlply(df1, .(taxa.raw, taxa.fine),
      function (x) modified.ttest(x$MAP, x$abundance, x[, 2:3])[c(1, 2, 4)])
dlply(df1, .(taxa.raw, taxa.fine),
      function (x) modified.ttest(x$SET, x$abundance, x[, 2:3])[c(1, 2, 4)])
dlply(df1, .(taxa.raw, taxa.fine),
      function (x) modified.ttest(x$PlantSR, x$abundance, x[, 2:3])[c(1, 2, 4)])

# End of the script
###########################################################
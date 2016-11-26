# =========================================================
# Title: Soil faunal community diversity analyses
#
# Author details: Author: Xin Jing, Contact details: jingxin0123@gmail.com
# Script info: This script performs statistical analyses on soil faunal community data.  
# Data info: Data consists of soil faunal community and environmental data.
# Data was collected on the Tibetan Plateau between July and September in 2011. 
#========================================================== 
# Removing objects from the currently active envrionment
rm(list = ls())

#==========================================================
# Loading library
library(plyr) # Calls: dlply
library(dplyr) # Calls: select
library(reshape2) # Calls: melt
library(ggplot2) # Calls: ggplot
library(SpatialPack) # Calls: modified.ttest
library(spdep) # Calls: dnearneigh
library(vegan) # Calls: cca
library(lavaan) # Calls: sem
#==========================================================
# Pair-plot matrix
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
# Loading files and checking data structure
fauna.dat <- read.table("soil_fauna_data_dryad.txt", header=TRUE)
names(fauna.dat)
str(fauna.dat)

clim.dat <- read.table("site_clim_worldClim.txt", header = TRUE)
names(clim.dat)
str(clim.dat)

# Data transformation
cols <- paste("bio", c(1:2, 4:11), sep = "")
clim.dat[cols] <- clim.dat[cols]/10

# Merging data
all.fauna.dat <- merge(fauna.dat, clim.dat, by = "Site")
names(all.fauna.dat)

###########################################################
# Table 1
df.table1 <- all.fauna.dat %>%
  select(Long, Lat, Ele.m, pH, SOC, SCNratio, bio12, bio15, bio1, bio4,
         plant.species.richness)
names(df.table1)
summary(df.table1)[c(4, 1, 6), ]

############################################################
# Figure S1
df.fig.s1 <- all.fauna.dat %>%
  select(pH, SOC, SCNratio, bio12, bio15, bio1, bio4, 
         plant.species.richness, faunal.richness, faunal.abundance, faunal.evenness)
names(df.fig.s1) <- c("pH", "SOC", "SCN", "MAP", "SEP", "MAT", "SET",
                      "PlantSR", "FaunaSR", "FaunaA", "FaunaE")
# pdf("supplementary_Fig1.pdf")
pairs(df.fig.s1, lower.panel=panel.smooth, 
      upper.panel=panel.cor, pch=19, cex=1.5,
      col=rgb(0,0,100,50, maxColorValue=255))
# dev.off()

###########################################################
# Figure S2
##PCA
df.fig.s2 <- all.fauna.dat %>%
  select(pH, SOC, SCNratio, bio12, bio15, bio1, bio4, plant.species.richness, 
         faunal.richness, faunal.abundance, faunal.evenness, Long, Lat)
names(df.fig.s2) <- c("pH", "SOC", "SCN", "MAP", "SEP", "MAT", "SET",
                      "PlantSR", "FaunaSR", "FaunaA", "FaunaE", "Long", "Lat")
pca.abiotic <- prcomp(df.fig.s2[1:8], center=T, scale=T)
ls(pca.abiotic)
pca.abiotic$sdev^2
pca.abiotic$rotation
screeplot(pca.abiotic, type="line", main="screeplot")
summary(pca.abiotic)
# pdf("supplementary_Fig2.pdf")
biplot(pca.abiotic, choices = c(1, 4), cex=1.2, cex.axis=1.2, cex.lab=1.5, 
       col=c("red", "blue"), xlab="PC1 (50.6%)", ylab="PC4 (6.99%)")
# dev.off()

#############################
# Table 3, Table S1, and Table S2
# Check normality for each varibale 
ndf.fig.s2 <- length(names(df.fig.s2))
for (i in 1:ndf.fig.s2){
  out <- normTest(df.fig.s2[,i])
  cat(names(df.fig.s2[i]), out, "\n", sep="+++")
}
cols <- c("SOC", "SCN", "SEP", "SET", "FaunaSR", "FaunaA")
df.fig.s2[cols] <- log(df.fig.s2[cols])

# Table S1
#Make a matrix of coordinates (X and Y coordinates)
xyLatLon <- cbind(df.fig.s2$Long, df.fig.s2$Lat)
coords <- as.matrix(xyLatLon)
#Define neighbourhood
nb200 <- dnearneigh(coords, 0, 200, longlat = T)
plot(nb200, coords)
summary(nb200) #Percentage nonzero weights: 19.61806 

nb400 <- dnearneigh(coords, 0, 400, longlat = T)
plot(nb400, coords)
summary(nb400) #Percentage nonzero weights: 50.52083

nb800 <- dnearneigh(coords, 0, 800, longlat = T)
plot(nb800, coords)
summary(nb800) #Percentage nonzero weights: 92.36111

for (i in 1:ndf.fig.s2){
  out1 <- moran.test(df.fig.s2[, i], nb2listw(nb200, style = "W"))
  print(i); print(out1[c(2, 3)])
  out2 <- moran.test(df.fig.s2[, i], nb2listw(nb400, style = "W"))
  print(out2[c(2, 3)])
  out3 <- moran.test(df.fig.s2[, i], nb2listw(nb800, style = "W"))
  print(out3[c(2, 3)])
}

#Spatial weights, illustrated with coding style "W" (row standardized)
nb200.w<-nb2listw(nb200, glist=NULL, style="W", zero.policy=FALSE)
nb400.w<-nb2listw(nb400, glist=NULL, style="W", zero.policy=FALSE)
nb800.w<-nb2listw(nb800, glist=NULL, style="W", zero.policy=FALSE)

# Table S2
names(df.fig.s2)
# 1, 2, 4, 7, 8
ndf.fig.s2 <- 8
for (i in 1:ndf.fig.s2){
  ols <- lm(df.fig.s2$FaunaSR ~ df.fig.s2[, i])
  res <- ols$residuals
  out1 <- moran.test(res, listw = nb200.w, zero.policy=T); print (i); print(out1)
  out2 <- moran.test(res, listw = nb400.w, zero.policy=T); print(out2)
  out3 <- moran.test(res, listw = nb800.w, zero.policy=T); print(out3)
}

for (i in 1:ndf.fig.s2){
  ols <- lm(df.fig.s2$FaunaA ~ df.fig.s2[, i])
  res <- ols$residuals
  out1 <- moran.test(res, listw = nb200.w, zero.policy=T); print (i); print(out1)
  out2 <- moran.test(res, listw = nb400.w, zero.policy=T); print(out2)
  out3 <- moran.test(res, listw = nb800.w, zero.policy=T); print(out3)
}

for (i in 1:ndf.fig.s2){
  ols <- lm(df.fig.s2$FaunaE ~ df.fig.s2[, i])
  res <- ols$residuals
  out1 <- moran.test(res, listw = nb200.w, zero.policy=T); print (i); print(out1)
  out2 <- moran.test(res, listw = nb400.w, zero.policy=T); print(out2)
  out3 <- moran.test(res, listw = nb800.w, zero.policy=T); print(out3)
}

# Table3 and Table S2
df.table3 <- cbind(df.fig.s2, pca.abiotic$x)

# model selection
# Faunal richness
full.mod <- lm(FaunaSR ~ PC1 + PC2 + PC3 + PC4 +PC5 + PC6 + PC7 + PC8, data = df.table3)
mod1 <- update(full.mod, ~. - PC8); summary(mod1)
mod2 <- update(mod1, ~. - PC7); summary(mod2)
mod3 <- update(mod2, ~. - PC6); summary(mod3)
mod4 <- update(mod3, ~. - PC5); summary(mod4)
mod5 <- update(mod4, ~. - PC4); summary(mod5)
mod6 <- update(mod5, ~. - PC3); summary(mod6)
mod7 <- update(mod6, ~. - PC2); summary(mod7); anova(mod7)
AIC(full.mod, mod1, mod2, mod3, mod4, mod5, mod6, mod7)
res <- mod7$residuals
moran.test(res, listw = nb200.w, zero.policy=T)
moran.test(res, listw = nb400.w, zero.policy=T)
moran.test(res, listw = nb800.w, zero.policy=T)

full.mod <- lm(FaunaSR ~ pH + I(pH^2) + SOC + I(SOC^2) + SCN + I(SCN^2) + 
                 MAP + I(MAP^2) + SEP + I(SEP^2) + MAT + I(MAT^2) + SET + 
                 I(SET^2) + PlantSR + I(PlantSR^2), data = df.table3)
summary(full.mod)
mod1 <- update(full.mod, ~. - I(PlantSR^2)); summary(mod1)
mod2 <- update(mod1, ~. - I(SET^2)); summary(mod2)
mod3 <- update(mod2, ~. - I(MAT^2)); summary(mod3)
mod4 <- update(mod3, ~. - I(pH^2)); summary(mod4)
mod5 <- update(mod4, ~. - MAT); summary(mod5)
mod6 <- update(mod5, ~. - I(SCN^2)); summary(mod6)
mod7 <- update(mod6, ~. - SET); summary(mod7)
mod8 <- update(mod7, ~. - I(MAP^2)); summary(mod8)
mod9 <- update(mod8, ~. - I(SOC^2)); summary(mod9)
mod10 <- update(mod9, ~. - I(SEP^2)); summary(mod10)
mod11 <- update(mod10, ~. - SEP); summary(mod11)
mod12 <- update(mod11, ~. - SCN); summary(mod12)
mod13 <- update(mod12, ~. - SOC); summary(mod13)
mod14 <- update(mod13, ~. - PlantSR); summary(mod14); anova(mod14)
res <- mod14$residuals
moran.test(res, listw = nb200.w, zero.policy=T)
moran.test(res, listw = nb400.w, zero.policy=T)
moran.test(res, listw = nb800.w, zero.policy=T)

# Faunal abundance
full.mod <- lm(FaunaA ~ PC1 + PC2 + PC3 + PC4 +PC5 + PC6 + PC7 + PC8, data = df.table3)
summary(full.mod)
mod1 <- update(full.mod, ~. - PC3); summary(mod1)
mod2 <- update(mod1, ~. - PC6); summary(mod2)
mod3 <- update(mod2, ~. - PC2); summary(mod3)
mod4 <- update(mod3, ~. - PC7); summary(mod4)
mod5 <- update(mod4, ~. - PC8); summary(mod5)
mod6 <- update(mod5, ~. - PC5); summary(mod6); anova(mod6)
res <- mod6$residuals
moran.test(res, listw = nb200.w, zero.policy=T)
moran.test(res, listw = nb400.w, zero.policy=T)
moran.test(res, listw = nb800.w, zero.policy=T)

full.mod <- lm(FaunaA ~ pH + I(pH^2) + SOC + I(SOC^2) + SCN + I(SCN^2) + 
                 MAP + I(MAP^2) + SEP + I(SEP^2) + MAT + I(MAT^2) + SET + 
                 I(SET^2) + PlantSR + I(PlantSR^2), data = df.table3)
summary(full.mod)
mod1 <- update(full.mod, ~. - I(SCN^2)); summary(mod1)
mod2 <- update(mod1, ~. - MAT); summary(mod2)
mod3 <- update(mod2, ~. - I(MAP^2)); summary(mod3)
mod4 <- update(mod3, ~. - SCN); summary(mod4)
mod5 <- update(mod4, ~. - I(SEP^2)); summary(mod5)
mod6 <- update(mod5, ~. - I(pH^2)); summary(mod6)
mod7 <- update(mod6, ~. - SOC); summary(mod7)
mod8 <- update(mod7, ~. - I(PlantSR^2)); summary(mod8)
mod9 <- update(mod8, ~. - PlantSR); summary(mod9)
mod10 <- update(mod9, ~. - SEP); summary(mod10)
mod11 <- update(mod10, ~. - pH); summary(mod11)
mod12 <- update(mod11, ~. - I(SET^2)); summary(mod12)
mod13 <- update(mod12, ~. - SET); summary(mod13)
mod14 <- update(mod13, ~. - I(MAT^2)); summary(mod14); anova(mod14)
res <- mod14$residuals
moran.test(res, listw = nb200.w, zero.policy=T)
moran.test(res, listw = nb400.w, zero.policy=T)
moran.test(res, listw = nb800.w, zero.policy=T)

# Faunal evenness
full.mod <- lm(FaunaE ~ PC1 + PC2 + PC3 + PC4 +PC5 + PC6 + PC7 + PC8, data = df.table3)
summary(full.mod)
mod1 <- update(full.mod, ~. - PC4); summary(mod1)
mod2 <- update(mod1, ~. - PC2); summary(mod2)
mod3 <- update(mod2, ~. - PC6); summary(mod3)
mod4 <- update(mod3, ~. - PC5); summary(mod4)
mod5 <- update(mod4, ~. - PC3); summary(mod5)
mod6 <- update(mod5, ~. - PC8); summary(mod6)
mod7 <- update(mod6, ~. - PC7); summary(mod7); anova(mod7)
res <- mod7$residuals
moran.test(res, listw = nb200.w, zero.policy=T)
moran.test(res, listw = nb400.w, zero.policy=T)
moran.test(res, listw = nb800.w, zero.policy=T)

full.mod <- lm(FaunaE ~ pH + I(pH^2) + SOC + I(SOC^2) + SCN + I(SCN^2) + 
                 MAP + I(MAP^2) + SEP + I(SEP^2) + MAT + I(MAT^2) + SET + 
                 I(SET^2) + PlantSR + I(PlantSR^2), data = df.table3)
summary(full.mod)
mod1 <- update(full.mod, ~. - I(PlantSR^2)); summary(mod1)
mod2 <- update(mod1, ~. - I(MAP^2)); summary(mod2)
mod3 <- update(mod2, ~. - PlantSR); summary(mod3)
mod4 <- update(mod3, ~. - I(SCN^2)); summary(mod4)
mod5 <- update(mod4, ~. - MAT); summary(mod5)
mod6 <- update(mod5, ~. - SET); summary(mod6)
mod7 <- update(mod6, ~. - I(MAT^2)); summary(mod7)
mod8 <- update(mod7, ~. - I(pH^2)); summary(mod8)
mod9 <- update(mod8, ~. - SOC); summary(mod9)
mod10 <- update(mod9, ~. - I(SOC^2)); summary(mod10)
mod11 <- update(mod10, ~. - SCN); summary(mod11)
mod12 <- update(mod11, ~. - pH); summary(mod12)
mod13 <- update(mod12, ~. - I(SET^2)); summary(mod13)
mod14 <- update(mod13, ~. - I(SEP^2)); summary(mod14)
mod15 <- update(mod14, ~. - SEP); summary(mod15); anova(mod15)
res <- mod15$residuals
moran.test(res, listw = nb200.w, zero.policy=T)
moran.test(res, listw = nb400.w, zero.policy=T)
moran.test(res, listw = nb800.w, zero.policy=T)

###########################################################
# Fig. 2
df.fig2 <- all.fauna.dat %>%
  select(Long, Lat, pH, SOC, bio12, bio4, plant.species.richness, faunal.richness, faunal.abundance, faunal.evenness)
names(df.fig2) <- c("Long", "Lat", "pH", "SOC", "MAP", "SET", "PlantSR", "FaunaSR", "FaunaA", "FaunaE")
cols <- c("SOC", "SET", "FaunaSR", "FaunaA")
df.fig2[cols] <- log(df.fig2[cols])

df.fig2.long <- melt(df.fig2, id.vars = c("Long", "Lat", "FaunaSR", "FaunaA", "FaunaE"))
names(df.fig2.long) <- c("Long", "Lat", "FaunaSR", "FaunaA", "FaunaE", "EnvVars", "EnvValues")
levels(df.fig2.long$EnvVars) <- c("Soil pH", "Soil organic carbon", "Mean annual precipitation",
                                   "Temperature seasonality", "Plant species richness")
df.fig2.long.long <- melt(df.fig2.long, id.vars = c("Long", "Lat", "EnvVars", "EnvValues"))
names(df.fig2.long.long)
levels(df.fig2.long.long$variable) <- c("Soil faunal richness", "Soil faunal abundance",
                                        "Soil faunal evenness")
ggplot(data = df.fig2.long.long, aes(x = EnvValues, y = value)) +
  geom_point(size = 1.5, colour = "blue") +
  geom_smooth(method = "lm", colour = "red", se = TRUE) +
  facet_grid(variable ~ EnvVars, scales = "free") +
  xlab("Biotic and abiotic variables") +
  ylab("Soil faunal diversity") +
  theme_bw(base_size = 14) +
  theme(panel.grid = element_blank(),
        legend.key = element_rect(colour = NA),
        strip.background = element_blank())
# A4 landscape fig2.pdf

## Modified t test "SpatialPack"
## Dutilleul, P. (1993), Modifying the t test for assessing the correlation between two spatial processes. Biometrics 49, 305–314.
dlply(df.fig2.long.long, .(variable, EnvVars),
      function (x) modified.ttest(x$EnvValues, x$value, x[, 1:2])[c(1, 2, 4)])

###########################################################
# Fig. S3
names(all.fauna.dat)
df.fig.s3 <- all.fauna.dat %>%
  select(Long, Lat, pH, SOC, bio12, bio4, plant.species.richness, Nematoda, Arthropoda)
names(df.fig.s3) <- c("Long", "Lat", "pH", "SOC", "MAP", "SET", "PlantSR", "Nematode", "Arthropod")
cols <- c("SOC", "SET")
df.fig.s3[cols] <- log(df.fig.s3[cols])
cols <- c("Nematode", "Arthropod")
df.fig.s3[cols] <- sqrt(df.fig.s3[cols])
df.fig.s3.long <- melt(df.fig.s3, id.vars = c("Long", "Lat", "Nematode", "Arthropod"))
names(df.fig.s3.long)
levels(df.fig.s3.long$variable) <- c("Soil pH", "Soil organic carbon", "Mean annual precipitation",
                                     "Temperature seasonality", "Plant species richness")
names(df.fig.s3.long) <- c("Long", "Lat", "Nematode", "Arthropod", "EnvVars", "EnvValues")
df.fig.s3.long.long <- melt(df.fig.s3.long, id.vars = c("Long", "Lat", "EnvVars", "EnvValues"))
names(df.fig.s3.long.long)
levels(df.fig.s3.long.long$variable) <- c("Nematode abundance", "Arthropod abundance")

ggplot(data = df.fig.s3.long.long, aes(x = EnvValues, y = value)) +
  geom_point(size = 1.5, colour = "blue") +
  geom_smooth(method = "lm", colour = "red", se = TRUE) +
  facet_grid(variable ~ EnvVars, scales = "free") +
  xlab("Biotic and abiotic variables") +
  ylab("Soil faunal abundance") +
  theme_bw(base_size = 14) +
  theme(panel.grid = element_blank(),
        legend.key = element_rect(colour = NA),
        strip.background = element_blank())
# size 6.27 * 11.69 .pdf
## Modified t test "SpatialPack"
## Dutilleul, P. (1993), Modifying the t test for assessing the correlation between two spatial processes. Biometrics 49, 305–314.
dlply(df.fig.s3.long.long, .(variable, EnvVars),
      function (x) modified.ttest(x$EnvValues, x$value, x[, 1:2])[c(1, 2, 4)])

###########################################################
# Fig. s4
# selecting environmental variables
df.fig.s4.env <- all.fauna.dat %>%
  select(Site, pH, SOC, bio12, bio4, plant.species.richness)
# data transformation
cols <- c("SOC", "bio4")
df.fig.s4.env[cols] <- log(df.fig.s4.env[cols])
df.fig.s4.env[-1] <- scale(df.fig.s4.env[-1], center = TRUE, scale = TRUE)
# selecting nematode and arthropod species
df.fig.s4.sp <- all.fauna.dat %>%
  select(Site, Tylenchida, Aphelenchida,
         Rhabditida, Enoplida, Dorylaimida,
         Coleoptera, Diptera, Lepidoptera,
         Hymenoptera, Acariformes, Hemiptera,
         Araneae, Isopoda, Collembola)
# merging envrionment and community data
df.fig.s4 <- merge(df.fig.s4.env, df.fig.s4.sp, by = "Site")
names(df.fig.s4)
# reshaping data
df.fig.s4.long <- melt(df.fig.s4, id.vars = c("Site", "pH", "SOC", "bio12", "bio4", "plant.species.richness"))
levels(df.fig.s4.long$variable)
# add a new colum
df.fig.s4.long$taxa.raw <- ifelse(df.fig.s4.long$variable %in% c("Tylenchida", "Aphelenchida",
                                                     "Rhabditida", "Enoplida", "Dorylaimida"),"Nematode", "Arthropod")
df.fig.s4.long$taxa.raw <- factor(df.fig.s4.long$taxa.raw, ordered = TRUE,
                                  levels = c("Nematode", "Arthropod"))
names(df.fig.s4.long) <- c("Site", "pH", "SOC", "bio12", "bio4", "plant.species.richness",
                           "taxa.fine", "abundance", "taxa.raw")
# reshaping data
df1 <- melt(df.fig.s4.long, id.vars = c("Site", "taxa.fine", "taxa.raw", "abundance"))
levels(df1$variable) <- c("Soil pH", "Soil organic carbon", "Mean annual precipitation",
                          "Temperature seasonality", "Plant species richness")
# ggplot
ggplot(df1, aes(x = abundance, y = value, shape = taxa.raw, colour = taxa.fine)) +
  geom_point(size = 2.5) +
  geom_smooth(method = "lm", se = FALSE) +
  facet_wrap(taxa.raw ~ variable, ncol = 5) +
  coord_flip() +
  xlab("Abundance (individuals/m2)") +
  ylab("Values of variables") +
  theme_bw(base_size = 12) +
  theme(panel.grid = element_blank(),
        legend.key = element_rect(colour = NA),
        strip.background=element_blank())
# size 6.27 * 11.69 .pdf

###########################################################
# Fig s5
#======================================
# CCA analysis using vegan packages
names(all.fauna.dat)
df.fig.s5.sp <- all.fauna.dat[c(13:31)]
df.fig.s5.env <- all.fauna.dat %>%
  select(pH, SOC, SCNratio, bio12, bio15, bio1, bio4, plant.species.richness)
names(df.fig.s5.env) <- c("pH", "SOC", "SCN", "MAP", "SEP", "MAT", "SET", "PlantSR")
mod <- cca(df.fig.s5.sp ~ pH + SOC + SCN + MAP + SEP + MAT +
             SET + PlantSR, df.fig.s5.env)
ls(mod)
summary(mod)
anova(mod)
# model selection
ordistep(mod, perm.max = 999) 
mod <- cca(df.fig.s5.sp ~ SCN + MAP, df.fig.s5.env)
summary(mod)
cca.res<-summary(mod)
cca.sites <-data.frame(cca.res$sites)
ord_df<-data.frame(CCA1=cca.sites$CCA1,CCA2=cca.sites$CCA2)
exp<-cca.res$concont
exp<-data.frame(exp$importance)
cca.species<-data.frame(cca.res$species)
cca.species<-data.frame(Cca1=cca.species$CCA1,Cca2=cca.species$CCA2,species=rownames(cca.species))
cca.benthos<-data.frame(cca.res$biplot)
cca.benthos<-data.frame(cca1=cca.benthos$CCA1,cca2=cca.benthos$CCA2, Species = rownames(cca.benthos))

ggplot(ord_df) +
  geom_point(mapping = aes(x=CCA1, y=CCA2),size = 5, 
             # color=rgb(199,124,255, maxColorValue = 255))+
             color="white")+
  xlim(-3, 3) + ylim(-3, 3) +
  geom_text(data = cca.species, 
            aes(x = Cca1, y = Cca2+0.2, label = species),
            size = 5, alpha=0.4, colour = "blue") +
  geom_point(data = cca.species, 
             aes(x = Cca1, y = Cca2, label = species),
             pch = 17, alpha=0.4, size = 5,colour = "blue") +
  geom_segment(data = cca.benthos,
               aes(x = 0, xend = 2.5*cca1, y = 0, yend = 2.5*cca2),
               arrow = arrow(length = unit(0.5, "cm")), size=0.6, colour = "red") +
  geom_text(data = cca.benthos, 
            aes(x = cca1*3, y = cca2*3, label = Species),
            size = 3, colour = "red") +
  labs(list(title = NULL, x = paste("CCA1 (",round(exp[2,1]*100,digits=1),"%)"), 
            y = paste("CCA2 (",round(exp[2,2]*100,digits=1),"%)"))) +
  geom_hline(yintercept = 0, linetype = 2, colour = "gray70") +
  geom_vline(xintercept = 0, linetype = 2, colour = "gray70") +
  theme_bw()+
  #theme(panel.border = element_rect(colour = "white", size=1))+
  theme(axis.title.x = element_text(face="bold", size=14))+
  theme(axis.title.y = element_text(face="bold", size=14))+
  theme(axis.text.x = element_text(size=14,color="black"))+
  theme(axis.text.y = element_text(size=14,color="black"))+
  theme(axis.line.x=element_blank())+
  theme(axis.line.y=element_blank())+
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank())
# ggsave("supplementary_Figs5.pdf", width = 6, height = 6)
###########################################################
# Fig. 3 Structural equation models
names(all.fauna.dat)
df.fig3 <- all.fauna.dat %>%
  select(pH, SOC, SCNratio, bio12, bio15, bio1, bio4, plant.species.richness, faunal.richness, faunal.abundance, faunal.evenness)
names(df.fig3) <- c("pH", "SOC", "SCN", "MAP", "SEP", "MAT", "SET", "PlantSR",
                    "FaunaSR", "FaunaA", "FaunaE")
# data transformation
cols <- c("SOC", "SET", "FaunaSR", "FaunaA")
df.fig3[cols] <- log(df.fig3[cols])
# df.fig3 <- scale(df.fig3, center = TRUE, scale = TRUE)
# df.fig3 <- as.data.frame(df.fig3)

#======================================
model1 <- '
# regressions
pH ~ MAP
# SOC ~ MAP + pH
PlantSR ~ MAP + SET
FaunaSR ~ MAP + PlantSR + pH
# residual correlations
MAP ~~ SET
'
#fit SEM
fit1 <- sem(model1, df.fig3, fixed.x = FALSE)
# fitting variables to the same scale
varTable(fit1)
df.fig3$PlantSR <- df.fig3$PlantSR/10
df.fig3$MAP <- df.fig3$MAP/100
df.fig3$SCN <- df.fig3$SCN/10
df.fig3$SET <- df.fig3$SET*10
fit1 <- sem(model1, df.fig3, fixed.x = FALSE)
summary(fit1, standardized = TRUE, fit.measures = TRUE, rsquare = TRUE)
fitMeasures(fit1)
#############################
model2 <- '
# regressions
pH ~ MAP
# SOC ~ MAP + pH
PlantSR ~ MAP + SET
FaunaA ~ MAP + PlantSR + pH
# residual correlations
MAP ~~ SET
'
#fit SEM
fit2 <- sem(model2, df.fig3, fixed.x = FALSE)
summary(fit2, standardized = TRUE, fit.measures = TRUE, rsquare = TRUE)
fitMeasures(fit2)
#############################
model3 <- '
# regressions
pH ~ MAP
PlantSR ~ MAP + SET
FaunaE ~ MAP + pH + PlantSR
# residual correlations
MAP ~~ SET
'
#fit SEM
fit3 <- sem(model3, df.fig3, fixed.x = FALSE)
summary(fit3, standardized = TRUE, fit.measures = TRUE, rsquare = TRUE)
fitMeasures(fit3)

model4 <- '
# measurement model
climate =~ SET + SEP
soil =~ pH + SOC + SCN
# regressions
PlantSR ~ climate
FaunaSR ~ soil + PlantSR
# residual correlations
climate ~~ soil
'
fit4 <- sem(model4, df.fig3, fixed.x = FALSE, ridge = 1e-05)
summary(fit4, standardized = TRUE, fit.measures = TRUE, rsquare = TRUE)

model5 <- '
# measurement model
climate =~ SET + SEP
soil =~ pH + SOC + SCN
# regressions
PlantSR ~ climate
FaunaA ~ soil + PlantSR
# residual correlations
climate ~~ soil
'
fit5 <- sem(model5, df.fig3, fixed.x = FALSE, ridge = 1e-05)
summary(fit5, standardized = TRUE, fit.measures = TRUE, rsquare = TRUE)

model6 <- '
# measurement model
climate =~ SET + SEP
soil =~ pH + SOC + SCN
# regressions
PlantSR ~ climate
FaunaE ~ soil + PlantSR
# residual correlations
climate ~~ soil
'
fit6 <- sem(model6, df.fig3, fixed.x = FALSE, ridge = 1e-05)
summary(fit6, standardized = TRUE, fit.measures = TRUE, rsquare = TRUE)

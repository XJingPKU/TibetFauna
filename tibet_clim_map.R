rm(list = ls())
###########################################################
library(raster)
# Read point locations
locs <- read.csv("./zhao_ke_map/tibet_site_xy.csv")
head(locs)
coordinates(locs) <- c("Long", "Lat")  # set spatial coordinates



world.clim.MAT <- raster("bio_30s_esri/bio_30s_bil/bio_1.bil")
world.clim.MAP <- raster("bio_30s_esri/bio_30s_bil/bio_12.bil")

# Define spatial projection
crs.geo <- CRS("+proj=longlat +ellps=WGS84 +datum=WGS84")  # geographical, datum WGS84
projection(world.clim.MAT) <- crs.geo
projection(world.clim.MAP) <- crs.geo

# Define spatial projection
proj4string(locs) <- crs.geo
summary(locs)
head(locs)

# crop
newext <- c(78, 103, 26, 40)
world.clim.MAT.c <- crop(world.clim.MAT, newext)
world.clim.MAT.c <- world.clim.MAT.c / 10

world.clim.MAP.c <- crop(world.clim.MAP, newext)

# plot it
pdf("./zhao_ke_map/site_clim.pdf", width = 8.27, height = 11.69)
par(mfcol = c(2, 1))
plot(world.clim.MAT.c, cex = 2.0,
     ylab = "Latitude (N)",
     xlab = "Longitude (E)",
     main = expression("Mean annual temperature ("~degree~C~")"))
# add site locations
plot(locs, pch = 2, add = T, col = "red")

plot(world.clim.MAP.c, cex = 2.0,
     ylab = "Latitude (N)",
     xlab = "Longitude (E)",
     main = expression(paste("Mean annual precipitation (mm yr"^"-1", ")")))
# add site locations
plot(locs, pch = 2, add = T, col = "red")
dev.off()

###########################################################
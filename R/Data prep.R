# Data preparation for habitat analysis

library(amt)








set.seed(106)
habitat <- readOGR(dsn = ".", layer = "HabitatLayer2")
r <- raster(habitat,resolution=c(5,5))
habitat$Category<- as.integer(habitat$Category)
ras <- rasterize(habitat,r,"Category", background=0)
ras[is.na(ras[])]<- 0
names(ras) <-"Category"
ras$Habitat <- as.factor(ras$Habitat, labels=c("Open","Urban","Woodland"))
plot(ras)


#remove within release pen movements
atlas2018 <- read.csv("Atlas2018egoallo.csv")
locs <- atlas2018[atlas2018$Day==T,] # daytime only
coordinates(locs) <- ~medianEast+medianNort
crs(locs) <- crs(habitat)
rp <- readOGR(dsn = ".", layer = "ReleasePen2")
a <-over(locs,rp)
locs <- locs[is.na(a$fid),] # removes points in the release pen






hr <- atlas2018[atlas2018$Day==T,c("medianEast","medianNort")] # daytime only
# ATLAS data
library(sp)
library(adehabitatHR)
library(adehabitatHS)
library(rgdal)
library(raster)
library(adehabitatMA)
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

#coordinates(locs) <- ~medianEast+medianNort
#coordinates(hr) <- ~medianEast+medianNort
#crs(hr) <- crs(habitat)

#raster over habitat polygons

r <- raster(habitat,resolution=c(5,5))
habitat$Category<- as.integer(habitat$Category)
ras <- rasterize(habitat,r,"Category", background=0)
ras[is.na(ras[])]<- 0
names(ras) <-"Category"
#ras$Habitat <- as.factor(ras$Habitat, labels=c("Open","Urban","Woodland"))
plot(ras)
image(ras)
#proj4string(locs) == proj4string(map)


# create HR of all inds combined to get an area to create random coordinates for each bird.

pc <- mcp(locs,percent=99)
plot(pc)

#writeOGR(obj=pc, dsn=".", layer="99thHR", driver="ESRI Shapefile")
str(pc)
rr <- mask(ras,pc)
#plot(rr)


#points(locs, col=as.numeric(slot(locs, "data")[,1]))
#availmap <- as(rr,"SpatialPixelsDataFrame")
map <- as(rr,"SpatialPixelsDataFrame")

#dat <- list(locs,map)
#locs$Bird <- factor(locs$Bird)
#cp <- count.points(locs, map)

hab <- slot(map, "data")[,1]
hab[is.na(hab)]<- 0
av <- factor(hab,labels=c("Open", "Urban","Woodland"))
(tav <- table(av))
slot(map,"data")[,1]<-as.numeric(av)


#make available dataset

sn <- data.frame(table(locs$Bird)) # number of obs for each bird

#permute 10000 times? to get actual average for available. 
avail<- NULL
for(i in 1:4){ # for each bird
    s <- NULL
    for(j in 1:10000){ # permute values
        samples <- data.frame(ID = sn[i,1], spsample(pc,sn[i,2],type="random"))
        coordinates(samples) <- ~x+y
        crs(samples) <- crs(map)
        samples$ID <- factor(samples$ID)
        samp.mapped <- join(samples, map)
        samp.mapped<- factor(samp.mapped, levels=c(1,2,3))
        d <- as.data.frame.matrix(table(slot(samples,"data")[,1],samp.mapped))
        s <- rbind(s,d)
        print(paste("bird: ", sn[i,1], "iteration: ", j))
    }
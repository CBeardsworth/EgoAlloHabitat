# ATLAS data
library(sp)
library(adehabitatHR)
library(adehabitatHS)
library(rgdal)
library(raster)
library(adehabitatMA)
library(tidyverse)
library(lattice)
library(lubridate)


setwd("D:/1_PhD/EgoAllo/EgoAlloHabitat")
library(raster)
library(lubridate)
library(amt)
library(tidyverse)
library(cowplot)

set.seed(106)

#remove within release pen movements
dat <- read.csv("Data/Atlas2018egoallo.csv")%>%
  filter(!is.na(`medianEast`)) %>%
  select(id ="Bird", sex = "Sex.x", mass = "Mass.25.07", x = "medianEast", y = "medianNort",
         t = "RealTime", strategy = "Strat")
dat$t <- as.POSIXct(dat$t, tz="GMT")
str(dat)

# Investigate 1 individual
dat_1 <- dat %>% filter(id == 1090)
dat_1 <- amt::make_track(dat_1, x, y, t, crs = sp::CRS("+init=epsg:27700"))


summarize_sampling_rate(dat_1)#check sampling rate of individuals, in case this differs between individuals we can resample (below)

stps <- amt::track_resample(dat_1, rate = minutes(10), tolerance = minutes(1)) %>% #resample as in a later example with more individuals, some are sampled at 10 minute intervals. Minutes() comes from lubridate and creates a class of 'Period' which can be passed to track_resample(). Different common time units (hour, second etc can be used)
  filter_min_n_burst(min_n = 3) %>% #choose to keep only bursts (subsets of the track with constant sampling rate within the specified tolerance)
  steps_by_burst() %>% #changes the data from a point representation into a step representation (step length, turn angles)
  time_of_day(include.crepuscule = FALSE) #wrapper around maptools::sunriset and maptools::crepuscule to determine day or night and also dawn/dusk if include.crepuscule = T

str(stps, width = 80, strict.width = "no", nchar.max = 80, give.attr = FALSE)# for each step start (x1_, y1_) and end (x2_, y2_) coordinates are given, start and end time (t1_, t2_), step length (sl_ (in CRS units)), turning angles (ta_ in degrees), time difference (dt_) and burst (burst_) to which the step belongs

habitat <- readOGR(dsn = "Data", layer = "HabitatLayer2")
habitat$Category<- as.integer(habitat$Category)
r <- raster(habitat,resolution=c(5,5))
ras <- rasterize(habitat,r,"Category", background=0)
ras[is.na(ras[])]<- 0
names(ras) <-"Category"

# Before modelling space use and habitat selection we should perform some exploratory data analysis based on step length and turning angles in different habitat types. 
eda1 <- stps %>% 
  extract_covariates(ras, where = "end") %>% #extract the covariate values at the start point of each step. Habitat selection may require 'end of step' choices and movement processes 'start'.  If covariates are extracted at the end of the step then they are typically included in the model as main effects to answer : how do covariates influence where the animal moves? If they are extracted at the beginning of the step then they are typically included as an interaction within the movement chars e.g. step length to ask e.g. :do animals move faster/more directed if they start in a given habitat? Covariate values at the start and end of a step can also be included in the model as an interaction with each other to test hypotheses of the type: are animals more likely to stay in a given habitat if they are already in that habitat?
  mutate(habitat = factor(Category, levels = c(0, 1, 2), labels = c("Open","Urban","Woodland"))) 

p1 <- eda1 %>% select(habitat, tod = tod_end_, sl_, ta_) %>%
    gather(key, val, -habitat, -tod) %>%
    filter(key == "sl_") %>%
    ggplot(., aes(val, group = tod, fill = tod)) + geom_density(alpha = 0.5) +
    facet_wrap(~ habitat, nrow = 2) +
    xlab("Step length [m]") + theme_light() +
    ylab("Density") +
    theme(legend.title = element_blank())

p2 <- eda1 %>% select(habitat, tod = tod_end_, sl_, ta_) %>%
    gather(key, val, -habitat, -tod) %>%
    filter(key == "ta_") %>%
    ggplot(., aes(val, group = tod, fill = tod)) + geom_density(alpha = 0.5) +
    facet_wrap(~ habitat, nrow = 2) +
    xlab("Turn angle") + theme_light() +
    theme(legend.title = element_blank(),
          axis.title.y = element_blank())

pg1 <- plot_grid(
  p1 + theme(legend.position = "none"),
  p2 + theme(legend.position = "none"), rel_widths = c(1, 1)
)

leg <- get_legend(p1)
plot_grid(pg1, leg, rel_widths = c(1, 0.1))

#### fit the iSSF model ####

#make the random steps either through sampling from the observed turn step-length and turn angle distribution to the observed step lengths (traditional SSF) or fitting parametric distribution to observed step lengths (either negative exponential, a half normal, a log normal or a gamme Avgar et al 2016 appendix 2) and turn angles (von mises) (iSSF)
m1 <- stps %>% amt::random_steps(n = 9) %>% #fit gamma distribution to step lengths and von mises distribution to turn angles using maximum  likelihood, use these distributions to generate 9 random steps to each observed step. The more steps, the lower the error but higher computational burden. Creates new column (step_id_ that identifies different strata (see model below))
  amt::extract_covariates(ras, where = "end") %>% #extract (environmental) covariates at the endpoint of each step
  amt::time_of_day(include.crepuscule = FALSE) %>%
  mutate(log_sl_ = log(sl_)) %>%
  amt::fit_issf(case_ ~ factor(Category) + log_sl_ + # include environmental covariate (wet) and log of the step length(log_sl_) as a modifier of the shape parameter of the underlying gamma distribution. 
                  strata(step_id_))#each step is paired with several control steps that form a stratum

#including cosines of the turning angles and their interaction with dat would modify the concentration parameter of the underlying von mises distribution for the turning angles and allow the degree of directional persistence to depend on time of day. (I.e. leave woodland in morning, return in the evening). Then use AIC to determine best model. 



# m1 <- d1 %>% amt::fit_issf(case_ ~ wet + sl_ + wet:tod_end_+ sl_:tod_end_ + strata(step_id_)) #wrapper to survival::clogit
# 
# m1 <- d1 %>% amt::fit_issf(case_ ~ wet + log_sl_ + sl_ + wet:tod_end_+ log_sl_:tod_end_ + sl_:tod_end_ + strata(step_id_))
# 
# AIC(m1$model)
summary(m1)

s <- summary(m1$model)$coefficients
s 

#This tells us :
# a) There is evidence to suggest that the animal prefers forested wetlands over other landuse classes, 
# b) there is no difference in Ricky's preference for wetlands between day and night, 
# c) There is evidence to modify the shape of the gamma distribution fit to the observed step lengths (through the log of the step length), 
# d) the modification of the shape parameter should be done separately for day and night, indicating that expected movement speeds differ between day and night
    
locs <- atlas2018[atlas2018$Day==T,] # daytime only
coordinates(locs) <- ~medianEast+medianNort
crs(locs) <- crs(habitat)
#rp <- readOGR(dsn = ".", layer = "ReleasePen2")
#a <-over(locs,rp)
#locs <- locs[is.na(a$fid),] # removes points in the release pen

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
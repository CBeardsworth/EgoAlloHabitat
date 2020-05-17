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
dat <- read.csv("Atlas2018egoallo.csv")%>%
    filter(!is.na(`medianEast`)) %>%
    select(x = "medianEast", y = "medianNort",
           t = "RealTime", strategy = "Strat")%>%
    filter(id == 1090)
dat$t <- as.POSIXct(dat$t)
str(dat)
dat <- amt::make_track(dat, x, y, t, crs = sp::CRS("+init=epsg:27700"))

amt::summarize_sampling_rate(dat)

stps <- amt::track_resample(dat, rate = lubridate::minutes(10), tolerance = minutes(1)) %>%
    amt::filter_min_n_burst(min_n = 3) %>% amt::steps_by_burst() %>%
    amt::time_of_day(include.crepuscule = TRUE)

str(stps, width = 80, strict.width = "no", nchar.max = 80, give.attr = FALSE)

eda1 <- stps %>% amt::extract_covariates(ras, where = "start") %>%
    mutate(ras = factor(Category, levels = c(0, 1,2), labels = c("Open","Urban","Woodland")))

p1 <- eda1 %>% select(ras, tod = tod_end_, sl_, ta_) %>%
    gather(key, val, -ras, -tod) %>%
    filter(key == "sl_") %>%
    ggplot(., aes(val, group = tod, fill = tod)) + geom_density(alpha = 0.5) +
    facet_wrap(~ ras, nrow = 2) +
    xlab("Step length [m]") + theme_light() +
    ylab("Density") +
    theme(legend.title = element_blank())

p2 <- eda1 %>% select(ras, tod = tod_end_, sl_, ta_) %>%
    gather(key, val, -ras, -tod) %>%
    filter(key == "ta_") %>%
    ggplot(., aes(val, group = tod, fill = tod)) + geom_density(alpha = 0.5) +
    facet_wrap(~ ras, nrow = 2) +
    xlab("Turn angle") + theme_light() +
    theme(legend.title = element_blank(),
          axis.title.y = element_blank())

library(cowplot)
pg1 <- plot_grid(
    p1 + theme(legend.position = "none"),
    p2 + theme(legend.position = "none"), rel_widths = c(1, 1)
)

leg <- get_legend(p1)
plot_grid(pg1, leg, rel_widths = c(1, 0.1))


m1 <-stps %>% amt::random_steps(n = 9) %>%
    amt::extract_covariates(ras) %>%
    amt::time_of_day(include.crepuscule = FALSE) %>%
    mutate(log_sl_ = log(sl_)) -> d1

m1 <- d1 %>% amt::fit_issf(case_ ~ Category + sl_ + Category:tod_end_+ sl_:tod_end_ + amt::strata(step_id_))
m1 <- d1 %>% amt::fit_issf(case_ ~ Category + log_sl_ + Category:tod_end_+ log_sl_:tod_end_ + amt::strata(step_id_))
m1 <- d1 %>% amt::fit_issf(case_ ~ Category + log_sl_ + sl_ + Category:tod_end_+ log_sl_:tod_end_ + sl_:tod_end_ + strata(step_id_))

AIC(m1$model)
summary(m1)


s <- summary(m1$model)$coefficients
s

print(xtable::xtable(s, digits = 4,
                     type = "latex",
                     caption.placement = "top"))

shape <- sl_shape(m1)
scale <- sl_scale(m1)

shape_adj_day <- amt::adjust_shape(shape, coef(m1)["log_sl_"])
shape_adj_night <- amt::adjust_shape(shape, coef(m1)["log_sl_"]) +
    coef(m1)["log_sl_:tod_end_night"]

scale_adj_day <- amt::adjust_scale(scale, coef(m1)["sl_"])
scale_adj_night <- amt::adjust_scale(scale, coef(m1)["sl_"]) +
    coef(m1)["sl_:tod_end_night"]

# speed
speed_day <- shape * scale_adj_day
speed_night <- shape * scale_adj_night

speed_day <- shape_adj_day * scale
speed_night <- shape_adj_night * scale

speed_day <- shape_adj_day * scale_adj_day
speed_night <- shape_adj_night * scale_adj_night

scale
shape

shape_adj_day
shape_adj_night

x <- seq(1, 500, 1)
plot(x, dgamma(x, shape = shape_adj_night, scale = scale_adj_night), type = "l")
lines(x, dgamma(x, shape = shape_adj_day, scale = scale_adj_day), type = "l")

# Bootstrap everything
mod_data <- stps %>% amt::random_steps(n = 9) %>%
    amt::extract_covariates(ras, where = "end") %>%
    amt::time_of_day(include.crepuscule = FALSE) %>%
    mutate(log_sl_ = log(sl_), cos_ta_ = cos(as_rad(ta_)))

strata <- unique(mod_data$step_id_)
n <- length(strata)

bt <- replicate(1000, {
    m_boot <- mod_data[mod_data$step_id_ %in% sample(strata, n, TRUE), ] %>%
        amt::fit_clogit(case_ ~ Category + log_sl_ + sl_ + Category:tod_end_ +
                            sl_:tod_end_ + log_sl_:tod_end_ + strata(step_id_))
    
    scale_adj_day <- amt::adjust_scale(shape, coef(m_boot)["sl_"])
    scale_adj_night <- amt::adjust_scale(shape, coef(m_boot)["sl_"]) +
        coef(m_boot)["sl_:tod_end_night"]
    
    shape_adj_day <- amt::adjust_shape(shape, coef(m_boot)["log_sl_"])
    
    shape_adj_night <- amt::adjust_shape(shape, coef(m_boot)["log_sl_"]) +
        coef(m_boot)["log_sl_:tod_end_night"]
    
    ## speed
    c(shape_adj_day * scale_adj_day,
      shape_adj_night * scale_adj_night)
})

bt2 <- data_frame(
    rep = 1:ncol(bt),
    day = bt[1, ],
    night = bt[2, ]
) %>% gather(key, val, -rep)


bt2 %>% group_by(key) %>% summarise(lq = quantile(val, 0.025),
                                    me = median(val),
                                    mean = mean(val),
                                    uq = quantile(val, 0.975))

# m/min
bt2 %>% group_by(key) %>% summarise(lq = quantile(val, 0.025) / 10,
                                    me = median(val) / 10,
                                    mean = mean(val) / 10,
                                    uq = quantile(val, 0.975) / 10)

## Simulate ud
wet_c <- crop(Category, amt::bbox(dat, spatial = TRUE, buff = 1e3))

mk <- amt::movement_kernel(scale, shape_adj_day, wet_c)
hk <- amt::habitat_kernel(list(wet = coef(m1)["Category"]), wet_c)


system.time(ssud_day <- amt::simulate_ud(
    mk, hk,
    as.numeric(stps[1, c("x1_", "y1_")]),
    n = 1e7))
plot(ssud_day)

system.time(tud_day <- amt::simulate_tud(mk, hk, as.numeric(stps[1501, c("x1_", "y1_")]), n = 72, n_rep = 5e3))
plot(tud_day)


# night
mk <- amt::movement_kernel(scale, shape_adj_night, wet_c)
hk <- amt::habitat_kernel(list(wet = coef(m1)["wet"] +
                                   coef(m1)["wet:tod_end_night"]), wet_c)

system.time(ssud_night <- amt::simulate_ud(
    mk, hk, as.numeric(stps[1, c("x1_", "y1_")]), n = 1e7))
plot(ssud_night)

system.time(tud_night <- amt::simulate_tud(mk, hk, as.numeric(stps[1501, c("x1_", "y1_")]), n = 72, n_rep = 5e3))
plot(tud_day)
plot(tud1 <- crop(tud_day, extent(c(1778000, 1782000, 2412000, 2415000))))
plot(tud2 <- crop(tud_night, extent(c(1778000, 1782000, 2412000, 2415000))))

pllog <- list(
    geom_raster(),
    coord_equal(),
    scale_fill_continuous(low = "white", high = "red",
                          tran = "log10", na.value = "white"),
    scale_y_continuous(expand = c(0, 0)),
    scale_x_continuous(expand = c(0, 0)),
    theme_light(),
    theme(legend.position = "none"))

pl <- list(
    geom_raster(),
    coord_equal(),
    scale_fill_continuous(low = "white", high = "red", na.value = "white"),
    scale_y_continuous(expand = c(0, 0)),
    scale_x_continuous(expand = c(0, 0)),
    theme_light(),
    theme(legend.position = "none"))

r1 <- data.frame(rasterToPoints(mk))
p1 <- ggplot(r1, aes(x, y, fill = d)) + pllog + ggtitle("Movement kernel (night)")

r2 <- data.frame(rasterToPoints(hk))
p2 <- ggplot(r2, aes(x, y, fill = layer)) + pl + ggtitle("Habitat kernel (night)")


r1 <- data.frame(rasterToPoints(tud1))
p3 <- ggplot(r1, aes(x, y, fill = layer)) + pllog + ggtitle("Transient UD (day)")

r2 <- data.frame(rasterToPoints(tud2))
p4 <- ggplot(r2, aes(x, y, fill = layer)) + pllog + ggtitle("Transient UD (night)")


r1 <- data.frame(rasterToPoints(ssud_day))
p5 <- ggplot(r1, aes(x, y, fill = layer)) + pl + ggtitle("Steady state UD (day)")

r2 <- data.frame(rasterToPoints(ssud_night))
p6 <- ggplot(r2, aes(x, y, fill = layer)) + pl + ggtitle("Steady state UD (night)")

cowplot::plot_grid(p1, p2, p3, p5, p4, p6, ncol = 2, labels = "AUTO")
#########
    
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
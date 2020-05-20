# Attempt to use mixed effect selection model. 

setwd("D:/1_PhD/EgoAllo/EgoAlloHabitat")
library(raster)
library(lubridate)
library(rgdal)
library(amt)
library(tidyverse)
library(cowplot)

set.seed(106)

dat <- read.csv("Data/Atlas2018egoallo.csv")%>%
    filter(!is.na(`medianEast`)) %>%
    select(id ="Bird", sex = "Sex.x", mass = "Mass.25.07", x = "medianEast", y = "medianNort",
           t = "RealTime", strategy = "Strat")
dat$t <- as.POSIXct(dat$t)
dat <- dat[order(dat$id,dat$t),]
dat_all <- dat %>% 
    nest(-c(id,sex,mass,strategy)) # from purrr this nest command can be used to nest data into a list column

dat_all <- dat_all %>%
    mutate(trk = lapply(data, function(d){
        amt::make_track(d, x, y, t, crs = sp::CRS("+init=epsg:27700"))}))


habitat <- readOGR(dsn = "Data", layer = "HabitatLayer2")
habitat$Category<- as.integer(habitat$Category)
r <- raster(habitat,resolution=c(5,5))
ras <- rasterize(habitat,r,"Category", background=0)
ras[is.na(ras[])]<- 0
names(ras) <- "Category"


dat_all %>%
    mutate(sr = lapply(trk, summarize_sampling_rate)) %>%
    select(id, sr) %>% 
    unnest_legacy()

dat1 <- dat_all %>% mutate(dat_clean = map(trk, ~ {
    .x %>% track_resample(rate = minutes(5), tolerance = seconds(120))
}))

dat_ssf <- dat1 %>% 
    mutate(stps = map(dat_clean, ~ .x %>% 
                          steps_by_burst() %>% 
                          random_steps() %>% 
                          extract_covariates(ras))) %>% 
    select(id, stps) %>% unnest_legacy() %>% 
    mutate(
        y = as.numeric(case_),
        id = as.numeric(factor(id)), 
        step_id = paste0(id, step_id_, sep = "-"))
dat_ssf


library(glmmTMB)
TMBStruc <- glmmTMB::glmmTMB(y ~ -1 + Category +
                                 (1|step_id) +(0 + Category | id),
                             family=poisson, data =dat_ssf, doFit=FALSE) 
TMBStruc$parameters$theta[1] <- log(1e3) 
TMBStruc$mapArg <- list(theta=factor(c(NA,1)))

glmm.TMB.random <- glmmTMB:::fitTMB(TMBStruc)
summary(glmm.TMB.random)
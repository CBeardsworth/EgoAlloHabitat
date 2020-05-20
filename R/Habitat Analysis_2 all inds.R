# attempt to run a multi-individual iSSA
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

#EITHER create a table with all the ssfs stored within the object
# m1 <- dat_all %>%
#     mutate(ssf = lapply(trk, function(x){ #create ssf variable in m1 using dat_all
#         x %>% amt::track_resample(rate = minutes(5), #resample as differences in sampling rate between individuals
#                                   tolerance = minutes(2)) %>%
#             amt::filter_min_n_burst(min_n = 3) %>% # get minimum length of burst
#             amt::steps_by_burst() %>% #changes the data from a point representation into a step representation (step length, turn angles)
#             amt::random_steps() %>% #default is n=10
#             amt::extract_covariates(ras, where = "end") %>% #extract (environmental) covariates at the endpoint of each step
#             mutate(landuse_end = factor(Category)) %>% #create landuse variable
#             amt::fit_issf(case_ ~ landuse_end:sex + landuse_end:strategy + strata(step_id_)) #run ssf for each track
#     }))


#OR run model then plot
m1 <- dat_all %>%
    mutate(steps = lapply(trk, function(x){ #create ssf variable in m1 using dat_all
        x %>% amt::track_resample(rate = minutes(5), #resample as differences in sampling rate between individuals
                                  tolerance = minutes(1)) %>%
            amt::filter_min_n_burst(min_n = 3) %>% # get minimum length of burst
            amt::steps_by_burst() %>% #changes the data from a point representation into a step representation (step length, turn angles)
            amt::random_steps() %>% #default is n=10
            amt::extract_covariates(ras, where = "end") %>% #extract (environmental) covariates at the endpoint of each step
            mutate(landuse_end = factor(Category))})) #create landuse variable

head(m1$steps,n=1)
#######
m2 <- m1 %>% mutate(fit = map(steps, ~ amt::fit_issf(., case_ ~ landuse_end +
                                                         strata(step_id_))))

d2 <- m2 %>% mutate(coef = map(fit, ~ broom::tidy(.x$model))) %>%
    select(id, sex, strategy, coef) %>% unnest_legacy %>%
    mutate(id = factor(id)) %>% group_by(term) %>%
    summarize(
        mean = mean(estimate),
        ymin = mean - 1.96 * sd(estimate),
        ymax = mean + 1.96 * sd(estimate)
    )

d2$x <- 1:nrow(d2)


p1 <- m2 %>% mutate(coef = map(fit, ~ broom::tidy(.x$model))) %>%
    select(id, strategy,sex, coef) %>% 
    unnest_legacy %>%
    mutate(id = factor(id)) %>%
    ggplot(., aes(x = term, y = estimate, group = id, col = sex, pch = strategy)) +
    geom_rect(mapping = aes(xmin = x - .4, xmax = x + .4, ymin = ymin,
                            ymax = ymax), data = d2, inherit.aes = FALSE,
              fill = "grey90") +
    geom_segment(mapping = aes(x = x - .4, xend = x + .4,
                               y = mean, yend = mean), data = d2, inherit.aes = FALSE,
                 size = 1) +
    geom_pointrange(aes(ymin = conf.low, ymax = conf.high),
                    position = position_dodge(width = 0.7), size = 0.8) +
    geom_hline(yintercept = 0, lty = 2) +
    labs(x = "Habitat", y = "Relative Selection Strength") +
    theme_light()+
    scale_x_discrete(labels = c("Urban", "Woodland"))

p1

# attempt to run a multi-individual iSSA
setwd("D:/1_PhD/EgoAllo/EgoAlloHabitat")
library(raster)
library(lubridate)
library(rgdal)
library(amt)
library(tidyverse)
library(cowplot)

set.seed(106)
habitat <- readOGR(dsn = "Data", layer = "HabitatLayer2")
habitat$Category<- as.integer(habitat$Category)
r <- raster(habitat,resolution=c(5,5))
ras <- rasterize(habitat,r,"Category", background=0)
ras[is.na(ras[])]<- 0
names(ras) <- "Category"

dat <- read.csv("Data/Atlas2018egoallo.csv")%>%
    filter(!is.na(`medianEast`) & Day==T) %>%
    select(id ="Bird", sex = "Sex.x", mass = "Mass.25.07", x = "medianEast", y = "medianNort",
           t = "RealTime", strategy = "Strat", score = "Diff")#%>%    filter(!id %in% c(1075,1095))
dat$t <- as.POSIXct(dat$t)
dat <- dat[dat$t > as.Date("2018-08-17"),]
dat <- dat[order(dat$id,dat$t),]
dat_all <- dat %>% 
    nest(-c(id,sex,strategy, mass,score)) # from purrr this nest command can be used to nest data into a list column

dat_all <- dat_all %>%
    mutate(trk = lapply(data, function(d){
        amt::make_track(d, x, y, t, crs = sp::CRS("+init=epsg:27700"))}))

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
#             amt::fit_issf(case_ ~ landuse_end + strata(step_id_)) #run ssf for each track
#     }))

##### Selection Model #####
#OR run model then plot
m1 <- dat_all %>%
    mutate(steps = lapply(trk, function(x){ #create ssf variable in m1 using dat_all
        x %>% amt::track_resample(rate = minutes(5), #resample as differences in sampling rate between individuals
                                  tolerance = minutes(1)) %>%
            amt::filter_min_n_burst(min_n = 3) %>% # get minimum length of burst (30min)
            amt::steps_by_burst() %>% #changes the data from a point representation into a step representation (step length, turn angles)
            amt::random_steps(n_control = 10) %>% #default is n=10
            amt::extract_covariates(ras, where = "end") %>% #extract (environmental) covariates at the endpoint of each step
            mutate(habitat = factor(Category, levels=c(2,0,1)),
                   log_sl = log(sl_))
        
        })) #create landuse variable

#######

m1 <- m1 %>% mutate(fit = map(steps, ~ amt::fit_issf(., case_ ~ habitat + sl_ +
                                                         strata(step_id_))))

order <- dat_all[order(dat_all$sex,dat_all$strategy),1:3]

inds <- m1 %>% mutate(coef = map(fit, ~ broom::tidy(.x$model))) %>%
    select(id, sex, strategy, coef) %>% 
    unnest_legacy %>%
    mutate(id = factor(id)) %>% 
    group_by(sex, strategy, term)

groups <- m1 %>% mutate(coef = map(fit, ~ broom::tidy(.x$model))) %>%
    select(id, sex, strategy, coef) %>% 
    unnest_legacy %>%
    mutate(id = factor(id)) %>% 
    group_by(sex, strategy, term)%>%
    summarize(
        mean = mean(estimate),
        ymin = mean - 1.96 * sd(estimate),
        ymax = mean + 1.96 * sd(estimate)
    )

d2 <- m1 %>% mutate(coef = map(fit, ~ broom::tidy(.x$model))) %>%
    select(id, sex, strategy, coef) %>% 
    unnest_legacy %>%
    mutate(id = factor(id)) %>% 
    group_by(term) %>%
    summarize(
        mean = mean(estimate),
        ymin = mean - 1.96 * sd(estimate),
        ymax = mean + 1.96 * sd(estimate)
        )

#m <- glm(data=d2[d2$term=="landuse_end0",], estimate ~ strategy + mass)
#summary(m)


d2$x <- 1:nrow(d2)
groups$x <- 1:nrow(groups)
inds <- inds[inds$term!="sl_",]
inds$group <- paste0(inds$sex, "-",inds$strategy)

ggplot(data= inds, aes(x=term, y=estimate, col=group, shape = sex))+
    geom_boxplot()+
    geom_pointrange(aes(ymin = conf.low, ymax = conf.high),
                    position = position_jitterdodge(0.2), size = 0.4, alpha=0.5) +
    theme_classic()+
    scale_color_manual(values=c("royalblue", "gold3", "royalblue","gold3"))+
    
    geom_hline(yintercept = 0, lty = 2) +
    labs(x = "Habitat", y = "iSSA Estimate (Reference = Woodland)") +
    scale_x_discrete(labels = c("Open", "Urban"))#+
    #guides(col=F,shape=F)

ggsave("iSSA_estimates.png", units="cm", width=16, height=16, dpi = 600)


p1 <- m1 %>% mutate(coef = map(fit, ~ broom::tidy(.x$model))) %>%
    select(id, strategy,sex,coef) %>% 
    unnest_legacy %>%
    mutate(id = factor(id,levels=order$id)) %>%
    filter(term!="sl_")%>%
    
    ggplot(., aes(x = term, y = estimate, group = id, col = strategy, pch = sex)) +
    geom_rect(mapping = aes(xmin = x - .4, xmax = x + .4, ymin = ymin,
                            ymax = ymax), data = groups, inherit.aes = FALSE,
              fill = "grey90") +
    geom_segment(mapping = aes(x = x - .4, xend = x + .4,
                               y = mean, yend = mean), data = groups, inherit.aes = FALSE,
                 size = 1) +
    geom_pointrange(aes(ymin = conf.low, ymax = conf.high),
                    position = position_dodge(width = 0.7), size = 0.8) +
    geom_hline(yintercept = 0, lty = 2) +
    labs(x = "Habitat", y = "iSSA Estimate (Reference = Woodland)") +
    theme_light()+
    scale_x_discrete(labels = c("Open", "Urban"))

p1


#### Movement model ####

m2 <- dat_all %>%
    mutate(steps = lapply(trk, function(x){ #create ssf variable in m1 using dat_all
        x %>% amt::track_resample(rate = minutes(5), #resample as differences in sampling rate between individuals
                                  tolerance = minutes(1)) %>%
            amt::filter_min_n_burst(min_n = 6) %>% # get minimum length of burst (30min)
            amt::steps_by_burst() %>% #changes the data from a point representation into a step representation (step length, turn angles)
            amt::random_steps(n_control = 40) %>% #default is n=10
            amt::extract_covariates(ras, where = "start") %>% #extract (environmental) covariates at the startpoint of each step and look at interactions with step length/turning angle to determine movement in the habitat
            mutate(landuse = factor(Category, levels=c(0,1,2)),
                   log_sl=log(sl_))
        
    })) #create landuse variable

#head(m1)

#######
m2 <- m2 %>% mutate(fit = map(steps, ~ amt::fit_issf(., case_ ~ landuse:ta_ + landuse:log_sl +
                                                         strata(step_id_))))


shape <- amt::sl_shape(m2)
scale <- amt::sl_scale(m2)

speed_day <- shape * scale

#order <- dat_all[order(dat_all$sex,dat_all$strategy),1:3]

d2 <- m2 %>% mutate(coef = map(fit, ~ broom::tidy(.x$model))) %>%
    select(id, sex, strategy, coef) %>% unnest_legacy %>%
    mutate(id = factor(id)) %>% 
    group_by(term) %>%
    summarize(
        mean = mean(estimate),
        ymin = mean - 1.96 * sd(estimate),
        ymax = mean + 1.96 * sd(estimate)
    )

d2$x <- 1:nrow(d2)


p2 <- m2 %>% mutate(coef = map(fit, ~ broom::tidy(.x$model))) %>%
    select(id, strategy,sex, coef) %>% 
    unnest_legacy %>%
    mutate(id = factor(id,levels=order$id)) %>%
    ggplot(., aes(x = term, y = estimate, group = id, col = strategy, pch = sex)) +
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
    theme_light()#+
    scale_x_discrete(labels = c("Urban", "Woodland"))

p2

# Attempt to use mixed effect selection model. 

setwd("D:/1_PhD/EgoAllo/EgoAlloHabitat")
library(raster)
library(lubridate)
library(rgdal)
library(amt)
library(tidyverse)
library(cowplot)

set.seed(106)

dat <- read.csv("Data/Atlas2018egoallo.csv") %>%
    filter(!is.na(`medianEast`) & Day==T) %>%
    select(id ="Bird", sex = "Sex.x", mass = "Mass.25.07", x = "medianEast", y = "medianNort",
           t = "RealTime", strategy = "Strat")
    #filter(id!=1075)%>%
    #mutate(sc.score= scale(score))
dat$t <- as.POSIXct(dat$t)
dat <- dat[dat$t > as.Date("2018-08-17"),]

#dat$day <- round(difftime(as.POSIXct(dat$t), as.POSIXct("2018-07-26 00:00:00", tz="GMT"), units="days"))
#dat <- dat[dat$day >= 30,]
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
                          filter_min_n_burst(min_n = 5)%>%
                          steps_by_burst() %>% 
                          random_steps(n=10) %>% 
                          extract_covariates(ras, where="end"))) %>% 
    select(id,sex,strategy, stps) %>% unnest_legacy() %>% 
    mutate(
        y = as.numeric(case_),
        id = as.numeric(factor(id)),
        habitat = factor(Category, levels=c(0,1,2), labels = c("open", "urban", "woodland")),
        step_id = paste0(id, "-", step_id_))

#check <- head(dat_ssf,50)


library(glmmTMB) #apparently this is a frequentist approximation
TMBStruc <- glmmTMB::glmmTMB(y ~ -1 + habitat + 
                                 (1|step_id) + # allows the intercept to vary for every true step + all random steps (stratum) for of each individual i.e. we fit stratum specific FIXED intercepts (in line after, these are fixed with a LARGE (1e3) variance) which makes inference of random-slope models possible with standard Bayesian techniques or with frequentist
                               (0 + habitat | id),
                             family=poisson, data =dat_ssf, doFit=FALSE) 
TMBStruc$parameters$theta[1] <- log(1e3) 
TMBStruc$mapArg <- list(theta=factor(c(NA,1)))

glmm.TMB.random <- glmmTMB:::fitTMB(TMBStruc)
summary(glmm.TMB.random)
confint(glmm.TMB.random)

#with sex in the model
TMBStruc2 <- glmmTMB::glmmTMB(y ~ -1 + habitat * sex + 
                                 (1|step_id) + (1|id),# +(0 + habitat | id),
                             family=poisson, data =dat_ssf, doFit=FALSE) 
TMBStruc2$parameters$theta[1] <- log(1e3) 
TMBStruc2$mapArg <- list(theta=factor(c(NA,1)))

glmm.TMB.random2 <- glmmTMB:::fitTMB(TMBStruc2)
summary(glmm.TMB.random2)
confint(glmm.TMB.random2)

#with strategy in the model
TMBStruc3 <- glmmTMB::glmmTMB(y ~ -1 + habitat *strategy + 
                                 (1|step_id) +(0 + habitat | id),
                             family=poisson, data =dat_ssf, doFit=FALSE) 
TMBStruc3$parameters$theta[1] <- log(1e3) 
TMBStruc3$mapArg <- list(theta=factor(c(NA,1)))

glmm.TMB.random3 <- glmmTMB:::fitTMB(TMBStruc3)
summary(glmm.TMB.random3)
confint(glmm.TMB.random3)

#with strategy in the model
TMBStruc4 <- glmmTMB::glmmTMB(y ~ -1 + habitat + habitat:strategy + habitat:sex +
                                  (1|step_id) +(0 + habitat | id),
                              family=poisson, data =dat_ssf, doFit=FALSE) 
TMBStruc4$parameters$theta[1] <- log(1e3) # set the value of the standad deviation of the first random effect (here (1|step_id))
TMBStruc4$mapArg <- list(theta = factor(c(NA,1)))

glmm.TMB.random4 <- glmmTMB:::fitTMB(TMBStruc4)
summary(glmm.TMB.random4)
confint(glmm.TMB.random4)


# 2StepCLogit other method
library(TwoStepCLogit)
r.Twostep <-  Ts.estim(formula = y ~ habitat:sex  + habitat:strategy + strata(step_id) + 
                           cluster(id), data = dat_ssf, random = ~ habitat,
                       all.m.1=F, D="UN(1)") 

r.Twostep$beta

r.Twostep$se

r.Twostep$D

#INLA
install.packages("INLA", repos=c(getOption("repos"), INLA="https://inla.r-inla-download.org/R/stable"), dep=TRUE)
library(INLA)
#' Set mean and precision for the priors of slope coefficients
mean.beta <- 0
prec.beta <- 1e-4  

#' The model formula for INLA, where we set the stratum-specific intercept variance to $10^6$ (or rather: the precision to $10^{-6}$), is given as follows:
formula.fixed <-  y ~  habitat * sex + habitat * strategy + 
  f(step_id, model="iid", hyper=list(theta=list(initial=log(1e-6),fixed=T))) 

#' Then fit the Poisson model
r.inla.fixed <- inla(formula.fixed, family ="Poisson", data=dat_ssf,
                     control.fixed = list(
                       mean = mean.beta,
                       prec = list(default = prec.beta)
                     )
)

#' The summary for the posterior distribution of the fixed effects:
r.inla.fixed$summary.fixed

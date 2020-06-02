# attempt to run a multi-individual iSSA
setwd("D:/1_PhD/EgoAllo/EgoAlloHabitat")

library(raster)
library(lubridate)
library(rgdal)
library(amt)
library(tidyverse)
library(cowplot)
library(beepr)
library(ggbeeswarm)
library(MuMIn)
library(lme4)
#### Maze Task ####

learning <- read.csv("mazeData.csv")

# Did pheasants learn the maze task (errors)

m1 <- glmer(NumErrors ~ TrialNumber + Sex + Group + (1|Bird), data=learning, family = poisson(link=log))
summary(m1)
output <- dredge(m1)
output
est.output<-model.avg(output, subset= delta < 2, revised.var = TRUE)
summary(est.output)

ggplot(data=learning)+
    geom_beeswarm(aes(x=TrialNumber, y=NumErrors),dodge.width = 0, cex=0.5, alpha=0.3)+
    geom_smooth(aes(x=TrialNumber, y=NumErrors), col="grey4", method="lm")+
    scale_x_continuous(name="Trial Number", limits=c(0.2,8.8), breaks=seq(1,8,1))+
    scale_y_continuous(name="Number of Errors", limits=c(0,45), breaks=seq(0,45,5))+
    guides(col=F)+
    theme_classic()

#ggsave("learning.png", units="cm", width=24, height=12, dpi = 600)

## Did birds vary in orientation strategy? ##
orient<- read.csv("mazeRotationResults.csv")

g1 <- glm(binStrat~Sex*Group, family=binomial(link="logit"),data=orient)
summary(g1)

dredge(g1)
#summary top model - treatment only
summary(glm(binStrat~Group, family=binomial(link="logit"),data=orient))

## Plots ##

ggplot(data=orient, aes(x=Group, fill=Strat))+
    geom_bar(position="dodge", size=0)+
    scale_fill_manual(values=c("royalblue","gold3"))+
    scale_y_continuous(expand=c(0,0), limits=c(0,50))+
    #facet_grid(cols=vars(Group))+
    xlab("Condition")+
    ylab("Number of Birds") + 
    theme_classic() +
    theme(legend.position="none")

ggsave("BarchartDiffsEgoAllo.jpeg", dpi=600, units="cm", width = 8, height = 10)

orient2 <- orient[orient$Group=="Experimental",]

ggplot(data=orient2, aes(x=Diff, fill=Strat))+
    geom_histogram(breaks=seq(-10,30,by=1), size=0)+
    scale_fill_manual(values=c("royalblue","gold3"))+
    scale_y_continuous(expand=c(0,0), limits=c(0,8), breaks=seq(0,8,1))+
    scale_x_continuous(expand=c(0,0), limits=c(-10,35), breaks=seq(-10,30,5))+
    xlab("Difference in Errors \n(Probe Trial - Final Training Trial)")+
    ylab("Number of Birds") + 
    theme_classic() +
    theme(legend.position="none")

ggsave("ExperimentalDiffsinErrors.jpeg", dpi=600, units="cm", width = 8, height = 10)

#### Habitat Selection ####
rm(list = ls(all.names = TRUE)) 
set.seed(106)

# Habitat
ras <- raster("habitat.grd")

#Bird Data

dat <- read.csv("Data/atlas2018-strategy.csv")

#t <- table(dat$id, dat$month) # look at how many fixes per bird in each month (for appendix)
#write.csv(t, "Number_Locs-bird-month.csv")

dat$t <- as.POSIXct(dat$t)
dat <- dat[dat$t > as.Date("2018-08-17") & dat$t < as.Date("2018-10-31"),]
dat <- dat[order(dat$id,dat$t),]
#table(dat$id, dat$month) #see updated fixes per bird. Sample size lowered but this is more accurate because they are no longer in the release pen and it is in the same 'season'.


dat_all <- dat %>% 
    nest(-c(id,sex,strategy)) # from purrr this nest command can be used to nest data into a list column

dat_all <- dat_all %>%
    mutate(trk = lapply(data, function(d){
        amt::make_track(d, x, y, t, crs = sp::CRS("+init=epsg:27700"))}))

dat_all %>%
    mutate(sr = lapply(trk, summarize_sampling_rate)) %>%
    select(id, sr) %>% 
    unnest_legacy()

#### Import bootstrapped model estimates ####

bootData <- read.csv("bootData.csv")

#### Bootstrap model estimates ####
#set.seed(123)
# bootData <- NULL
# 
# for(i in 1:1000){
#     
#     print(paste("iteration:", i)) #count iterations to give estimate of time left.
#     
#     ## Prepare Data ##
#     mData <- dat_all %>%
#         mutate(steps = lapply(trk, function(x){ #create ssf variable in m1 using dat_all
#             x %>% amt::track_resample(rate = minutes(5), #resample as differences in sampling rate between individuals
#                                       tolerance = minutes(1)) %>%
#                 amt::filter_min_n_burst(min_n = 3) %>% # get minimum length of burst (30min)
#                 amt::steps_by_burst() %>% #changes the data from a point representation into a step representation (step length, turn angles)
#                 amt::random_steps(n_control = 10) %>% #default is n=10
#                 amt::extract_covariates(ras, where = "end") %>% #extract (environmental) covariates at the endpoint of each step
#                 mutate(habitat = factor(Category, levels=c(2,0,1))) #create landuse variable
#             
#         })) 
#     
#     ## Run Model ##
#     tryCatch({
#         m1 <- mData %>% mutate(fit = map(steps, ~ amt::fit_issf(., case_ ~ habitat + 
#                                                                 strata(step_id_))))}, error=function(e){cat("ERROR:", conditionMessage(e), "\n")}) 
#     
#     ## Extract model coefficients ##
#     inds <- m1 %>% mutate(coef = map(fit, ~ broom::tidy(.x$model))) %>%
#         select(id, sex, strategy, coef) %>% 
#         unnest_legacy %>%
#         mutate(id = factor(id)) %>% 
#         group_by(sex, strategy, term)
#     
#     inds$boot <- i # assign iteration number to new column
#     bootData <- rbind(bootData, inds) # attach to dataset that includes coefs from all iterations
#     
#     if(i == 1000){
#         beep() # tells you when the loop is over :)
#     }}

##### Calculate Coefficient Averages ####

bootMean <- bootData %>%
    group_by(id, mass, sex, strategy, term) %>% #get mean per individual, keep columns for sex strategy and term (habitat type) in the dataframe
    summarize(mean= mean(estimate),
              low.ci= quantile(estimate, 0.025),
              high.ci= quantile(estimate, 0.975))

#### Model selection (Woodland vs Open habitat). Criteria = lowest AIC 
woodVopen <- bootMean[bootMean$term=="habitat0",]
options(na.action="na.fail") #change for dredge function to run


model1 <- glm(mean ~ sex + strategy, data=woodVopen)
summary(model1)

model2 <- glm(mean ~ sex, data=woodVopen)
summary(model2) 

model3 <- glm(mean ~ mass + strategy, data=woodVopen)
summary(model3) 

model4 <- glm(mean ~ mass,data=woodVopen)
summary(model4) 

model5 <- glm(mean ~ strategy, data=woodVopen)
summary(model5)

model6 <- glm(data=woodVopen, mean ~ 1)
summary(model6) 

output <- model.sel(model1,model2,model3, model4,model5,model6)
output
est.output<-model.avg(output, subset= delta < 2, revised.var = TRUE)
summary(est.output)

#### Model selection (Woodland vs Urban habitat) = lowest AIC
woodVurban <- bootMean[bootMean$term=="habitat1",]

model1 <- glm(mean ~ sex + strategy, data=woodVurban)
summary(model1) #AIC =  16.608

model2 <- glm(mean ~ sex, data=woodVurban)
summary(model2) #AIC = 14.778

model3 <- glm(mean ~ mass + strategy, data=woodVurban)
summary(model3) # AIC = 13.146

model4 <- glm(mean ~ mass,data=woodVurban)
summary(model4) #AIC = 11.526 ###Best model

model5 <- glm(mean ~ strategy, data=woodVurban)
summary(model5) #AIC = 15.316

model6 <- glm(data=woodVurban, mean ~ 1)
summary(model6) #AIC = 15.529

#model7 <- glm(mean ~ sex*mass, data=woodVurban) # since sex and mass are in the top models, ensure that a model together (including an interaction) will not influence the results. It does not as this is delta 3.98

output <- model.sel(model1,model2,model3, model4,model5,model6)
output
est.output<-model.avg(output, subset= delta < 2, revised.var = TRUE)
summary(est.output)

#### Plot Results ####

bootMean$group <- paste0(bootMean$sex, "-",bootMean$strategy) #make new variable for groupings in the plot

popMean <- bootData %>% # Get population averages for paper results section and to show on plot
    group_by(sex, strategy, term) %>%
    summarize(mean= mean(estimate),
              low.ci= quantile(estimate, 0.025),
              high.ci= quantile(estimate, 0.975))

popMean$group <- paste0(popMean$sex, "-",popMean$strategy)

levels(bootMean$sex) <- c("Female", "Male") #rename from F and M for plot
levels(bootMean$term) <- c("Open", "Urban") #rename for plot
levels(bootMean$strategy) <- c("Allocentric", "Other") #rename for plot


ggplot(data= bootMean, aes(x=mass, y=mean,col=strategy, shape=sex))+
    facet_grid(.~term) +
    geom_point(size=3, alpha=0.1)+
    geom_pointrange(aes(ymin = low.ci, ymax = high.ci),
                    size = 0.6, alpha=0.7) +
    scale_color_manual(values=c("royalblue", "gold3"))+
    theme_classic()+
    geom_hline(yintercept = 0, lty = 2) +
    labs(x = "Mass (g)", y = "iSSA Estimate") 

ggsave("iSSA_estimates.png", units="cm", width=16, height=14, dpi = 600)

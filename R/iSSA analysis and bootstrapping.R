
#### Habitat Selection #### --------------------------------------------------------------------------------

# Habitat
ras <- raster("Data/habitat.grd")

#Bird Data

dat <- read.csv("Data/atlas2018-strategy.csv")

#t <- table(dat$id, dat$month) # look at how many fixes per bird in each month (for appendix)
#write.csv(t, "Number_Locs-bird-month.csv")

dat$t <- as.POSIXct(dat$t)
dat <- dat[dat$t > as.Date("2018-08-17") & dat$t < as.Date("2018-10-31"),]
dat <- dat[order(dat$id,dat$t),]
#table(dat$id, dat$month) #see updated fixes per bird. Sample size lowered but this is more accurate because they are no longer in the release pen and it is in the same 'season'.

#remove within release pen movements
dat$t <- as.POSIXct(dat$t, tz="GMT")
str(dat)
table(dat$id)

#find mean values to set in RSS calcs
feeders <- read.csv("Data/FeederCoords2017_27700.csv")%>%
    mutate(x = coords.x1, y = coords.x2)

dat_1 <- amt::make_track(dat, x, y, t, crs = sp::CRS("+init=epsg:27700"))%>%
    track_resample(rate = minutes(5), tolerance = minutes(1)) %>% 
    steps_by_burst()%>% #changes the data from a point representation into a step representation (step length, turn angles)
    mutate(log_sl_ = log(sl_), cos_ta_ = cos(ta_), dist_feedr = NA)
for(j in 1:nrow(dat_1)){
    dat_1$dist_feedr[j] <- min(sqrt((dat_1$x2_[j]-feeders$x)^2 + (dat_1$y2_[j]-feeders$y)^2)) #calculate distance to closest feeder
}

m_dist_feedr <- median(dat_1$dist_feedr)
m_sl <- median(dat_1$sl_, na.rm=T)
m_cos_ta <- median(dat_1$cos_ta_, na.rm=T)
m_log_sl <- median(dat_1$log_sl_, na.rm=T)


# Investigate Habitat (using end of steps) and movement (using start of steps)-------------------------------------------

birds <- unique(dat$id)
#test models
all_aic <- NULL

for(i in birds){
    
    print(paste("Modelling bird: ", i))
    dat_1 <- dat %>% filter(id == i)
    dat_1 <- amt::make_track(dat_1, x, y, t, crs = sp::CRS("+init=epsg:27700"))
    
    summarize_sampling_rate(dat_1)
    m_data <- amt::track_resample(dat_1, rate = minutes(5), tolerance = minutes(1)) %>% #resample, shouldn't be necessary with this dataset as all are at 5 min intervals anyway
        time_of_day(include.crepuscule = FALSE)%>% #wrapper around maptools::sunriset and maptools::crepuscule to determine day or night and also dawn/dusk if include.crepuscule = T
        filter(tod_ == "day") %>%
        #filter_min_n_burst(min_n = 3) %>% #choose to keep only bursts (subsets of the track with constant sampling rate within the specified tolerance)
        add_count(burst_)%>%#filter_min_n_burst does not work so use different function to remove bursts with only 2 points. 
        filter(n>=3)%>%
        steps_by_burst()%>% #changes the data from a point representation into a step representation (step length, turn angles)
        amt::random_steps(n = 10) %>% #fit gamma distribution to step lengths and von mises distribution to turn angles using maximum-likelihood, use these distributions to generate 10 random steps to each observed step. 
        amt::extract_covariates(ras, where = "both") %>% #extract (environmental) covariates at the endpoint of each step
        amt::time_of_day(include.crepuscule = FALSE) %>%
        mutate(log_sl_ = log(sl_), cos_ta_ = cos(ta_), 
               Habitat_start= factor(Category_start,levels=c(2,0,1),labels = c("wood", "other", "other")),
               Habitat_end = factor(Category_end,levels=c(2,0,1), labels = c("wood", "other", "other")), # some birds do not enter urban habitat enough to run the model. Combine open and urban
               #labels = c("wood", "open", "urb")), #original classifications
               dist_feedr = NA, log_dist_feedr = NA)
    
    for(j in 1:nrow(m_data)){
        m_data$dist_feedr[j] <- min(sqrt((m_data$x2_[j]-feeders$x)^2 + (m_data$y2_[j]-feeders$y)^2)) #calculate distance to closest feeder
    }
    
    m_data$log_dist_feedr <- log(m_data$dist_feedr)
    
    strata <- unique(m_data$step_id_)
    n <- length(strata)   
    bt <- replicate(10, {
        
        mod_data <- m_data[m_data$step_id_ %in% sample(strata, n, TRUE), ] 
        
        m1<- mod_data %>%
            amt::fit_clogit(case_ ~ Habitat_end + log_sl_ + sl_ + cos_ta_ + strata(step_id_), model=T)%>%
            AIC()
        m2<- mod_data %>%
            amt::fit_clogit(case_ ~ Habitat_end + log_sl_ + sl_ + cos_ta_ + dist_feedr + strata(step_id_), model=T)%>%
            AIC()
        m3<- mod_data %>%
            amt::fit_clogit(case_ ~ Habitat_end + log_sl_ + sl_ + cos_ta_ + log_dist_feedr + strata(step_id_), model=T)%>%
            AIC()
        m4<- mod_data %>%
            amt::fit_clogit(case_ ~ Habitat_start * (log_sl_ + sl_ + cos_ta_) + strata(step_id_), model=T)%>%
            AIC()
        m5<- mod_data %>%
            amt::fit_clogit(case_ ~ Habitat_start * (log_sl_ + sl_ + cos_ta_) + dist_feedr + strata(step_id_), model=T)%>%
            AIC()
        m6<- mod_data %>%
            amt::fit_clogit(case_ ~ Habitat_start * (log_sl_ + sl_ + cos_ta_) + log_dist_feedr + strata(step_id_), model=T)%>%
            AIC()
        m7<- mod_data %>%
            amt::fit_clogit(case_ ~ Habitat_end + Habitat_start * (log_sl_ + sl_ + cos_ta_) + strata(step_id_), model=T)%>%
            AIC()
        m8<- mod_data %>%
            amt::fit_clogit(case_ ~ Habitat_end + Habitat_start * (log_sl_ + sl_ + cos_ta_) + dist_feedr + strata(step_id_), model=T)%>%
            AIC()
        m9<- mod_data %>%
            amt::fit_clogit(case_ ~ Habitat_end + Habitat_start * (log_sl_ + sl_ + cos_ta_) + log_dist_feedr + strata(step_id_), model=T)%>%
            AIC()
        
        aics <- c(m1,m2,m3,m4,m5,m6,m7,m8,m9)
        
    })
    
    boot_aic <- data_frame(rep = 1:ncol(bt),
                           m1 = bt[1,],
                           m2 = bt[2,],
                           m3 = bt[3,], 
                           m4 = bt[4,],
                           m5 = bt[5,], 
                           m6 = bt[6,],
                           m7 = bt[7,],
                           m8 = bt[8,],
                           m9 = bt[9,]) %>% 
        gather(key, val, -rep)
    
    
    boot_aic <- boot_aic %>% 
        group_by(key) %>% 
        summarise(lq = quantile(val, 0.025),
                  me = median(val),
                  mean = mean(val),
                  uq = quantile(val, 0.975),
                  n = length(val),
                  sd= sd(val),
                  se = sd(val)/sqrt(length(val)),
                  inv_se = 1/(se^2)) %>%
        mutate(AICc = mean-min(mean), min = ifelse(AICc==0,1,0))
    
    
    boot_aic$id <- i
    all_aic <- rbind(all_aic, boot_aic)
}

table(all_aic$key,all_aic$min)

# M8 is the best model, use this formulation for both selection and movement.

all_avail <- NULL
all_RSS <- NULL
all_boot <- NULL
for(i in birds){
    
    print(paste("Modelling bird: ", i))
    dat_1 <- dat %>% filter(id == i)
    dat_1 <- amt::make_track(dat_1, x, y, t, crs = sp::CRS("+init=epsg:27700"))
    
    
    
    summarize_sampling_rate(dat_1)
    m_data <- amt::track_resample(dat_1, rate = minutes(5), tolerance = minutes(1)) %>% #resample, shouldn't be necessary with this dataset as all are at 5 min intervals anyway
        time_of_day(include.crepuscule = FALSE)%>% #wrapper around maptools::sunriset and maptools::crepuscule to determine day or night and also dawn/dusk if include.crepuscule = T
        filter(tod_ == "day") %>%
        #filter_min_n_burst(min_n = 3) %>% #choose to keep only bursts (subsets of the track with constant sampling rate within the specified tolerance)
        add_count(burst_)%>%#filter_min_n_burst does not work so use different function to remove bursts with only 2 points. 
        filter(n>=3)%>%
        steps_by_burst()%>% #changes the data from a point representation into a step representation (step length, turn angles)
        amt::random_steps(n = 10) %>% #fit gamma distribution to step lengths and von mises distribution to turn angles using maximum-likelihood, use these distributions to generate 10 random steps to each observed step. 
        amt::extract_covariates(ras, where = "both") %>% #extract (environmental) covariates at the endpoint of each step
        amt::time_of_day(include.crepuscule = FALSE) %>%
        mutate(log_sl_ = log(sl_), cos_ta_ = cos(ta_), 
               Habitat_start= factor(Category_start,levels=c(2,0,1),labels = c("wood", "other", "other")),
               Habitat_end = factor(Category_end,levels=c(2,0,1), labels = c("wood", "other", "other")), # some birds do not enter urban habitat enough to run the model. Combine open and urban
               #labels = c("wood", "open", "urb")), #original classifications
               dist_feedr = NA, log_dist_feedr = NA)
    
    for(j in 1:nrow(m_data)){
        m_data$dist_feedr[j] <- min(sqrt((m_data$x2_[j]-feeders$x)^2 + (m_data$y2_[j]-feeders$y)^2)) #calculate distance to closest feeder
    }
    
    m_data$log_dist_feedr <- log(m_data$dist_feedr)
    
    #get availability for summary tables
    avail <- as.data.frame(table(m_data[m_data$case_==FALSE,]$Habitat_end))
    avail$Var1 <- paste("Avail_Habitat", avail$Var1, sep="_")
    avail <- pivot_wider(avail, names_from = Var1, values_from = "Freq")
    av <- data.frame(id=i, avail)
    all_avail <- rbind(all_avail, av)
    
    # log rss----------------------------------------------------------------------------------------------
    
    # data.frame of x1s
    x1 <- data.frame(Habitat_end = factor(c("other")),
                     Habitat_start = factor(c("wood")),
                     dist_feedr = m_dist_feedr,
                     sl_ = m_sl,
                     cos_ta_ = m_cos_ta, 
                     log_sl_ = m_log_sl) 
    
    
    # data.frame of x2
    x2 <- data.frame(Habitat_end = factor("wood", levels = levels(m_end$Habitat)),
                     Habitat_start = factor("wood", levels = levels(m_end$Habitat)),
                     dist_feedr = m_dist_feedr,
                     sl_ = m_sl,
                     cos_ta_ = m_cos_ta, 
                     log_sl_ = m_log_sl) 
    
    x3 <- data.frame(Habitat_end = factor(c("other")),
                     Habitat_start = factor(c("other")),
                     dist_feedr = m_dist_feedr,
                     sl_ = m_sl,
                     cos_ta_ = m_cos_ta, 
                     log_sl_ = m_log_sl) 
    
    
    # data.frame of x4
    x4 <- data.frame(Habitat_end = factor("wood", levels = levels(m_end$Habitat)),
                     Habitat_start = factor("other", levels = levels(m_end$Habitat)),
                     dist_feedr = m_dist_feedr,
                     sl_ = m_sl,
                     cos_ta_ = m_cos_ta, 
                     log_sl_ = m_log_sl) 
    
    #model_best <- m_end %>% amt::fit_issf(case_ ~ Habitat * (log_sl_ + sl_ + cos_ta_) + dist_feedr + strata(step_id_), model=TRUE)
    #logRSS <- log_rss(object = model_best, x1 = x1, x2 = x2, ci = "boot", n_boot=10) # could use log_rss from amt package but this currently does not give se of resulting distribution which we need to put in the glm as inverse variance to account for uncertainty.
    
    # Bootstrap
    mod_data <- m_data
    strata <- unique(mod_data$step_id_)
    n <- length(strata)    
    
    bt <- replicate(1000, {
        
        m_boot <- mod_data[mod_data$step_id_ %in% sample(strata, n, TRUE), ] %>% 
            amt::fit_clogit(case_ ~ Habitat_end + Habitat_start * (log_sl_ + sl_ + cos_ta_) + dist_feedr + strata(step_id_), model=T)
        
        uncenter <- sum(coef(m_boot$model) * m_boot$model$means, na.rm=T)
        x1_dummy <- x1
        x2_dummy <- x2
        x3_dummy <- x3
        x4_dummy <- x4
        
        x1_dummy$step_id_ = m_boot$model$model$`strata(step_id_)`[1] #predict function won't run without the step id column
        x2_dummy$step_id_ = m_boot$model$model$`strata(step_id_)`[1]
        x3_dummy$step_id_ = m_boot$model$model$`strata(step_id_)`[1] #predict function won't run without the step id column
        x4_dummy$step_id_ = m_boot$model$model$`strata(step_id_)`[1]
        
        p1 <- predict(m_boot$model, newdata=x1_dummy, type= "lp", reference= "sample", se.fit=T) #this comes from the survival package. see brian smith's explanation of how this is a useful way of calculating RSS https://bsmity13.github.io/log_rss/#proof_of_predict()_method
        p2 <- predict(m_boot$model, newdata=x2_dummy, type= "lp", reference= "sample", se.fit=T)
        p3 <- predict(m_boot$model, newdata=x2_dummy, type= "lp", reference= "sample", se.fit=T)
        p4 <- predict(m_boot$model, newdata=x2_dummy, type= "lp", reference= "sample", se.fit=T)
        
        y_x1 <- p1$fit + uncenter
        y_x2 <- p2$fit + uncenter
        y_x3 <- p3$fit + uncenter
        y_x4 <- p4$fit + uncenter
        
        log_rss1 <- unname(y_x1 - y_x2) #Calculate log_rss wood2other
        log_rss2 <- unname(y_x3 - y_x4) #Calculate log_rss other2other
        
        # movement ---------------
        
        shp_ref <- as.numeric(m_boot$sl_$params$shape + m_boot$model$coefficients["log_sl_"])
        scl_ref <- as.numeric((m_boot$sl_$params$scale^-1) - m_boot$model$coefficients["sl_"])
        mean_speed_ref <- shp_ref/scl_ref
        
        #other
        shp_other <- as.numeric(m_boot$sl_$params$shape + (m_boot$model$coefficients["log_sl_"]+ m_boot$model$coefficients["Habitat_startother:log_sl_"]))
        scl_other <- as.numeric((m_boot$sl_$params$scale^-1) - (m_boot$model$coefficients["sl_"]+ m_boot$model$coefficients["Habitat_startother:sl_"]))
        mean_speed_other <- shp_other/scl_other
        
        #beta coefs
        habitat_end <- m_boot$model$coefficients["Habitat_endother"]
        costa <- m_boot$model$coefficients["cos_ta_"]
        sl <- m_boot$model$coefficients["sl_"]
        logsl <- m_boot$model$coefficients["log_sl_"]
        distfeedr <-  m_boot$model$coefficients["dist_feedr"]
        habitatsl <-  m_boot$model$coefficients["sl_"] + m_boot$model$coefficients["Habitat_startother:sl_"]
        habitatlogsl <-  m_boot$model$coefficients["log_sl_"] + m_boot$model$coefficients["Habitat_startother:log_sl_"]
        habitatcosta <- m_boot$model$coefficients["cos_ta_"] + m_boot$model$coefficients["Habitat_startother:cos_ta_"]
        
        #output
        booty <- c(mean_speed_ref,mean_speed_other, habitat_end, costa, habitatcosta ,sl, habitatsl, logsl, habitatlogsl, distfeedr, log_rss1,log_rss2)
        
    })
    boot_RSS <- data_frame(rep = 1:ncol(bt),
                           mean_speed_wood = bt[1,],
                           mean_speed_other = bt[2,], 
                           habitat_end = bt[3,],
                           cos_ta_ = bt[4,], 
                           habitat_start_cos_ta = bt[5,],
                           sl_ = bt[6,],
                           habitat_start_sl_ = bt[7,], 
                           log_sl_ = bt[8,], 
                           habitat_log_sl_ = bt[9,],
                           dist_feedr = bt[10,], 
                           log_RSS_wood2other= bt[11,],
                           log_RSS_other2other= bt[12,]) %>% 
        gather(key, val, -rep)
    all_boot$id <- i
    all_boot <- rbind(all_boot, boot_RSS)
    
    boot_RSS_sum <- boot_RSS %>% 
        group_by(key) %>% 
        summarise(lq = quantile(val, 0.025),
                  me = median(val),
                  mean = mean(val),
                  uq = quantile(val, 0.975),
                  n = length(val),
                  sd= sd(val),
                  se = sd(val)/sqrt(length(val)),
                  inv_se = 1/(se^2))
    
    boot_RSS_sum$id <- i
    all_RSS <- rbind(all_RSS, boot_RSS_sum)
    
    
}

all_RSS <- merge(all_RSS, unique(dat[,c(1,2,7)]), by = "id") # add sex and strategy
all_boot$id <- rep(birds, each =12000)
all_boot <- merge(all_boot, unique(dat[,c(1,2,7)]), by = "id")
#all_RSS <- read.csv("habitatOrientation_coefs.csv")
#avail <- read.csv("habitatOrientation_avail.csv")

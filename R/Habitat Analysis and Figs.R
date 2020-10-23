library(MuMIn)
library(amt)
library(lme4)
library(extrafont)
require(gridExtra)
setwd("D:/OneDrive - NIOZ/99_PhD/1_Manuscripts/EgoAllo Habitat/Data")

all_RSS <- read.csv("habitatOrientation_coefs.csv")
all_RSS$strategy <- factor(all_RSS$strategy, levels= c("Egocentric", "Allocentric")) # change reference level to get allocentric estimates instead of mixed/ego birds. Clearer to explain. 
all_avail <- read.csv("habitatOrientation_avail.csv")

# selection  (RSS)

rss_data <- all_RSS[all_RSS$key=="log_RSS_wood2other",]
rss_data <- merge(rss_data, all_avail, by = "id")
rss_data$fig_strategy <- factor(rss_data$strategy, levels=c("Allocentric", "Egocentric")) #change levels around so in same order as previous figs.

hab <- ggplot(data= rss_data, aes(x=fig_strategy, y=mean, col=strategy)) +
    geom_boxplot(show.legend = F, outlier.shape = NA)+
    geom_pointrange(aes(ymin = lq, ymax = uq, alpha = inv_se), size=0.3,
                    position=position_jitterdodge(), show.legend = F) +
    scale_color_manual(values=c( "gold3","royalblue"))+
    theme_classic()+
    theme(text=element_text(family="Calibri"), axis.text.x = element_text(size=7))+
    scale_x_discrete(labels=c("Allocentric" = "Allocentric", "Egocentric" = "Mixed/Egocentric"))+
    geom_hline(yintercept = 0, lty = 2) +
    labs(tag="a",x = "Strategy", y = "log-RSS")

hab
#ggsave("Fig3a_log_RSS.tiff", dpi=600, width= 5, height = 10, units = "cm")
habitat_glm <- glm(mean ~ sex * strategy + Avail_Habitat_other, data = rss_data, weight = inv_se, na.action="na.fail")
habitat_output <- dredge(habitat_glm, fixed = c("Avail_Habitat_other"))
habitat_output
habitat_est.output <- model.avg(habitat_output, subset= delta < 2, revised.var = TRUE)
summary(habitat_est.output)

#remove outlier
rss_data2 <- rss_data[rss_data$mean<0,]

ggplot(data= rss_data2, aes(x=fig_strategy, y=mean, col=strategy)) +
    geom_boxplot(show.legend = F, outlier.shape = NA)+
    geom_pointrange(aes(ymin = lq, ymax = uq, alpha = inv_se), size=0.3,
                    position=position_jitterdodge(), show.legend = F) +
    scale_color_manual(values=c( "gold3","royalblue"))+
    theme_classic()+
    theme(text=element_text(family="Calibri"), axis.text.x = element_text(size=7))+
    scale_x_discrete(labels=c("Allocentric" = "Allocentric", "Egocentric" = "Mixed/Egocentric"))+
    geom_hline(yintercept = 0, lty = 2) +
    labs(x = "Strategy", y = "log-RSS")


habitat_glm <- glm(mean ~ sex * strategy + Avail_Habitat_other, data = rss_data2, weight = inv_se, na.action="na.fail")
habitat_output <- dredge(habitat_glm, fixed = c("Avail_Habitat_other"))
habitat_output
habitat_est.output <- model.avg(habitat_output, subset= delta < 2, revised.var = TRUE)
summary(habitat_est.output)


# speed
speed_data <- all_RSS[all_RSS$key%in% c("mean_speed_wood", "mean_speed_other"),]
speed_data$key <- factor(speed_data$key, levels=c("mean_speed_wood", "mean_speed_other"), labels = c("Wood", "Other"))


mov1 <- ggplot(data= speed_data, aes(x=sex, y=mean, col=key)) +
    geom_boxplot(outlier.shape = NA)+
    geom_pointrange(aes(ymin = lq, ymax = uq, alpha = inv_se),
                    size = 0.3, position=position_jitterdodge(), show.legend = F) +
    scale_color_manual(values=c("forestgreen", "purple"))+
    theme_classic()+
    theme(legend.position = "none",text=element_text(family="Calibri"))+
    scale_y_continuous(limits=c(0,60), expand= c(0,0))+
    guides(col = guide_legend(title="Habitat"))+
    labs(tag="b",x = "Sex", y = "Mean Displacement Distance (m/5min)")
mov1

#ggsave("Fig3b_speed.tiff", dpi=600, width= 5, height = 10, units = "cm")

speed_glm <- glm(mean~ sex * strategy * key, data=speed_data, weight= inv_se,na.action="na.fail")
speed_output <- dredge(speed_glm)
speed_output #only 1 model
speed_glm <- glm(mean~ sex * key, data=speed_data, weight= inv_se,na.action="na.fail")
summary(speed_glm)
#speed_est.output<-model.avg(speed_output, subset= delta <= 2, revised.var = TRUE)
#summary(speed_est.output)
#all birds move faster in open, males move faster than females

# directionality 

dir_data <- all_RSS[all_RSS$key%in% c("cos_ta_", "habitat_start_cos_ta"),]
dir_data$key <- factor(dir_data$key, levels=c("cos_ta_", "habitat_start_cos_ta"), labels = c("Wood", "Other"))


mov2 <- ggplot(data= dir_data, aes(x=sex, y=mean, col=key)) +
    geom_boxplot(outlier.shape = NA)+
    geom_pointrange(aes(ymin = lq, ymax = uq, alpha = inv_se),
                    size = 0.3, position=position_jitterdodge(), show.legend = F) +
    scale_color_manual(values=c("forestgreen", "purple"))+
    theme_classic()+
    theme(legend.position = "none",text=element_text(family="Calibri"))+
    guides(col = guide_legend(title="Habitat"))+
    geom_hline(yintercept = 0, lty = 2) +
    labs(tag="c",x = "Sex", y = expression(paste("iSSA cos(turning angle) ", beta ," coefficient")))
mov2

#ggsave("Fig3c_directionality.tiff", dpi=600, width= 5, height = 10, units = "cm")

dir_glm <- glm(mean~ strategy * key * sex, data=dir_data, weights = inv_se,na.action="na.fail")
dir_output <- dredge(dir_glm)
dir_output
dir_est.output<-model.avg(dir_output, subset= delta < 2, revised.var = TRUE)
summary(dir_est.output)

mve <- grid.arrange(hab,mov1,mov2, ncol=3)
#setwd("D:/OneDrive - NIOZ/99_PhD/1_Manuscripts/EgoAllo Habitat/Submitted files_EcologyLetters3/Figures")
ggsave("Fig3.tiff",mve, dpi=600, width= 15, height = 10, units = "cm")


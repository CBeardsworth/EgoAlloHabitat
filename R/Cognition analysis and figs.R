library(MuMIn)
library(lme4)
library(extrafont)
require(gridExtra)
library(ggbeeswarm)
# Cognition analysis # ----------------------------------------------------------------------------------
setwd("D:/OneDrive - NIOZ/99_PhD/1_Manuscripts/EgoAllo Habitat/Data")

learning <- read.csv("mazeData.csv")

# Did pheasants learn the maze task (errors)

m1 <- glmer(NumErrors ~ TrialNumber + Sex + Group + (1|Bird), data=learning, family = poisson(link=log), na.action="na.fail")
summary(m1)
output <- dredge(m1)
output
est.output<-model.avg(output, subset= delta < 2, revised.var = TRUE)
summary(est.output)

ggplot(data=learning)+
    geom_beeswarm(aes(x=TrialNumber,y=NumErrors),size=0.01, cex=0.1)+
    geom_smooth(aes(x=TrialNumber, y=NumErrors), col="grey4", method="lm")+
    scale_x_continuous(expand=c(0,0),name="Trial Number", limits=c(0,9), breaks=seq(1,8,1))+
    scale_y_continuous(name="Number of Errors", limits=c(-2,max(learning$NumErrors)),breaks=seq(0,45,5), expand=c(0,0))+
    theme(text=element_text(family="Calibri"))+
    guides(col=F)+
    theme_classic()

ggsave("learning.tiff", units="cm", width=15, height=10, dpi = 600)

## Did birds vary in orientation strategy? ##
orient<- read.csv("mazeRotationResults.csv")

g1 <- glm(binStrat~Sex*Group, family=binomial(link="logit"),data=orient, na.action="na.fail")
summary(g1)

dredge(g1)
#summary top model - treatment only
summary(glm(binStrat~Group, family=binomial(link="logit"),data=orient))

## Plots ##

or1 <- ggplot(data=orient, aes(x=Group, fill=Strat))+
    geom_bar(position="dodge", size=0)+
    scale_fill_manual(values=c("royalblue","gold3"))+
    scale_y_continuous(expand=c(0,0), limits=c(0,50))+
    #facet_grid(cols=vars(Group))+
    xlab("Condition \n ")+
    ylab("Number of Birds") + 
    theme_classic() +
    theme(legend.position="none",text=element_text(family="Calibri"))+
    labs(tag="a")
or1
#ggsave("BarchartDiffsEgoAllo.jpeg", dpi=600, units="cm", width = 8, height = 10)

orient2 <- orient[orient$Group=="Experimental",]

or2 <- ggplot(data=orient2, aes(x=Diff, fill=Strat))+
    geom_histogram(breaks=seq(-10,30,by=1), size=0)+
    scale_fill_manual(values=c("royalblue","gold3"))+
    scale_y_continuous(expand=c(0,0), limits=c(0,8), breaks=seq(0,8,1))+
    scale_x_continuous(expand=c(0,0), limits=c(-10,35), breaks=seq(-10,30,5))+
    xlab("Difference in Errors \n(Probe Trial - Final Training Trial)")+
    ylab("Number of Birds") + 
    theme_classic() +
    theme(legend.position="none",text=element_text(family="Calibri"))+
    labs(tag="b")
or2
or <- grid.arrange(or1,or2, ncol=2)
ggsave("Fig2_allo-egomix.tiff", or, dpi=600, width=15, height = 10, units = "cm")

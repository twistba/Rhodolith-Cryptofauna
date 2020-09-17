# --- Crypotofuana Rhodolith paper - NIWA ---
## Author:       Brenton Twist ##
## Date Created: 2020-06-04 


# This script contains code for Running analyses on data collected on cryptofauna from two rhodolith species in northen NZ
# This script examines differences between sites for the same rhodolith species  -
# looking at differences in species richness of cryptofuana

# --- Packages ---
library(tidyverse)
library(lme4)
library(lmerTest)
library(ggsignif)
library(ggthemes)



# --- Load data ---
#   The following code loads the datasets and creating subsets

wd<-getwd()# gets working directory 
rich<-read.csv(paste(wd,"/data/cryptofauna_richness.csv",sep=""))

rich.l<-rich[rich$species == 'Lithothamnion',] #subset for Lithothamnion
rich.l$site<-as.factor(as.character(rich.l$site)) #remove sites that are no longer present
rich.l$transect<-as.factor(as.character(rich.l$transect))#remove transects that are no longer present

rich.s<-rich[rich$species == 'Sporolithon',] #subset for Sporolithon
rich.s$site<-as.factor(as.character(rich.s$site)) #remove sites that are no longer present
rich.s$transect<-as.factor(as.character(rich.s$transect))#remove transects that are no longer present


## --- Analysing data ----

# -Cryptofauna species richness -

#ALL
all.m1<-lmer(richness~species + (1|site/transect),data=rich)
plot(all.m1)#  homogeneity of variance -not normally distributed - log transform
shapiro.test(residuals(all.m1)) # significant not normal distributed

all.m2<-lmer(log(richness+1)~species + (1|site/transect),data=rich)#singular fit try removing transect random variable
all.m2<-lmer(log(richness+1)~species + (1|site),data=rich)
shapiro.test(residuals(all.m2))# normally distributed
summary(all.m2,ddf = "Kenward-Roger")
anova(all.m2,ddf = "Kenward-Roger")

#motile
mot.m1<-lmer(richness_motile~species+ (1|site/transect),data=rich)
plot(mot.m1)#  homogeneity of variance -not normally distributed - log transform
shapiro.test(residuals(mot.m1)) # significant not normally distributed

mot.m2<-lmer(log(richness_motile+1)~species + (1|site/transect),data=rich)
shapiro.test(residuals(mot.m2))# normally distributed
summary(mot.m2,ddf = "Kenward-Roger")
anova(mot.m2,ddf = "Kenward-Roger")

#sessile
ses.m2<-lmer(log(richness_sessile+1)~species+ (1|site/transect),data=rich)
summary(ses.m2,ddf = "Kenward-Roger")
anova(ses.m2,ddf = "Kenward-Roger")



# -Lithothamnion- Cryptofauna species richness -
  # Use Log + 1 transformations as all previous analyses used this

#ALL
all.l.m<-lmer(log(richness+1)~site + (1|transect),data=rich.l)
shapiro.test(residuals(all.l.m))# normally distributed
summary(all.l.m,ddf = "Kenward-Roger")
anova(all.l.m,ddf = "Kenward-Roger")

#mobile
mot.l.m<-lmer(log(richness_motile+1)~site + (1|transect),data=rich.l)
shapiro.test(residuals(mot.l.m))# normally distributed
summary(mot.l.m,ddf = "Kenward-Roger")
anova(mot.l.m,ddf = "Kenward-Roger")

#sessile
ses.l.m<-lmer(log(richness_sessile+1)~site + (1|transect),data=rich.l)
shapiro.test(residuals(ses.l.m))# normally distributed
summary(ses.l.m,ddf = "Kenward-Roger")
anova(ses.l.m,ddf = "Kenward-Roger")


# -Sporolithon- Cryptofauna species richness -
 # Use Log + 1 transformations as all previous analyses used this

#ALL
all.s.m<-lmer(log(richness+1)~site + (1|transect),data=rich.s)
shapiro.test(residuals(all.s.m))# normally distributed
summary(all.s.m,ddf = "Kenward-Roger")
anova(all.s.m,ddf = "Kenward-Roger")

#mobile
mot.s.m<-lmer(log(richness_motile+1)~site + (1|transect),data=rich.s)
shapiro.test(residuals(mot.s.m))# normally distributed
summary(mot.s.m,ddf = "Kenward-Roger")
anova(mot.s.m,ddf = "Kenward-Roger")

#sessile
ses.s.m<-lmer(log(richness_sessile+1)~site + (1|transect),data=rich.s)
shapiro.test(residuals(ses.s.m))# normally distributed
summary(ses.s.m,ddf = "Kenward-Roger")
anova(ses.s.m,ddf = "Kenward-Roger")

## --- END Analysing data ----


## --- Plotting data ----

# SE function
se <- function(x, na.rm=FALSE){
  if(na.rm){
    se <- sd(x, na.rm=TRUE)/sqrt(length(na.omit(x)))
  } else {
    se <- sd(x)/sqrt(length(x))
  }
  return(se)
}


# -Cryptofauna species richness: Figure 2 -

#summarising data
grouping <- group_by(rich, species)
p.all <- summarise(grouping, value=as.factor("Total"),mean=mean(richness),
                      se=se(richness))
p.mot <- summarise(grouping,value=as.factor("Motile"), mean=mean(richness_motile),
                      se=se(richness_motile))
p.ses <- summarise(grouping,value=as.factor("Sessile"), mean=mean(richness_sessile),
                      se=se(richness_sessile))
plot_data<-rbind(p.all,p.mot,p.ses)

#plotting data with ggplot
p1<-ggplot(data = plot_data, 
           aes(x=value,
               y= mean, 
               ymin=mean-se, 
               ymax=mean+se,       
               fill=species))
p1<- p1 + geom_bar(position="dodge", stat = "identity",colour="black")
p1<- p1 + geom_signif(annotations = c("***", "***","**"),
                      y_position = c(14,9,6), xmin=c(0.8,1.8,2.8), xmax=c(1.2,2.2,3.2))
p1<- p1 + scale_fill_manual(values=c("grey", "white"))
p1<- p1 + geom_errorbar( position = position_dodge(width=0.9),aes(width=0.3), colour="black")
#p1<- p1 + scale_x_discrete(labels=xlabs)
p1<- p1 + ggthemes::theme_few()
p1<- p1 + labs(y="Richness (number of species)", x=NULL,fill="Species")
p1<- p1 +theme(legend.spacing.y = unit(2, "mm"), 
               panel.border = element_rect(colour = "black", fill=NA),
               #axis.text.x = element_text(vjust=-2.5),
               legend.background = element_blank(),
               legend.box.background = element_rect(colour = "white"))# change to black for border
p1 <- p1 + theme(legend.position = c(0.85, 0.85))
p1

ggsave(p1, width = 21, height = 29.7/2, units = "cm", filename = "./figures/Fig2_Richness_species.png")



# -Lithothamnion- Cryptofauna species richness: Figure 3 -

#summarising data
grouping.l <- group_by(rich.l, site)
p.all.l <- summarise(grouping.l, value=as.factor("Total"),mean=mean(richness),
                   se=se(richness))
p.mot.l <- summarise(grouping.l,value=as.factor("Motile"), mean=mean(richness_motile),
                   se=se(richness_motile))
p.ses.l <- summarise(grouping.l,value=as.factor("Sessile"), mean=mean(richness_sessile),
                   se=se(richness_sessile))
plot_data_lith<-rbind(p.all.l,p.mot.l,p.ses.l)
plot_data_lith$site <- factor(plot_data_lith$site, levels = c("TMR","MTA"))#re-order levels

#plotting data with ggplot
p2<-ggplot(data = plot_data_lith, 
           aes(x=value,
               y= mean, 
               ymin=mean-se, 
               ymax=mean+se,       
               fill=site))
p2<- p2 + geom_bar(position="dodge", stat = "identity",colour="black")
p2<- p2 + geom_signif(annotations = c("ns", "ns","ns"),
                      y_position = c(15,10,6), xmin=c(0.8,1.8,2.8), xmax=c(1.2,2.2,3.2))
#p2<- p2 + scale_fill_manual(values=c("grey", "white"))
p2<- p2 + scale_fill_grey(start = 0.2, end = 0.8)
p2<- p2 + geom_errorbar( position = position_dodge(width=0.9),aes(width=0.3), colour="black")
#p2<- p2 + scale_x_discrete(labels=xlabs)
p2<- p2 + ggthemes::theme_few()
p2<- p2 + labs(y="Richness (number of species)", x=NULL,fill="Site")
p2<- p2 +theme(legend.spacing.y = unit(2, "mm"), 
               panel.border = element_rect(colour = "black", fill=NA),
               #axis.text.x = element_text(vjust=-2.5),
               legend.background = element_blank(),
               legend.box.background = element_rect(colour = "white"))# change to black for border
p2 <- p2 + theme(legend.position = c(0.85, 0.85),legend.title=element_text(size=16),legend.text=element_text(size=14),
                 axis.text=element_text(size=14),axis.title=element_text(size=16))
p2

ggsave(p2, width = 21, height = 29.7/2, units = "cm", dpi= 1200,
       filename = "./figures/Fig3_Lithothamnion_sites_richness.png")



# -Sporolithon- Cryptofauna species richness:Figure 4 -

#summarising data
grouping.s <- group_by(rich.s, site)
p.all.s <- summarise(grouping.s, value=as.factor("Total"),mean=mean(richness),
                     se=se(richness))
p.mot.s <- summarise(grouping.s,value=as.factor("Motile"), mean=mean(richness_motile),
                     se=se(richness_motile))
p.ses.s <- summarise(grouping.s,value=as.factor("Sessile"), mean=mean(richness_sessile),
                     se=se(richness_sessile))
plot_data_spor<-rbind(p.all.s,p.mot.s,p.ses.s)

#plotting data with ggplot
p3<-ggplot(data = plot_data_spor, 
           aes(x=value,
               y= mean, 
               ymin=mean-se, 
               ymax=mean+se,       
               fill=site))
p3<- p3 + geom_bar(position="dodge", stat = "identity",colour="black")
p3<- p3 + geom_signif(annotations = c("**", "*","*"),
                      y_position = c(11,5.5,6.25), xmin=c(0.8,1.8,2.8), xmax=c(1.2,2.2,3.2))
#p3<- p3 + scale_fill_manual(values=c("grey", "white"))
p3<- p3 + scale_fill_grey(start = 0.4, end = 0.8)
p3<- p3 + geom_errorbar( position = position_dodge(width=0.9),aes(width=0.3), colour="black")
#p3<- p3 + scale_x_discrete(labels=xlabs)
p3<- p3 + ggthemes::theme_few()
p3<- p3 + labs(y="Richness (number of species)", x=NULL,fill="Site")
p3<- p3 +theme(legend.spacing.y = unit(2, "mm"), 
               panel.border = element_rect(colour = "black", fill=NA),
               #axis.text.x = element_text(vjust=-2.5),
               legend.background = element_blank(),
               legend.box.background = element_rect(colour = "white"))# change to black for border
p3 <- p3 + theme(legend.position = c(0.85, 0.85),legend.title=element_text(size=16),legend.text=element_text(size=14),
                 axis.text=element_text(size=14),axis.title=element_text(size=16))
p3

ggsave(p3, width = 21, height = 29.7/2, units = "cm",dpi= 1200,
       filename = "./figures/Fig4_Sporolithon_sites_richness.png")


# --- Crypotofuana Rhodolith paper - NIWA ---
## Author:       Brenton Twist ##
## Date Created: 2020-06-04 


# This script contains code for Running analyses on data collected on cryptofauna from two rhodolith species in northen NZ
# This script examines differences between sites for the same rhodolith species  -
# looking at differences in multivariate community compostion between species and sites

# --- Packages ---
library(vegan)
library(ggplot2)
library(ggthemes)
library(indicspecies)



# --- Load data ---
#   The following code loads the datasets and creating subsets

wd<-getwd()# gets working directory
species<-read.csv(paste(wd,"/data/cryptofauna_species_abundance.csv",sep=""))

sp<-subset(species,select = -c(id,site,transect,species,rhodolith_number,species_site,
                               species_site_combined))#remove metadata to get species abundance matrix

sp.sqr<-sqrt(sqrt(sp))# fourth root transforms data

phy.bray<-vegdist(sp.sqr, method="bray")# fourth root transformed to Bray-Curtis dissimilarity matrix



## --- Plotting data ----

#produce nMDS plot
set.seed(1234)
sol<- metaMDS(sp.sqr, trymax=500) 

#Obtain plot coordinates
MDS1<-sol$point[,1]# stores x nMDS coordinates
MDS2<-sol$point[,2]# stores y nMDS coordinates
NMDS<-data.frame(species=species$species,site=species$site,
                 species_site=paste(species$species,species$site,sep="-"),
                 MDS1,MDS2)# dataframe with x y coordinates
NMDS$species_site <- factor(NMDS$species_site, 
                            levels = c("Lithothamnion-TMR", "Lithothamnion-MTA", 
                                       "Sporolithon-KWB" ,  "Sporolithon-MTA" ))


#plotting data with ggplot
p1 <- ggplot(NMDS, aes(x=MDS1, y=MDS2)) 
p1 <- p1 + geom_point(aes(shape=species_site),size=3)
p1 <- p1 + scale_shape_manual(values=c(1,0,17,15) )
#p1 <- p1 + geom_text_repel(aes(label=site),size=3,force=0.5)
#p1 <- p1 + geom_path(data=df_ell, aes(x=NMDS1, y=NMDS2,linetype=species), size=1)
#p1 <- p1 +  scale_linetype_manual(values=c("solid", "solid","twodash", "dotted", "solid"),guide=FALSE)
p1 <- p1 + theme_few()
p1 <- p1+ labs(shape='Species-Site') 
p1<- p1 + expand_limits(y = c(1.9),x=c(-1.9))
p1 <- p1 + theme(axis.title=element_blank(),axis.text=element_blank(),axis.ticks=element_blank())
p1 <- p1 + annotate("text", x = 0.8, y = 1.8, label = paste("2D Stress: ",round(sol$stress,digits=2)),size=5)
p1 <- p1 +  theme(legend.spacing.y = unit(0, "mm"), 
                  panel.border = element_rect(colour = "black", fill=NA),
                  legend.background = element_blank(),
                  legend.box.background = element_rect(colour = "white"))
p1 <- p1 + theme(legend.position = c(0.15, 0.85),legend.title=element_text(size=16),
                 legend.text=element_text(size=14))
p1

ggsave(p1, width = 21, height = 29.7/2,  units = "cm",dpi= 1200, 
       filename = "./figures/Fig5_nMDS_cryptofauna.png")                      

## --- END Plotting data ----


## --- Analysing data ----

#PERMANOVA
mod.all<- adonis(sp.sqr~species$site*species$species,method="bray")
mod.all


#Post-hoc test for PERMANOVA


# Creates a dataframe comparing each site species
set.seed(1238)
phy_sp<-subset(species,select = -c(id,transect,rhodolith_number,species_site_combined))
site.lst<-as.character(unique(phy_sp$species_site))
P<-c()
R2<-c()
MeanSqs<-c()
Df<-c()
comparison<-c()
F.Model<-c()
j=1
for(i in 1:4){
  for(n in 1:4){
    if(n>i){
      combin<-droplevels(subset(phy_sp,species_site==site.lst[i]|species_site==site.lst[n]))
      combin.df<-subset(combin,select=-c(site,species,species_site))
      combin.df.sq<-sqrt(sqrt(combin.df))
      permanova<- adonis(combin.df.sq~as.factor(combin$species_site) )
      F.Model[j]<- permanova$aov.tab$F.Model[1] 
      R2[j]<- permanova$aov.tab$R2[1] 
      P[j]<-permanova$aov.tab$`Pr(>F)`[1]
      MeanSqs[j]<-permanova$aov.tab$MeanSqs[1]
      Df[j]<-permanova$aov.tab$Df[1]
      comparison[j]<-paste(site.lst[i],site.lst[n],sep=" - ")
      j=j+1
    }                
    
  }
} 
adjusted.P<-p.adjust(P,method = "bonferroni", n= length(P))
PERM.comp<-data.frame(comparison,Df,MeanSqs,F.Model,R2,P,adjusted.P)
PERM.comp 

#indicator species analysis
indval <- multipatt(sp.sqr,species$species_site_combined,
                    control = how(nperm=999), duleg=TRUE) 
summary(indval)


# Write out the summary info for the model
sink("./output/PERMANOVA_species.txt")
cat("PERMANOVA combined site and species:\n")
print(mod.all)
cat("\n----\n\n")

cat("Post Hoc comparisons:\n")
print(PERM.comp)
cat("\n----\n\n")

cat("Indicator analysis:\n")
print(summary(indval))
cat("\n----\n\n")
sink()


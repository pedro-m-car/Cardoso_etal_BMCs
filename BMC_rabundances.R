##Reading data and setting libraries:

bmc_rabun <- read.csv("/Users/cardospm/Downloads/BMC experiment 2019/BMC_16s/BMC_rabundances.csv", header=TRUE, stringsAsFactors=FALSE)
bmc_rabun_t0 <- bmc_rabun[which(bmc_rabun$Timepoint == "T0"),names(bmc_rabun) %in% c("X","Treatment","Temperature","Timepoint","Group","Inoculation","Zotu42","Zotu145","Zotu1")]
bmc_rabun_t1 <- bmc_rabun[which(bmc_rabun$Timepoint == "T1"),names(bmc_rabun) %in% c("X","Treatment","Temperature","Timepoint","Group","Inoculation","Zotu42","Zotu145","Zotu1")]
bmc_rabun_t2 <- bmc_rabun[which(bmc_rabun$Timepoint == "T2"),names(bmc_rabun) %in% c("X","Treatment","Temperature","Timepoint","Group","Inoculation","Zotu42","Zotu145","Zotu1")]
bmc_rabun_t3 <- bmc_rabun[which(bmc_rabun$Timepoint == "T3"),names(bmc_rabun) %in% c("X","Treatment","Temperature","Timepoint","Group","Inoculation","Zotu42","Zotu145","Zotu1")]
bmc_rabun_t3_heat <- bmc_rabun_t3[which(bmc_rabun_t3$Temperature == "H"),names(bmc_rabun_t3) %in% c("X","Treatment","Temperature","Timepoint","Group","Inoculation","Zotu42","Zotu145","Zotu1")]
bmc_rabun_t3_con <- bmc_rabun_t3[which(bmc_rabun_t3$Temperature == "C"),names(bmc_rabun_t3) %in% c("X","Treatment","Temperature","Timepoint","Group","Inoculation","Zotu42","Zotu145","Zotu1")]

require(lattice)
require(visreg)
require(ggplot2)
library(dplyr)
library("ggpubr")
library("rcompanion")

## Creating Boxplots with relative abundances

rabun_42 <- ggplot(bmc_rabun,aes(Timepoint,Zotu42,fill=Treatment))+geom_boxplot(outlier.shape = NA)+labs(x="Timepoint", y=expression(paste("Relative Abundance of ", italic("Halomonas"), " sp. (%)")))+ylim(0,0.75) + geom_point(position = position_jitterdodge(jitter.width = 0.1), shape= 21, size =3) + stat_summary(fun.y=mean, geom="point", shape=4, size=4, aes(group=Treatment),position=position_dodge(0.75))+scale_fill_manual(values = c("#339966","#FF9900","#3399FF","#990000"))

rabun_42 + theme(
  legend.title = element_text(size=18, face="bold"),
  legend.text = element_text(size=18),
  axis.title.x = element_text(size=20),
  axis.title.y = element_text(size=20),
  axis.text = element_text(size=18),
  panel.background = element_rect(fill = "white",
                                  colour = "white",
                                  size = 0.5, linetype = "solid"),
  panel.grid.minor = element_line(size = 0.25, linetype = 'solid',
                                  colour = "grey"),
  axis.line.x = element_line(size = 1, linetype = "solid", colour = "black"),
  axis.line.y = element_line(size = 1, linetype = "solid", colour = "black"),
)

rabun_145 <- ggplot(bmc_rabun,aes(Timepoint,Zotu145,fill=Treatment))+geom_boxplot(outlier.shape = NA)+labs(x="Timepoint", y=expression(paste("Relative Abundance of ", italic("Cobetia"), " sp. (%)")))+ylim(0,0.6) + geom_point(position = position_jitterdodge(jitter.width = 0.1), shape= 21, size =3) + stat_summary(fun.y=mean, geom="point", shape=4, size=4, aes(group=Treatment),position=position_dodge(0.75))+scale_fill_manual(values = c("#339966","#FF9900","#3399FF","#990000"))

rabun_145 + theme(
  legend.title = element_text(size=18, face="bold"),
  legend.text = element_text(size=18),
  axis.title.x = element_text(size=20),
  axis.title.y = element_text(size=20),
  axis.text = element_text(size=18),
  panel.background = element_rect(fill = "white",
                                  colour = "white",
                                  size = 0.5, linetype = "solid"),
  panel.grid.minor = element_line(size = 0.25, linetype = 'solid',
                                  colour = "grey"),
  axis.line.x = element_line(size = 1, linetype = "solid", colour = "black"),
  axis.line.y = element_line(size = 1, linetype = "solid", colour = "black"),
)

rabun_1 <- ggplot(bmc_rabun,aes(Timepoint,Zotu1,fill=Treatment))+geom_boxplot(outlier.shape = NA)+labs(x="Timepoint", y=expression(paste("Relative Abundance of ", italic("Pseudoalteromonas"), " sp. (%)")))+ylim(0,100) + geom_point(position = position_jitterdodge(jitter.width = 0.1), shape= 21, size =3) + stat_summary(fun.y=mean, geom="point", shape=4, size=4, aes(group=Treatment),position=position_dodge(0.75))+scale_fill_manual(values = c("#339966","#FF9900","#3399FF","#990000"))

rabun_1 + theme(
  legend.title = element_text(size=18, face="bold"),
  legend.text = element_text(size=18),
  axis.title.x = element_text(size=20),
  axis.title.y = element_text(size=20),
  axis.text = element_text(size=18),
  panel.background = element_rect(fill = "white",
                                  colour = "white",
                                  size = 0.5, linetype = "solid"),
  panel.grid.minor = element_line(size = 0.25, linetype = 'solid',
                                  colour = "grey"),
  axis.line.x = element_line(size = 1, linetype = "solid", colour = "black"),
  axis.line.y = element_line(size = 1, linetype = "solid", colour = "black"),
)

## Testing data normality:

ggqqplot(bmc_rabun$Zotu42)

ggqqplot(bmc_rabun$Zotu145)

shapiro.test(bmc_rabun$Zotu42)

shapiro.test(bmc_rabun$Zotu145)

## Testing Hypothesis:

###Zotu42

###T0

library(vegan)

betadisper_BMC_T0_zotu42<-betadisper(d=vegdist(bmc_rabun_t0$Zotu42, method="euclidean"), group=bmc_rabun_t0$Group, type="median") #before using PERMANOVA homoscedasticity of groups must be tested
betadisper_BMC_T0_zotu42

permutest(x=betadisper_BMC_T0_zotu42, pairwise=TRUE, permutations=9999, p.method="bonf") 


perm.betadisper_BMC_T0_zotu42<-adonis(bmc_rabun_t0$Zotu42~Temperature*Inoculation, data=bmc_rabun_t0, permutations=9999, method="euclidean", p.method="bonf") #PERMANOVA
perm.betadisper_BMC_T0_zotu42

library(RVAideMemoire)

pairwise.BMC_T0_zotu42_TEMP<-pairwise.perm.manova(resp=vegdist(bmc_rabun_t0$Zotu42, "euclidean"), fact=bmc_rabun_t0$Temperature, test=NULL, nperm=9999, p.method="bonf")
pairwise.BMC_T0_zotu42_TEMP

pairwise.BMC_T0_zotu42_TREAT<-pairwise.perm.manova(resp=vegdist(bmc_rabun_t0$Zotu42, "euclidean"), fact=bmc_rabun_t0$Inoculation, test=NULL, nperm=9999, p.method="bonf")
pairwise.BMC_T0_zotu42_TREAT

pairwise.BMC_T0_zotu42<-pairwise.perm.manova(resp=vegdist(bmc_rabun_t0$Zotu42, "euclidean"), fact=bmc_rabun_t0$Treatment, test=NULL, nperm=9999, p.method="bonf")
pairwise.BMC_T0_zotu42

###T1

library(vegan)

betadisper_BMC_T1_zotu42<-betadisper(d=vegdist(bmc_rabun_t1$Zotu42, method="euclidean"), group=bmc_rabun_t1$Group, type="median") #before using PERMANOVA homoscedasticity of groups must be tested
betadisper_BMC_T1_zotu42

permutest(x=betadisper_BMC_T1_zotu42, pairwise=TRUE, permutations=9999, p.method="bonf") 


perm.betadisper_BMC_T1_zotu42<-adonis(bmc_rabun_t1$Zotu42~Temperature*Inoculation, data=bmc_rabun_t1, permutations=9999, method="euclidean", p.method="bonf") #PERMANOVA
perm.betadisper_BMC_T1_zotu42

library(RVAideMemoire)

pairwise.BMC_T1_zotu42_TEMP<-pairwise.perm.manova(resp=vegdist(bmc_rabun_t1$Zotu42, "euclidean"), fact=bmc_rabun_t1$Temperature, test=NULL, nperm=9999, p.method="bonf")
pairwise.BMC_T1_zotu42_TEMP

pairwise.BMC_T1_zotu42_TREAT<-pairwise.perm.manova(resp=vegdist(bmc_rabun_t1$Zotu42, "euclidean"), fact=bmc_rabun_t1$Inoculation, test=NULL, nperm=9999, p.method="bonf")
pairwise.BMC_T1_zotu42_TREAT

pairwise.BMC_T1_zotu42<-pairwise.perm.manova(resp=vegdist(bmc_rabun_t1$Zotu42, "euclidean"), fact=bmc_rabun_t1$Treatment, test=NULL, nperm=9999, p.method="bonf")
pairwise.BMC_T1_zotu42

###T2

library(vegan)

betadisper_BMC_T2_zotu42<-betadisper(d=vegdist(bmc_rabun_t2$Zotu42, method="euclidean"), group=bmc_rabun_t2$Group, type="median") #before using PERMANOVA homoscedasticity of groups must be tested
betadisper_BMC_T2_zotu42

permutest(x=betadisper_BMC_T2_zotu42, pairwise=TRUE, permutations=9999, p.method="bonf") 


perm.betadisper_BMC_T2_zotu42<-adonis(bmc_rabun_t2$Zotu42~Temperature*Inoculation, data=bmc_rabun_t2, permutations=9999, method="euclidean", p.method="bonf") #PERMANOVA
perm.betadisper_BMC_T2_zotu42

library(RVAideMemoire)

pairwise.BMC_T2_zotu42_TEMP<-pairwise.perm.manova(resp=vegdist(bmc_rabun_t2$Zotu42, "euclidean"), fact=bmc_rabun_t2$Temperature, test=NULL, nperm=9999, p.method="bonf")
pairwise.BMC_T2_zotu42_TEMP

pairwise.BMC_T2_zotu42_TREAT<-pairwise.perm.manova(resp=vegdist(bmc_rabun_t2$Zotu42, "euclidean"), fact=bmc_rabun_t2$Inoculation, test=NULL, nperm=9999, p.method="bonf")
pairwise.BMC_T2_zotu42_TREAT

pairwise.BMC_T2_zotu42<-pairwise.perm.manova(resp=vegdist(bmc_rabun_t2$Zotu42, "euclidean"), fact=bmc_rabun_t2$Treatment, test=NULL, nperm=9999, p.method="bonf")
pairwise.BMC_T2_zotu42

##Samples not homoscedastic, need to test with non-parametric:
library(FSA)

wilcox.test(bmc_rabun_t2_heat$Zotu42 ~ sample_data(bmc_rabun_t2_heat)$Treatment)
wilcox.test(bmc_rabun_t2_con$Zotu42 ~ sample_data(bmc_rabun_t2_con)$Treatment)

kruskal.test(bmc_rabun_t2$Zotu42 ~ sample_data(bmc_rabun_t2)$Treatment)

pairwise.wilcox.test(bmc_rabun_t2$Zotu42, bmc_rabun_t2$Treatment, p.adjust.method = 'BH')

dunnTest(Zotu42 ~ Treatment,
         data = bmc_rabun_t2,
         method = "bh"
)

scheirerRayHare(Zotu42 ~ Temperature + Inoculation,
                data = bmc_rabun_t2)

###T3

library(vegan)

betadisper_BMC_T3_zotu42<-betadisper(d=vegdist(bmc_rabun_t3$Zotu42, method="euclidean"), group=bmc_rabun_t3$Group, type="median") #before using PERMANOVA homoscedasticity of groups must be tested
betadisper_BMC_T3_zotu42

permutest(x=betadisper_BMC_T3_zotu42, pairwise=TRUE, permutations=9999, p.method="bonf") 


perm.betadisper_BMC_T3_zotu42<-adonis(bmc_rabun_t3$Zotu42~Temperature*Inoculation, data=bmc_rabun_t3, permutations=9999, method="euclidean", p.method="bonf") #PERMANOVA
perm.betadisper_BMC_T3_zotu42

library(RVAideMemoire)

pairwise.BMC_T3_zotu42_TEMP<-pairwise.perm.manova(resp=vegdist(bmc_rabun_t3$Zotu42, "euclidean"), fact=bmc_rabun_t3$Temperature, test=NULL, nperm=9999, p.method="bonf")
pairwise.BMC_T3_zotu42_TEMP

pairwise.BMC_T3_zotu42_TREAT<-pairwise.perm.manova(resp=vegdist(bmc_rabun_t3$Zotu42, "euclidean"), fact=bmc_rabun_t3$Inoculation, test=NULL, nperm=9999, p.method="bonf")
pairwise.BMC_T3_zotu42_TREAT

pairwise.BMC_T3_zotu42<-pairwise.perm.manova(resp=vegdist(bmc_rabun_t3$Zotu42, "euclidean"), fact=bmc_rabun_t3$Treatment, test=NULL, nperm=9999, p.method="bonf")
pairwise.BMC_T3_zotu42

###Zotu145

###T0

library(vegan)

betadisper_BMC_T0_zotu145<-betadisper(d=vegdist(bmc_rabun_t0$Zotu145, method="euclidean"), group=bmc_rabun_t0$Group, type="median") #before using PERMANOVA homoscedasticity of groups must be tested
betadisper_BMC_T0_zotu145

permutest(x=betadisper_BMC_T0_zotu145, pairwise=TRUE, permutations=9999, p.method="bonf") 


perm.betadisper_BMC_T0_zotu145<-adonis(bmc_rabun_t0$Zotu145~Temperature*Inoculation, data=bmc_rabun_t0, permutations=9999, method="euclidean", p.method="bonf") #PERMANOVA
perm.betadisper_BMC_T0_zotu145

library(RVAideMemoire)

pairwise.BMC_T0_zotu145_TEMP<-pairwise.perm.manova(resp=vegdist(bmc_rabun_t0$Zotu145, "euclidean"), fact=bmc_rabun_t0$Temperature, test=NULL, nperm=9999, p.method="bonf")
pairwise.BMC_T0_zotu145_TEMP

pairwise.BMC_T0_zotu145_TREAT<-pairwise.perm.manova(resp=vegdist(bmc_rabun_t0$Zotu145, "euclidean"), fact=bmc_rabun_t0$Inoculation, test=NULL, nperm=9999, p.method="bonf")
pairwise.BMC_T0_zotu145_TREAT

pairwise.BMC_T0_zotu145<-pairwise.perm.manova(resp=vegdist(bmc_rabun_t0$Zotu145, "euclidean"), fact=bmc_rabun_t0$Treatment, test=NULL, nperm=9999, p.method="bonf")
pairwise.BMC_T0_zotu145

###T1

library(vegan)

betadisper_BMC_T1_zotu145<-betadisper(d=vegdist(bmc_rabun_t1$Zotu145, method="euclidean"), group=bmc_rabun_t1$Group, type="median") #before using PERMANOVA homoscedasticity of groups must be tested
betadisper_BMC_T1_zotu145

permutest(x=betadisper_BMC_T1_zotu145, pairwise=TRUE, permutations=9999, p.method="bonf") 


perm.betadisper_BMC_T1_zotu145<-adonis(bmc_rabun_t1$Zotu145~Temperature*Inoculation, data=bmc_rabun_t1, permutations=9999, method="euclidean", p.method="bonf") #PERMANOVA
perm.betadisper_BMC_T1_zotu145

library(RVAideMemoire)

pairwise.BMC_T1_zotu145_TEMP<-pairwise.perm.manova(resp=vegdist(bmc_rabun_t1$Zotu145, "euclidean"), fact=bmc_rabun_t1$Temperature, test=NULL, nperm=9999, p.method="bonf")
pairwise.BMC_T1_zotu145_TEMP

pairwise.BMC_T1_zotu145_TREAT<-pairwise.perm.manova(resp=vegdist(bmc_rabun_t1$Zotu145, "euclidean"), fact=bmc_rabun_t1$Inoculation, test=NULL, nperm=9999, p.method="bonf")
pairwise.BMC_T1_zotu145_TREAT

pairwise.BMC_T1_zotu145<-pairwise.perm.manova(resp=vegdist(bmc_rabun_t1$Zotu145, "euclidean"), fact=bmc_rabun_t1$Treatment, test=NULL, nperm=9999, p.method="bonf")
pairwise.BMC_T1_zotu145

###T2

library(vegan)

betadisper_BMC_T2_zotu145<-betadisper(d=vegdist(bmc_rabun_t2$Zotu145, method="euclidean"), group=bmc_rabun_t2$Group, type="median") #before using PERMANOVA homoscedasticity of groups must be tested
betadisper_BMC_T2_zotu145

permutest(x=betadisper_BMC_T2_zotu145, pairwise=TRUE, permutations=9999, p.method="bonf") 


perm.betadisper_BMC_T2_zotu145<-adonis(bmc_rabun_t2$Zotu145~Temperature*Inoculation, data=bmc_rabun_t2, permutations=9999, method="euclidean", p.method="bonf") #PERMANOVA
perm.betadisper_BMC_T2_zotu145

library(RVAideMemoire)

pairwise.BMC_T2_zotu145_TEMP<-pairwise.perm.manova(resp=vegdist(bmc_rabun_t2$Zotu145, "euclidean"), fact=bmc_rabun_t2$Temperature, test=NULL, nperm=9999, p.method="bonf")
pairwise.BMC_T2_zotu145_TEMP

pairwise.BMC_T2_zotu145_TREAT<-pairwise.perm.manova(resp=vegdist(bmc_rabun_t2$Zotu145, "euclidean"), fact=bmc_rabun_t2$Inoculation, test=NULL, nperm=9999, p.method="bonf")
pairwise.BMC_T2_zotu145_TREAT

pairwise.BMC_T2_zotu145<-pairwise.perm.manova(resp=vegdist(bmc_rabun_t2$Zotu145, "euclidean"), fact=bmc_rabun_t2$Treatment, test=NULL, nperm=9999, p.method="bonf")
pairwise.BMC_T2_zotu145

###T3

library(vegan)

betadisper_BMC_T3_zotu145<-betadisper(d=vegdist(bmc_rabun_t3$Zotu145, method="euclidean"), group=bmc_rabun_t3$Group, type="median") #before using PERMANOVA homoscedasticity of groups must be tested
betadisper_BMC_T3_zotu145

permutest(x=betadisper_BMC_T3_zotu145, pairwise=TRUE, permutations=9999, p.method="bonf") 


perm.betadisper_BMC_T3_zotu145<-adonis(bmc_rabun_t3$Zotu145~Temperature*Inoculation, data=bmc_rabun_t3, permutations=9999, method="euclidean", p.method="bonf") #PERMANOVA
perm.betadisper_BMC_T3_zotu145

library(RVAideMemoire)

pairwise.BMC_T3_zotu145_TEMP<-pairwise.perm.manova(resp=vegdist(bmc_rabun_t3$Zotu145, "euclidean"), fact=bmc_rabun_t3$Temperature, test=NULL, nperm=9999, p.method="bonf")
pairwise.BMC_T3_zotu145_TEMP

pairwise.BMC_T3_zotu145_TREAT<-pairwise.perm.manova(resp=vegdist(bmc_rabun_t3$Zotu145, "euclidean"), fact=bmc_rabun_t3$Inoculation, test=NULL, nperm=9999, p.method="bonf")
pairwise.BMC_T3_zotu145_TREAT

pairwise.BMC_T3_zotu145<-pairwise.perm.manova(resp=vegdist(bmc_rabun_t3$Zotu145, "euclidean"), fact=bmc_rabun_t3$Treatment, test=NULL, nperm=9999, p.method="bonf")
pairwise.BMC_T3_zotu145

##Samples not homoscedastic, need to test with non-parametric:
library(FSA)

wilcox.test(bmc_rabun_t3_heat$Zotu145 ~ sample_data(bmc_rabun_t3_heat)$Treatment)
wilcox.test(bmc_rabun_t3_con$Zotu145 ~ sample_data(bmc_rabun_t3_con)$Treatment)

kruskal.test(bmc_rabun_t3$Zotu145 ~ sample_data(bmc_rabun_t3)$Treatment)

pairwise.wilcox.test(bmc_rabun_t3$Zotu145, bmc_rabun_t3$Treatment, p.adjust.method = 'BH')

dunnTest(Zotu145 ~ Treatment,
         data = bmc_rabun_t3,
         method = "bh"
)

scheirerRayHare(Zotu145 ~ Temperature + Inoculation,
                data = bmc_rabun_t3)

###zotu1

###T0

library(vegan)

betadisper_BMC_T0_zotu1<-betadisper(d=vegdist(bmc_rabun_t0$Zotu1, method="euclidean"), group=bmc_rabun_t0$Group, type="median") #before using PERMANOVA homoscedasticity of groups must be tested
betadisper_BMC_T0_zotu1

permutest(x=betadisper_BMC_T0_zotu1, pairwise=TRUE, permutations=9999, p.method="bonf") 


perm.betadisper_BMC_T0_zotu1<-adonis(bmc_rabun_t0$Zotu1~Temperature*Inoculation, data=bmc_rabun_t0, permutations=9999, method="euclidean", p.method="bonf") #PERMANOVA
perm.betadisper_BMC_T0_zotu1

library(RVAideMemoire)

pairwise.BMC_T0_zotu1_TEMP<-pairwise.perm.manova(resp=vegdist(bmc_rabun_t0$Zotu1, "euclidean"), fact=bmc_rabun_t0$Temperature, test=NULL, nperm=9999, p.method="bonf")
pairwise.BMC_T0_zotu1_TEMP

pairwise.BMC_T0_zotu1_TREAT<-pairwise.perm.manova(resp=vegdist(bmc_rabun_t0$Zotu1, "euclidean"), fact=bmc_rabun_t0$Inoculation, test=NULL, nperm=9999, p.method="bonf")
pairwise.BMC_T0_zotu1_TREAT

pairwise.BMC_T0_zotu1<-pairwise.perm.manova(resp=vegdist(bmc_rabun_t0$Zotu1, "euclidean"), fact=bmc_rabun_t0$Treatment, test=NULL, nperm=9999, p.method="bonf")
pairwise.BMC_T0_zotu1

###T1

library(vegan)

betadisper_BMC_T1_zotu1<-betadisper(d=vegdist(bmc_rabun_t1$Zotu1, method="euclidean"), group=bmc_rabun_t1$Group, type="median") #before using PERMANOVA homoscedasticity of groups must be tested
betadisper_BMC_T1_zotu1

permutest(x=betadisper_BMC_T1_zotu1, pairwise=TRUE, permutations=9999, p.method="bonf") 


perm.betadisper_BMC_T1_zotu1<-adonis(bmc_rabun_t1$Zotu1~Temperature*Inoculation, data=bmc_rabun_t1, permutations=9999, method="euclidean", p.method="bonf") #PERMANOVA
perm.betadisper_BMC_T1_zotu1

library(RVAideMemoire)

pairwise.BMC_T1_zotu1_TEMP<-pairwise.perm.manova(resp=vegdist(bmc_rabun_t1$Zotu1, "euclidean"), fact=bmc_rabun_t1$Temperature, test=NULL, nperm=9999, p.method="bonf")
pairwise.BMC_T1_zotu1_TEMP

pairwise.BMC_T1_zotu1_TREAT<-pairwise.perm.manova(resp=vegdist(bmc_rabun_t1$Zotu1, "euclidean"), fact=bmc_rabun_t1$Inoculation, test=NULL, nperm=9999, p.method="bonf")
pairwise.BMC_T1_zotu1_TREAT

pairwise.BMC_T1_zotu1<-pairwise.perm.manova(resp=vegdist(bmc_rabun_t1$Zotu1, "euclidean"), fact=bmc_rabun_t1$Treatment, test=NULL, nperm=9999, p.method="bonf")
pairwise.BMC_T1_zotu1

###T2

library(vegan)

betadisper_BMC_T2_zotu1<-betadisper(d=vegdist(bmc_rabun_t2$Zotu1, method="euclidean"), group=bmc_rabun_t2$Group, type="median") #before using PERMANOVA homoscedasticity of groups must be tested
betadisper_BMC_T2_zotu1

permutest(x=betadisper_BMC_T2_zotu1, pairwise=TRUE, permutations=9999, p.method="bonf") 


perm.betadisper_BMC_T2_zotu1<-adonis(bmc_rabun_t2$Zotu1~Temperature*Inoculation, data=bmc_rabun_t2, permutations=9999, method="euclidean", p.method="bonf") #PERMANOVA
perm.betadisper_BMC_T2_zotu1

library(RVAideMemoire)

pairwise.BMC_T2_zotu1_TEMP<-pairwise.perm.manova(resp=vegdist(bmc_rabun_t2$Zotu1, "euclidean"), fact=bmc_rabun_t2$Temperature, test=NULL, nperm=9999, p.method="bonf")
pairwise.BMC_T2_zotu1_TEMP

pairwise.BMC_T2_zotu1_TREAT<-pairwise.perm.manova(resp=vegdist(bmc_rabun_t2$Zotu1, "euclidean"), fact=bmc_rabun_t2$Inoculation, test=NULL, nperm=9999, p.method="bonf")
pairwise.BMC_T2_zotu1_TREAT

pairwise.BMC_T2_zotu1<-pairwise.perm.manova(resp=vegdist(bmc_rabun_t2$Zotu1, "euclidean"), fact=bmc_rabun_t2$Treatment, test=NULL, nperm=9999, p.method="bonf")
pairwise.BMC_T2_zotu1

###T3

library(vegan)

betadisper_BMC_T3_zotu1<-betadisper(d=vegdist(bmc_rabun_t3$Zotu1, method="euclidean"), group=bmc_rabun_t3$Group, type="median") #before using PERMANOVA homoscedasticity of groups must be tested
betadisper_BMC_T3_zotu1

permutest(x=betadisper_BMC_T3_zotu1, pairwise=TRUE, permutations=9999, p.method="bonf") 


perm.betadisper_BMC_T3_zotu1<-adonis(bmc_rabun_t3$Zotu1~Temperature*Inoculation, data=bmc_rabun_t3, permutations=9999, method="euclidean", p.method="bonf") #PERMANOVA
perm.betadisper_BMC_T3_zotu1

library(RVAideMemoire)

pairwise.BMC_T3_zotu1_TEMP<-pairwise.perm.manova(resp=vegdist(bmc_rabun_t3$Zotu1, "euclidean"), fact=bmc_rabun_t3$Temperature, test=NULL, nperm=9999, p.method="bonf")
pairwise.BMC_T3_zotu1_TEMP

pairwise.BMC_T3_zotu1_TREAT<-pairwise.perm.manova(resp=vegdist(bmc_rabun_t3$Zotu1, "euclidean"), fact=bmc_rabun_t3$Inoculation, test=NULL, nperm=9999, p.method="bonf")
pairwise.BMC_T3_zotu1_TREAT

pairwise.BMC_T3_zotu1<-pairwise.perm.manova(resp=vegdist(bmc_rabun_t3$Zotu1, "euclidean"), fact=bmc_rabun_t3$Treatment, test=NULL, nperm=9999, p.method="bonf")
pairwise.BMC_T3_zotu1

##Reading data and setting libraries:

proxydata <- read.csv("/Users/cardospm/Downloads/BMC experiment 2019/BMC_exp2019_tables_Parameters_Pedro.csv", header=TRUE, stringsAsFactors=FALSE)
require(lme4)
require(lattice)
require(visreg)
require(ggplot2)
library(dplyr)
library("ggpubr")

##Eliminating Vibrio treatments:

proxydata <- proxydata[which(proxydata$Treatment == "BMC" | proxydata$Treatment == "Saline"),names(proxydata) %in% c("Aquarium","Temperature","Treatment","Net_PP","Respiration","Gross_PP","Calcification","Replicate","Timepoint","Days")]

##Dividing data by timepoints:

proxy_t3 <- proxydata[which(proxydata$Timepoint == "T3"),names(proxydata) %in% c("Aquarium","Temperature","Treatment","Net_PP","Respiration","Gross_PP","Calcification","Replicate","Timepoint","Days")]

proxy_t2 <- proxydata[which(proxydata$Timepoint == "T2"),names(proxydata) %in% c("Aquarium","Temperature","Treatment","Net_PP","Respiration","Gross_PP","Calcification","Replicate","Timepoint","Days")]

proxy_t1 <- proxydata[which(proxydata$Timepoint == "T1"),names(proxydata) %in% c("Aquarium","Temperature","Treatment","Net_PP","Respiration","Gross_PP","Calcification","Replicate","Timepoint","Days")]

proxy_t0 <- proxydata[which(proxydata$Timepoint == "T0"),names(proxydata) %in% c("Aquarium","Temperature","Treatment","Net_PP","Respiration","Gross_PP","Calcification","Replicate","Timepoint","Days")]

##Generating Boxplots for Gross Primary Productivity:

gt0 <- ggplot(proxy_t0,aes(Replicate,Gross_PP,fill=Replicate))+geom_boxplot()+labs(x="Group", y=expression("Gross Primary Productivity (µmol O"['2']* " cm"^"-2"*" h"^"-1"*")"))+ylim(-0.1,1.2) + geom_point(position = "jitter", shape= 21, size =3) + stat_summary(fun.y=mean, geom="point", shape=4, size=6)+scale_fill_manual(values = c("#339966","#FF9900","#3399FF","#990000"))

gt0 + theme(
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

gt1 <- ggplot(proxy_t1,aes(Replicate,Gross_PP,fill=Replicate))+geom_boxplot()+labs(x="Group", y=expression("Gross Primary Productivity (µmol O"['2']* " cm"^"-2"*" h"^"-1"*")"))+ylim(-0.1,1.2) + geom_point(position = "jitter", shape= 21, size =3) + stat_summary(fun.y=mean, geom="point", shape=4, size=6)+scale_fill_manual(values = c("#339966","#FF9900","#3399FF","#990000"))

gt1 + theme(
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


gt2 <- ggplot(proxy_t2,aes(Replicate,Gross_PP,fill=Replicate))+geom_boxplot()+labs(x="Group", y=expression("Gross Primary Productivity (µmol O"['2']* " cm"^"-2"*" h"^"-1"*")"))+ylim(-0.1,1.2) + geom_point(position = "jitter", shape= 21, size =3) + stat_summary(fun.y=mean, geom="point", shape=4, size=6)+scale_fill_manual(values = c("#339966","#FF9900","#3399FF","#990000"))

gt2 + theme(
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


gt3 <- ggplot(proxy_t3,aes(Replicate,Gross_PP,fill=Replicate))+geom_boxplot()+labs(x="Group", y=expression("Gross Primary Productivity (µmol O"['2']* " cm"^"-2"*" h"^"-1"*")"))+ylim(-0.1,1.2) + geom_point(position = "jitter", shape= 21, size =3) + stat_summary(fun.y=mean, geom="point", shape=4, size=6)+scale_fill_manual(values = c("#339966","#FF9900","#3399FF","#990000"))

gt3 + theme(
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

gt <- ggplot(proxydata,aes(Timepoint,Gross_PP,fill=Replicate))+geom_boxplot(outlier.shape = NA)+labs(x="Timepoint", y=expression("Gross Primary Productivity (µmol O"['2']* " cm"^"-2"*" h"^"-1"*")"))+ylim(-0.1,1.2) + geom_point(position = position_jitterdodge(jitter.width = 0.1), shape= 21, size =3) + stat_summary(fun.y=mean, geom="point", shape=4, size=4, aes(group=Replicate),position=position_dodge(0.75))+scale_fill_manual(values = c("#339966","#FF9900","#3399FF","#990000"))

gt + theme(
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

##Generating Boxplots for Net Primary Productivity:

nt0 <- ggplot(proxy_t0,aes(Replicate,Net_PP,fill=Replicate))+geom_boxplot()+labs(x="Group", y=expression("Net Primary Productivity (µmol O"['2']* " cm"^"-2"*" h"^"-1"*")")) + geom_point(position = "jitter", shape= 21, size =3) + stat_summary(fun.y=mean, geom="point", shape=4, size=6)+scale_fill_manual(values = c("#339966","#FF9900","#3399FF","#990000"))

nt0 + theme(
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

nt1 <- ggplot(proxy_t1,aes(Replicate,Net_PP,fill=Replicate))+geom_boxplot()+labs(x="Group", y=expression("Net Primary Productivity (µmol O"['2']* " cm"^"-2"*" h"^"-1"*")")) + geom_point(position = "jitter", shape= 21, size =3) + stat_summary(fun.y=mean, geom="point", shape=4, size=6)+scale_fill_manual(values = c("#339966","#FF9900","#3399FF","#990000"))

nt1 + theme(
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


nt2 <- ggplot(proxy_t2,aes(Replicate,Net_PP,fill=Replicate))+geom_boxplot()+labs(x="Group", y=expression("Net Primary Productivity (µmol O"['2']* " cm"^"-2"*" h"^"-1"*")"))+ylim(-0.1,1.2) + geom_point(position = "jitter", shape= 21, size =3) + stat_summary(fun.y=mean, geom="point", shape=4, size=6)+scale_fill_manual(values = c("#339966","#FF9900","#3399FF","#990000"))

nt2 + theme(
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


nt3 <- ggplot(proxy_t3,aes(Replicate,Net_PP,fill=Replicate))+geom_boxplot()+labs(x="Group", y=expression("Net Primary Productivity (µmol O"['2']* " cm"^"-2"*" h"^"-1"*")"))+ylim(-0.1,1.2) + geom_point(position = "jitter", shape= 21, size =3) + stat_summary(fun.y=mean, geom="point", shape=4, size=6)+scale_fill_manual(values = c("#339966","#FF9900","#3399FF","#990000"))

nt3 + theme(
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

nt <- ggplot(proxydata,aes(Timepoint,Net_PP,fill=Replicate))+geom_boxplot(outlier.shape = NA)+labs(x="Timepoint", y=expression("Net Primary Productivity (µmol O"['2']* " cm"^"-2"*" h"^"-1"*")")) + geom_point(position = position_jitterdodge(jitter.width = 0.1), shape= 21, size =3) + stat_summary(fun.y=mean, geom="point", shape=4, size=4, aes(group=Replicate),position=position_dodge(0.75))+scale_fill_manual(values = c("#339966","#FF9900","#3399FF","#990000"))

nt + theme(
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

##Generating Boxplots for Respiration:

rt0 <- ggplot(proxy_t0,aes(Replicate,Respiration,fill=Replicate))+geom_boxplot()+labs(x="Group", y=expression("Respiration (µmol O"['2']* " cm"^"-2"*" h"^"-1"*")")) + geom_point(position = "jitter", shape= 21, size =3) + stat_summary(fun.y=mean, geom="point", shape=4, size=6)+scale_fill_manual(values = c("#339966","#FF9900","#3399FF","#990000"))

rt0 + theme(
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

rt1 <- ggplot(proxy_t1,aes(Replicate,Respiration,fill=Replicate))+geom_boxplot()+labs(x="Group", y=expression("Respiration (µmol O"['2']* " cm"^"-2"*" h"^"-1"*")")) + geom_point(position = "jitter", shape= 21, size =3) + stat_summary(fun.y=mean, geom="point", shape=4, size=6)+scale_fill_manual(values = c("#339966","#FF9900","#3399FF","#990000"))

rt1 + theme(
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


rt2 <- ggplot(proxy_t2,aes(Replicate,Respiration,fill=Replicate))+geom_boxplot()+labs(x="Group", y=expression("Respiration (µmol O"['2']* " cm"^"-2"*" h"^"-1"*")"))+ylim(-0.1,1.2) + geom_point(position = "jitter", shape= 21, size =3) + stat_summary(fun.y=mean, geom="point", shape=4, size=6)+scale_fill_manual(values = c("#339966","#FF9900","#3399FF","#990000"))

rt2 + theme(
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


rt3 <- ggplot(proxy_t3,aes(Replicate,Respiration,fill=Replicate))+geom_boxplot()+labs(x="Group", y=expression("Respiration (µmol O"['2']* " cm"^"-2"*" h"^"-1"*")"))+ylim(-0.1,1.2) + geom_point(position = "jitter", shape= 21, size =3) + stat_summary(fun.y=mean, geom="point", shape=4, size=6)+scale_fill_manual(values = c("#339966","#FF9900","#3399FF","#990000"))

rt3 + theme(
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

rt <- ggplot(proxydata,aes(Timepoint,Respiration,fill=Replicate))+geom_boxplot(outlier.shape = NA)+labs(x="Timepoint", y=expression("Respiration (µmol O"['2']* " cm"^"-2"*" h"^"-1"*")")) + geom_point(position = position_jitterdodge(jitter.width = 0.1), shape= 21, size =3) + stat_summary(fun.y=mean, geom="point", shape=4, size=4, aes(group=Replicate),position=position_dodge(0.75))+scale_fill_manual(values = c("#339966","#FF9900","#3399FF","#990000"))

rt + theme(
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

##Generating Boxplots for Calcification:

ct <- ggplot(proxydata,aes(Timepoint,Calcification,fill=Replicate))+geom_boxplot(outlier.shape = NA)+labs(x="Timepoint", y=expression("Calcification (µmol CaCO"['3']* " cm"^"-2"*" h"^"-1"*")"))+scale_fill_manual(values = c("#339966","#FF9900","#3399FF","#990000"))+ylim(-0.3,0.8) + geom_point(position = position_jitterdodge(jitter.width = 0.1), shape= 21, size =3) + stat_summary(fun.y=mean, geom="point", shape=4, size=4, aes(group=Replicate),position=position_dodge(0.75))

ct + theme(
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

ggqqplot(proxydata$Calcification)

ggqqplot(proxydata$Gross_PP)

shapiro.test(proxydata$Calcification)

shapiro.test(proxydata$Gross_PP)

## Testing Hypothesis:

#PERMANOVA

##Calcification

###T0

library(vegan)

betadisper_BMC_T0_calc<-betadisper(d=vegdist(proxy_t0$Calcification, method="euclidean"), group=proxy_t0$Replicate, type="median") #before using PERMANOVA homoscedasticity of groups must be tested
betadisper_BMC_T0_calc

permutest(x=betadisper_BMC_T0_calc, pairwise=TRUE, permutations=9999, p.method="bonf") 


perm.BMC_T0_calc<-adonis(proxy_t0$Calcification~Temperature*Treatment, data=proxy_t0, permutations=9999, method="euclidean", p.method="bonf") #PERMANOVA
perm.BMC_T0_calc

library(RVAideMemoire)

pairwise.BMC_T0_calc_TEMP<-pairwise.perm.manova(resp=vegdist(proxy_t0$Calcification, "euclidean"), fact=proxy_t0$Temperature, test=NULL, nperm=9999, p.method="bonf")
pairwise.BMC_T0_calc_TEMP

pairwise.BMC_T0_calc_TREAT<-pairwise.perm.manova(resp=vegdist(proxy_t0$Calcification, "euclidean"), fact=proxy_t0$Treatment, test=NULL, nperm=9999, p.method="bonf")
pairwise.BMC_T0_calc_TREAT

pairwise.BMC_T0_calc<-pairwise.perm.manova(resp=vegdist(proxy_t0$Calcification, "euclidean"), fact=proxy_t0$Replicate, test=NULL, nperm=9999, p.method="bonf")
pairwise.BMC_T0_calc

###T1

betadisper_BMC_T1_calc<-betadisper(d=vegdist(proxy_t1$Calcification, method="euclidean"), group=proxy_t1$Replicate, type="median") #before using PERMANOVA homoscedasticity of groups must be tested
betadisper_BMC_T1_calc

permutest(x=betadisper_BMC_T1_calc, pairwise=TRUE, permutations=9999, p.method="bonf") 


perm.BMC_T1_calc<-adonis(proxy_t1$Calcification~Temperature*Treatment, data=proxy_t1, permutations=9999, method="euclidean", p.method="bonf") #PERMANOVA
perm.BMC_T1_calc

library(RVAideMemoire)

pairwise.BMC_T1_calc_TEMP<-pairwise.perm.manova(resp=vegdist(proxy_t1$Calcification, "euclidean"), fact=proxy_t1$Temperature, test=NULL, nperm=9999, p.method="bonf")
pairwise.BMC_T1_calc_TEMP

pairwise.BMC_T1_calc_TREAT<-pairwise.perm.manova(resp=vegdist(proxy_t1$Calcification, "euclidean"), fact=proxy_t1$Treatment, test=NULL, nperm=9999, p.method="bonf")
pairwise.BMC_T1_calc_TREAT

pairwise.BMC_T1_calc<-pairwise.perm.manova(resp=vegdist(proxy_t1$Calcification, "euclidean"), fact=proxy_t1$Replicate, test=NULL, nperm=9999, p.method="bonf")
pairwise.BMC_T1_calc

###T2

betadisper_BMC_T2_calc<-betadisper(d=vegdist(proxy_t2$Calcification, method="euclidean"), group=proxy_t2$Replicate, type="median") #before using PERMANOVA homoscedasticity of groups must be tested
betadisper_BMC_T2_calc

permutest(x=betadisper_BMC_T2_calc, pairwise=TRUE, permutations=9999, p.method="bonf") 


perm.BMC_T2_calc<-adonis(proxy_t2$Calcification~Temperature*Treatment, data=proxy_t2, permutations=9999, method="euclidean", p.method="bonf") #PERMANOVA
perm.BMC_T2_calc

library(RVAideMemoire)

pairwise.BMC_T2_calc_TEMP<-pairwise.perm.manova(resp=vegdist(proxy_t2$Calcification, "euclidean"), fact=proxy_t2$Temperature, test=NULL, nperm=9999, p.method="bonf")
pairwise.BMC_T2_calc_TEMP

pairwise.BMC_T2_calc_TREAT<-pairwise.perm.manova(resp=vegdist(proxy_t2$Calcification, "euclidean"), fact=proxy_t2$Treatment, test=NULL, nperm=9999, p.method="bonf")
pairwise.BMC_T2_calc_TREAT

pairwise.BMC_T2_calc<-pairwise.perm.manova(resp=vegdist(proxy_t2$Calcification, "euclidean"), fact=proxy_t2$Replicate, test=NULL, nperm=9999, p.method="bonf")
pairwise.BMC_T2_calc

###T3

betadisper_BMC_T3_calc<-betadisper(d=vegdist(proxy_t3$Calcification, method="euclidean"), group=proxy_t3$Replicate, type="median") #before using PERMANOVA homoscedasticity of groups must be tested
betadisper_BMC_T3_calc

permutest(x=betadisper_BMC_T3_calc, pairwise=TRUE, permutations=9999, p.method="bonf") 


perm.BMC_T3_calc<-adonis(proxy_t3$Calcification~Temperature*Treatment, data=proxy_t3, permutations=9999, method="euclidean", p.method="bonf") #PERMANOVA
perm.BMC_T3_calc

library(RVAideMemoire)

pairwise.BMC_T3_calc_TEMP<-pairwise.perm.manova(resp=vegdist(proxy_t3$Calcification, "euclidean"), fact=proxy_t3$Temperature, test=NULL, nperm=9999, p.method="bonf")
pairwise.BMC_T3_calc_TEMP

pairwise.BMC_T3_calc_TREAT<-pairwise.perm.manova(resp=vegdist(proxy_t3$Calcification, "euclidean"), fact=proxy_t3$Treatment, test=NULL, nperm=9999, p.method="bonf")
pairwise.BMC_T3_calc_TREAT

pairwise.BMC_T3_calc<-pairwise.perm.manova(resp=vegdist(proxy_t3$Calcification, "euclidean"), fact=proxy_t3$Replicate, test=NULL, nperm=9999, p.method="bonf")
pairwise.BMC_T3_calc

##Gross Primary Productivity

###T0

library(vegan)

betadisper_BMC_T0_gpp<-betadisper(d=vegdist(proxy_t0$Gross_PP, method="euclidean"), group=proxy_t0$Replicate, type="median") #before using PERMANOVA homoscedasticity of groups must be tested
betadisper_BMC_T0_gpp

permutest(x=betadisper_BMC_T0_gpp, pairwise=TRUE, permutations=9999, p.method="bonf") 


perm.BMC_T0_gpp<-adonis(proxy_t0$Gross_PP~Temperature*Treatment, data=proxy_t0, permutations=9999, method="euclidean", p.method="bonf") #PERMANOVA
perm.BMC_T0_gpp

library(RVAideMemoire)

pairwise.BMC_T0_gpp_TEMP<-pairwise.perm.manova(resp=vegdist(proxy_t0$Gross_PP, "euclidean"), fact=proxy_t0$Temperature, test=NULL, nperm=9999, p.method="bonf")
pairwise.BMC_T0_gpp_TEMP

pairwise.BMC_T0_gpp_TREAT<-pairwise.perm.manova(resp=vegdist(proxy_t0$Gross_PP, "euclidean"), fact=proxy_t0$Treatment, test=NULL, nperm=9999, p.method="bonf")
pairwise.BMC_T0_gpp_TREAT

pairwise.BMC_T0_gpp<-pairwise.perm.manova(resp=vegdist(proxy_t0$Gross_PP, "euclidean"), fact=proxy_t0$Replicate, test=NULL, nperm=9999, p.method="bonf")
pairwise.BMC_T0_gpp

###T1

betadisper_BMC_T1_gpp<-betadisper(d=vegdist(proxy_t1$Gross_PP, method="euclidean"), group=proxy_t1$Replicate, type="median") #before using PERMANOVA homoscedasticity of groups must be tested
betadisper_BMC_T1_gpp

permutest(x=betadisper_BMC_T1_gpp, pairwise=TRUE, permutations=9999, p.method="bonf") 


perm.BMC_T1_gpp<-adonis(proxy_t1$Gross_PP~Temperature*Treatment, data=proxy_t1, permutations=9999, method="euclidean", p.method="bonf") #PERMANOVA
perm.BMC_T1_gpp

library(RVAideMemoire)

pairwise.BMC_T1_gpp_TEMP<-pairwise.perm.manova(resp=vegdist(proxy_t1$Gross_PP, "euclidean"), fact=proxy_t1$Temperature, test=NULL, nperm=9999, p.method="bonf")
pairwise.BMC_T1_gpp_TEMP

pairwise.BMC_T1_gpp_TREAT<-pairwise.perm.manova(resp=vegdist(proxy_t1$Gross_PP, "euclidean"), fact=proxy_t1$Treatment, test=NULL, nperm=9999, p.method="bonf")
pairwise.BMC_T1_gpp_TREAT

pairwise.BMC_T1_gpp<-pairwise.perm.manova(resp=vegdist(proxy_t1$Gross_PP, "euclidean"), fact=proxy_t1$Replicate, test=NULL, nperm=9999, p.method="bonf")
pairwise.BMC_T1_gpp

###T2

betadisper_BMC_T2_gpp<-betadisper(d=vegdist(proxy_t2$Gross_PP, method="euclidean"), group=proxy_t2$Replicate, type="median") #before using PERMANOVA homoscedasticity of groups must be tested
betadisper_BMC_T2_gpp

permutest(x=betadisper_BMC_T2_gpp, pairwise=TRUE, permutations=9999, p.method="bonf") 


perm.BMC_T2_gpp<-adonis(proxy_t2$Gross_PP~Temperature*Treatment, data=proxy_t2, permutations=9999, method="euclidean", p.method="bonf") #PERMANOVA
perm.BMC_T2_gpp

library(RVAideMemoire)

pairwise.BMC_T2_gpp_TEMP<-pairwise.perm.manova(resp=vegdist(proxy_t2$Gross_PP, "euclidean"), fact=proxy_t2$Temperature, test=NULL, nperm=9999, p.method="bonf")
pairwise.BMC_T2_gpp_TEMP

pairwise.BMC_T2_gpp_TREAT<-pairwise.perm.manova(resp=vegdist(proxy_t2$Gross_PP, "euclidean"), fact=proxy_t2$Treatment, test=NULL, nperm=9999, p.method="bonf")
pairwise.BMC_T2_gpp_TREAT

pairwise.BMC_T2_gpp<-pairwise.perm.manova(resp=vegdist(proxy_t2$Gross_PP, "euclidean"), fact=proxy_t2$Replicate, test=NULL, nperm=9999, p.method="bonf")
pairwise.BMC_T2_gpp

##Samples not homoscedastic, need to test with non-parametric:
library(FSA)

wilcox.test(bmc_rabun_t3_heat$Zotu42 ~ sample_data(bmc_rabun_t3_heat)$Treatment)
wilcox.test(bmc_rabun_t3_con$Zotu42 ~ sample_data(bmc_rabun_t3_con)$Treatment)

kruskal.test(bmc_rabun_t3$Zotu42 ~ sample_data(bmc_rabun_t3)$Treatment)

pairwise.wilcox.test(bmc_rabun_t3$Zotu42, bmc_rabun_t3$Treatment, p.adjust.method = 'BH')

dunnTest(Zotu42 ~ Treatment,
         data = bmc_rabun_t2,
         method = "bh"
)

###T3

betadisper_BMC_T3_gpp<-betadisper(d=vegdist(proxy_t3$Gross_PP, method="euclidean"), group=proxy_t3$Replicate, type="median") #before using PERMANOVA homoscedasticity of groups must be tested
betadisper_BMC_T3_gpp

permutest(x=betadisper_BMC_T3_gpp, pairwise=TRUE, permutations=9999, p.method="bonf") 


perm.BMC_T3_gpp<-adonis(proxy_t3$Gross_PP~Temperature*Treatment, data=proxy_t3, permutations=9999, method="euclidean", p.method="bonf") #PERMANOVA
perm.BMC_T3_gpp

library(RVAideMemoire)

pairwise.BMC_T3_gpp_TEMP<-pairwise.perm.manova(resp=vegdist(proxy_t3$Gross_PP, "euclidean"), fact=proxy_t3$Temperature, test=NULL, nperm=9999, p.method="bonf")
pairwise.BMC_T3_gpp_TEMP

pairwise.BMC_T3_gpp_TREAT<-pairwise.perm.manova(resp=vegdist(proxy_t3$Gross_PP, "euclidean"), fact=proxy_t3$Treatment, test=NULL, nperm=9999, p.method="bonf")
pairwise.BMC_T3_gpp_TREAT

pairwise.BMC_T3_gpp<-pairwise.perm.manova(resp=vegdist(proxy_t3$Gross_PP, "euclidean"), fact=proxy_t3$Replicate, test=NULL, nperm=9999, p.method="bonf")
pairwise.BMC_T3_gpp

##Net Primary Productivity

###T0

library(vegan)

betadisper_BMC_T0_npp<-betadisper(d=vegdist(proxy_t0$Net_PP, method="euclidean"), group=proxy_t0$Replicate, type="median") #before using PERMANOVA homoscedasticity of groups must be tested
betadisper_BMC_T0_npp

permutest(x=betadisper_BMC_T0_npp, pairwise=TRUE, permutations=9999, p.method="bonf") 


perm.BMC_T0_npp<-adonis(proxy_t0$Net_PP~Temperature*Treatment, data=proxy_t0, permutations=9999, method="euclidean", p.method="bonf") #PERMANOVA
perm.BMC_T0_npp

library(RVAideMemoire)

pairwise.BMC_T0_npp_TEMP<-pairwise.perm.manova(resp=vegdist(proxy_t0$Net_PP, "euclidean"), fact=proxy_t0$Temperature, test=NULL, nperm=9999, p.method="bonf")
pairwise.BMC_T0_npp_TEMP

pairwise.BMC_T0_npp_TREAT<-pairwise.perm.manova(resp=vegdist(proxy_t0$Net_PP, "euclidean"), fact=proxy_t0$Treatment, test=NULL, nperm=9999, p.method="bonf")
pairwise.BMC_T0_npp_TREAT

pairwise.BMC_T0_npp<-pairwise.perm.manova(resp=vegdist(proxy_t0$Net_PP, "euclidean"), fact=proxy_t0$Replicate, test=NULL, nperm=9999, p.method="bonf")
pairwise.BMC_T0_npp

###T1

betadisper_BMC_T1_npp<-betadisper(d=vegdist(proxy_t1$Net_PP, method="euclidean"), group=proxy_t1$Replicate, type="median") #before using PERMANOVA homoscedasticity of groups must be tested
betadisper_BMC_T1_npp

permutest(x=betadisper_BMC_T1_npp, pairwise=TRUE, permutations=9999, p.method="bonf") 


perm.BMC_T1_npp<-adonis(proxy_t1$Net_PP~Temperature*Treatment, data=proxy_t1, permutations=9999, method="euclidean", p.method="bonf") #PERMANOVA
perm.BMC_T1_npp

library(RVAideMemoire)

pairwise.BMC_T1_npp_TEMP<-pairwise.perm.manova(resp=vegdist(proxy_t1$Net_PP, "euclidean"), fact=proxy_t1$Temperature, test=NULL, nperm=9999, p.method="bonf")
pairwise.BMC_T1_npp_TEMP

pairwise.BMC_T1_npp_TREAT<-pairwise.perm.manova(resp=vegdist(proxy_t1$Net_PP, "euclidean"), fact=proxy_t1$Treatment, test=NULL, nperm=9999, p.method="bonf")
pairwise.BMC_T1_npp_TREAT

pairwise.BMC_T1_npp<-pairwise.perm.manova(resp=vegdist(proxy_t1$Net_PP, "euclidean"), fact=proxy_t1$Replicate, test=NULL, nperm=9999, p.method="bonf")
pairwise.BMC_T1_npp

###T2

betadisper_BMC_T2_npp<-betadisper(d=vegdist(proxy_t2$Net_PP, method="euclidean"), group=proxy_t2$Replicate, type="median") #before using PERMANOVA homoscedasticity of groups must be tested
betadisper_BMC_T2_npp

permutest(x=betadisper_BMC_T2_npp, pairwise=TRUE, permutations=9999, p.method="bonf") 


perm.BMC_T2_npp<-adonis(proxy_t2$Net_PP~Temperature*Treatment, data=proxy_t2, permutations=9999, method="euclidean", p.method="bonf") #PERMANOVA
perm.BMC_T2_npp

library(RVAideMemoire)

pairwise.BMC_T2_npp_TEMP<-pairwise.perm.manova(resp=vegdist(proxy_t2$Net_PP, "euclidean"), fact=proxy_t2$Temperature, test=NULL, nperm=9999, p.method="bonf")
pairwise.BMC_T2_npp_TEMP

pairwise.BMC_T2_npp_TREAT<-pairwise.perm.manova(resp=vegdist(proxy_t2$Net_PP, "euclidean"), fact=proxy_t2$Treatment, test=NULL, nperm=9999, p.method="bonf")
pairwise.BMC_T2_npp_TREAT

pairwise.BMC_T2_npp<-pairwise.perm.manova(resp=vegdist(proxy_t2$Net_PP, "euclidean"), fact=proxy_t2$Replicate, test=NULL, nperm=9999, p.method="bonf")
pairwise.BMC_T2_npp

###T3

betadisper_BMC_T3_npp<-betadisper(d=vegdist(proxy_t3$Net_PP, method="euclidean"), group=proxy_t3$Replicate, type="median") #before using PERMANOVA homoscedasticity of groups must be tested
betadisper_BMC_T3_npp

permutest(x=betadisper_BMC_T3_npp, pairwise=TRUE, permutations=9999, p.method="bonf") 


perm.BMC_T3_npp<-adonis(proxy_t3$Net_PP~Temperature*Treatment, data=proxy_t3, permutations=9999, method="euclidean", p.method="bonf") #PERMANOVA
perm.BMC_T3_npp

library(RVAideMemoire)

pairwise.BMC_T3_npp_TEMP<-pairwise.perm.manova(resp=vegdist(proxy_t3$Net_PP, "euclidean"), fact=proxy_t3$Temperature, test=NULL, nperm=9999, p.method="bonf")
pairwise.BMC_T3_npp_TEMP

pairwise.BMC_T3_npp_TREAT<-pairwise.perm.manova(resp=vegdist(proxy_t3$Net_PP, "euclidean"), fact=proxy_t3$Treatment, test=NULL, nperm=9999, p.method="bonf")
pairwise.BMC_T3_npp_TREAT

pairwise.BMC_T3_npp<-pairwise.perm.manova(resp=vegdist(proxy_t3$Net_PP, "euclidean"), fact=proxy_t3$Replicate, test=NULL, nperm=9999, p.method="bonf")
pairwise.BMC_T3_npp

##Respiration

###T0

library(vegan)

betadisper_BMC_T0_res<-betadisper(d=vegdist(proxy_t0$Respiration, method="euclidean"), group=proxy_t0$Replicate, type="median") #before using PERMANOVA homoscedasticity of groups must be tested
betadisper_BMC_T0_res

permutest(x=betadisper_BMC_T0_res, pairwise=TRUE, permutations=9999, p.method="bonf") 


perm.BMC_T0_res<-adonis(proxy_t0$Respiration~Temperature*Treatment, data=proxy_t0, permutations=9999, method="euclidean", p.method="bonf") #PERMANOVA
perm.BMC_T0_res

library(RVAideMemoire)

pairwise.BMC_T0_res_TEMP<-pairwise.perm.manova(resp=vegdist(proxy_t0$Respiration, "euclidean"), fact=proxy_t0$Temperature, test=NULL, nperm=9999, p.method="bonf")
pairwise.BMC_T0_res_TEMP

pairwise.BMC_T0_res_TREAT<-pairwise.perm.manova(resp=vegdist(proxy_t0$Respiration, "euclidean"), fact=proxy_t0$Treatment, test=NULL, nperm=9999, p.method="bonf")
pairwise.BMC_T0_res_TREAT

pairwise.BMC_T0_res<-pairwise.perm.manova(resp=vegdist(proxy_t0$Respiration, "euclidean"), fact=proxy_t0$Replicate, test=NULL, nperm=9999, p.method="bonf")
pairwise.BMC_T0_res

###T1

betadisper_BMC_T1_res<-betadisper(d=vegdist(proxy_t1$Respiration, method="euclidean"), group=proxy_t1$Replicate, type="median") #before using PERMANOVA homoscedasticity of groups must be tested
betadisper_BMC_T1_res

permutest(x=betadisper_BMC_T1_res, pairwise=TRUE, permutations=9999, p.method="bonf") 


perm.BMC_T1_res<-adonis(proxy_t1$Respiration~Temperature*Treatment, data=proxy_t1, permutations=9999, method="euclidean", p.method="bonf") #PERMANOVA
perm.BMC_T1_res

library(RVAideMemoire)

pairwise.BMC_T1_res_TEMP<-pairwise.perm.manova(resp=vegdist(proxy_t1$Respiration, "euclidean"), fact=proxy_t1$Temperature, test=NULL, nperm=9999, p.method="bonf")
pairwise.BMC_T1_res_TEMP

pairwise.BMC_T1_res_TREAT<-pairwise.perm.manova(resp=vegdist(proxy_t1$Respiration, "euclidean"), fact=proxy_t1$Treatment, test=NULL, nperm=9999, p.method="bonf")
pairwise.BMC_T1_res_TREAT

pairwise.BMC_T1_res<-pairwise.perm.manova(resp=vegdist(proxy_t1$Respiration, "euclidean"), fact=proxy_t1$Replicate, test=NULL, nperm=9999, p.method="bonf")
pairwise.BMC_T1_res

###T2

betadisper_BMC_T2_res<-betadisper(d=vegdist(proxy_t2$Respiration, method="euclidean"), group=proxy_t2$Replicate, type="median") #before using PERMANOVA homoscedasticity of groups must be tested
betadisper_BMC_T2_res

permutest(x=betadisper_BMC_T2_res, pairwise=TRUE, permutations=9999, p.method="bonf") 


perm.BMC_T2_res<-adonis(proxy_t2$Respiration~Temperature*Treatment, data=proxy_t2, permutations=9999, method="euclidean", p.method="bonf") #PERMANOVA
perm.BMC_T2_res

library(RVAideMemoire)

pairwise.BMC_T2_res_TEMP<-pairwise.perm.manova(resp=vegdist(proxy_t2$Respiration, "euclidean"), fact=proxy_t2$Temperature, test=NULL, nperm=9999, p.method="bonf")
pairwise.BMC_T2_res_TEMP

pairwise.BMC_T2_res_TREAT<-pairwise.perm.manova(resp=vegdist(proxy_t2$Respiration, "euclidean"), fact=proxy_t2$Treatment, test=NULL, nperm=9999, p.method="bonf")
pairwise.BMC_T2_res_TREAT

pairwise.BMC_T2_res<-pairwise.perm.manova(resp=vegdist(proxy_t2$Respiration, "euclidean"), fact=proxy_t2$Replicate, test=NULL, nperm=9999, p.method="bonf")
pairwise.BMC_T2_res

###T3

betadisper_BMC_T3_res<-betadisper(d=vegdist(proxy_t3$Respiration, method="euclidean"), group=proxy_t3$Replicate, type="median") #before using PERMANOVA homoscedasticity of groups must be tested
betadisper_BMC_T3_res

permutest(x=betadisper_BMC_T3_res, pairwise=TRUE, permutations=9999, p.method="bonf") 


perm.BMC_T3_res<-adonis(proxy_t3$Respiration~Temperature*Treatment, data=proxy_t3, permutations=9999, method="euclidean", p.method="bonf") #PERMANOVA
perm.BMC_T3_res

library(RVAideMemoire)

pairwise.BMC_T3_res_TEMP<-pairwise.perm.manova(resp=vegdist(proxy_t3$Respiration, "euclidean"), fact=proxy_t3$Temperature, test=NULL, nperm=9999, p.method="bonf")
pairwise.BMC_T3_res_TEMP

pairwise.BMC_T3_res_TREAT<-pairwise.perm.manova(resp=vegdist(proxy_t3$Respiration, "euclidean"), fact=proxy_t3$Treatment, test=NULL, nperm=9999, p.method="bonf")
pairwise.BMC_T3_res_TREAT

pairwise.BMC_T3_res<-pairwise.perm.manova(resp=vegdist(proxy_t3$Respiration, "euclidean"), fact=proxy_t3$Replicate, test=NULL, nperm=9999, p.method="bonf")
pairwise.BMC_T3_res
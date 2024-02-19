## BMC experiment bacteria community analysis

# Preparing the data
# Using the data file "pedro_data.txt" and one of the taxonomy files, we will go through phyloseq to look at the data and make figures. First, make sure you are in the right folder by running:
getwd()
setwd("/Users/cardospm/Downloads/BMC experiment 2019/BMC_16s")

# read OTU data, calling this 'dat' 
# use the 'stringsAsFactors' argument to read OTU names as character
dat <- read.delim('pedro_data.txt', stringsAsFactors=F)

dat$timepoint <- factor(dat$timepoint, levels=c('T0','T1','T2', 'T3', 'none'))
dat$inoculation <- factor(dat$inoculation, levels=c('Control', 'BMC', 'none'))
dat$temperature <- factor(dat$temperature, levels=c('26C','31C','none'))
dat$treatment <- factor(dat$treatment, levels=c('CS','CB','HS','HB','none'))
dat$Sample_or_Control <- factor(dat$Sample_or_Control, levels=c("Control Sample","True Sample"))
str(dat[, 1:6])

# read in the taxonomy, use 'stringsAsFactors=F' to keep text as 'character'
tax <- read.delim('pedro_taxonomy_SILVA_cleaned.txt', stringsAsFactors=F)

## Creating phyloseq object
# We need to specify the different components to create a [phyloseq](https://github.com/joey711/phyloseq) object.

# load phyloseq, dplyr(to manipulate data)
library(phyloseq)
library(dplyr)
library(ggplot2)
library(microbiome)
library(BiocManager)
library(DESeq2)
library(decontam)

# extract the OTU data from 'dat' and convert to matrix format
# using 'starts_with' with 'select' allows us to chose all columns containing OTU counts

otus <- as.matrix(select(dat, starts_with('Zotu')))

# rows need to be named with the sample names from 'dat', which we can do directly because they are in the same order
rownames(otus) <- dat$XOTUID

# extract the sample data from from 'dat' and keep as dataframe format (mix of character and numeric variables)
samps <- select(dat, XOTUID:treatment | inoculation | temperature | timepoint | Sample_or_Control)

# rows need to be named with the sample names from 'dat', which we can do directly because they are in the same order
rownames(samps) <- dat$XOTUID

# extract the taxonomy info for each OTU from 'tax' and convert to matrix format
taxonomy <- as.matrix(select(tax, kingdom:species))

# rows need to be named with the OTU names from 'tax', which we can do directly because they are in the same order
rownames(taxonomy) <- tax$XOTUID

# merge the three objects into a single phyloseq object using 'merge_phyloseq'
# each function nexted within the call to 'merge_phyloseq' creates a special object for that type of data
# because of how the OTU matrix is oriented, we have to specify 'taxa_are_rows=F' 
phy <- merge_phyloseq(otu_table(otus, taxa_are_rows=F), 
                sample_data(samps),
                tax_table(taxonomy))

# remove extra objects to keep workspace tidy
rm(otus, samps, taxonomy)

# decontaminate contaminant zOTUs

sample_data(phy)$is.neg <- sample_data(phy)$Sample_or_Control == "Control Sample"
contamdf.prev <- isContaminant(phy, method="prevalence", neg="is.neg")
table(contamdf.prev$contaminant)
phy.noncontam <- prune_taxa(!contamdf.prev$contaminant, phy)
phy.noncontam
phy = phy.noncontam


# remove the Blank sample
phy <- prune_samples(sample_names(phy) != "B1_S140", phy)
phy <- prune_samples(sample_names(phy) != "B2_S141", phy)

# remove archaea/eukarya
phy <- subset_taxa(phy, kingdom=="Bacteria")

# remove OTUs with abundance zero
phy <- prune_samples(sample_sums(phy)>0, phy)

# remove OTUs that hit to Chloroplast or Mitochondria (because they are showing eukarya)
phy <- subset_taxa(phy, order != "o_Chloroplast")
phy <- subset_taxa(phy, family != "f_Mitochondria")

##exporting phyloseq data to df

OTUdf <- as.data.frame(otu_table(phy))
TAXdf <- as.data.frame(tax_table(phy))

#write.csv(OTUdf,file='OTUdf.csv',row.names=TRUE)
#write.csv(TAXdf,file='TAXdf.csv',row.names=TRUE)

## Analyzing the data
# Now that we have the phyloseq object, we can look at the different aspects of the data and make figures.

# first load the ggplot2 library
library(ggplot2)

#------------------------------------------------------------------------------------------------

# we can look at the alpha diversity in each of the samples

min_lib <- min(sample_sums(phy))

nsamp = nsamples(phy)
trials = 100

set.seed(3)

r <- rarefy_even_depth(phy, sample.size = min_lib, verbose = FALSE, replace = TRUE)

phy.t0 = subset_samples(r, timepoint == "T0")
phy.t0.heat = subset_samples(phy.t0, temperature == "31C")
phy.t0.con = subset_samples(phy.t0, temperature == "26C")
phy.t0.bmc = subset_samples(phy.t0, inoculation == "BMC")
phy.t0.sal = subset_samples(phy.t0, inoculation == "Control")

phy.t1 = subset_samples(r, timepoint == "T1")
phy.t1.heat = subset_samples(phy.t1, temperature == "31C")
phy.t1.con = subset_samples(phy.t1, temperature == "26C")
phy.t1.bmc = subset_samples(phy.t1, inoculation == "BMC")
phy.t1.sal = subset_samples(phy.t1, inoculation == "Control")

phy.t2 = subset_samples(r, timepoint == "T2")
phy.t2.heat = subset_samples(phy.t2, temperature == "31C")
phy.t2.con = subset_samples(phy.t2, temperature == "26C")
phy.t2.bmc = subset_samples(phy.t2, inoculation == "BMC")
phy.t2.sal = subset_samples(phy.t2, inoculation == "Control")

phy.t3 = subset_samples(r, timepoint == "T3")
phy.t3.heat = subset_samples(phy.t3, temperature == "31C")
phy.t3.con = subset_samples(phy.t3, temperature == "26C")
phy.t3.bmc = subset_samples(phy.t3, inoculation == "BMC")
phy.t3.sal = subset_samples(phy.t3, inoculation == "Control")

rich<-estimate_richness(r)
rich.t0.con<-estimate_richness(phy.t0.con)
rich.t1.con<-estimate_richness(phy.t1.con)
rich.t2.con<-estimate_richness(phy.t2.con)
rich.t3.con<-estimate_richness(phy.t3.con)

rich.t0.heat<-estimate_richness(phy.t0.heat)
rich.t1.heat<-estimate_richness(phy.t1.heat)
rich.t2.heat<-estimate_richness(phy.t2.heat)
rich.t3.heat<-estimate_richness(phy.t3.heat)

rich.t0.bmc<-estimate_richness(phy.t0.bmc)
rich.t1.bmc<-estimate_richness(phy.t1.bmc)
rich.t2.bmc<-estimate_richness(phy.t2.bmc)
rich.t3.bmc<-estimate_richness(phy.t3.bmc)

rich.t0.sal<-estimate_richness(phy.t0.sal)
rich.t1.sal<-estimate_richness(phy.t1.sal)
rich.t2.sal<-estimate_richness(phy.t2.sal)
rich.t3.sal<-estimate_richness(phy.t3.sal)

# we can plot at the alpha diversity in each of the samples
l<-plot_richness(r, x="treatment", measures=c("Shannon"), color="inoculation") +
  geom_boxplot(fill=NA) +
  facet_wrap(~timepoint) +
  theme_bw() +
  theme(axis.text.x=element_text(angle=45,hjust=1,vjust=1)) +
  scale_y_continuous(name="Shannon diversity", limits=c(0,6.5), breaks = seq(0,6,1))
print(l)

# run statistical test to determine the significant difference between groups
hist(rich$Shannon)
shapiro.test(rich$Shannon)

#two-way ANOVA if samples are normal
res.aov <- aov(rich$Shannon ~ sample_data(phy)$timepoint + sample_data(phy)$treatment)
summary(res.aov)
TukeyHSD(res.aov)

#one-way
res.aov.2 <- aov(rich.t2$Shannon ~ sample_data(phy.t2)$treatment)
summary(res.aov.2)
TukeyHSD(res.aov.2)

#Wilcoxon test for non-normal samples

wilcox.test(rich.t0.con$Shannon ~ sample_data(phy.t0.con)$treatment)
wilcox.test(rich.t1.con$Shannon ~ sample_data(phy.t1.con)$treatment)
wilcox.test(rich.t2.con$Shannon ~ sample_data(phy.t2.con)$treatment)
wilcox.test(rich.t3.con$Shannon ~ sample_data(phy.t3.con)$treatment)

wilcox.test(rich.t0.heat$Shannon ~ sample_data(phy.t0.heat)$treatment)
wilcox.test(rich.t1.heat$Shannon ~ sample_data(phy.t1.heat)$treatment)
wilcox.test(rich.t2.heat$Shannon ~ sample_data(phy.t2.heat)$treatment)
wilcox.test(rich.t3.heat$Shannon ~ sample_data(phy.t3.heat)$treatment)

wilcox.test(rich.t0.bmc$Shannon ~ sample_data(phy.t0.bmc)$treatment)
wilcox.test(rich.t1.bmc$Shannon ~ sample_data(phy.t1.bmc)$treatment)
wilcox.test(rich.t2.bmc$Shannon ~ sample_data(phy.t2.bmc)$treatment)
wilcox.test(rich.t3.bmc$Shannon ~ sample_data(phy.t3.bmc)$treatment)

wilcox.test(rich.t0.sal$Shannon ~ sample_data(phy.t0.sal)$treatment)
wilcox.test(rich.t1.sal$Shannon ~ sample_data(phy.t1.sal)$treatment)
wilcox.test(rich.t2.sal$Shannon ~ sample_data(phy.t2.sal)$treatment)
wilcox.test(rich.t3.sal$Shannon ~ sample_data(phy.t3.sal)$treatment)

#------------------------------------------------------------------------------------------------

## Next, we can make a barplot with the data
# create a dummy variable with both the categories to make a relative abundance plot
sample_data(phy)$dummy <- paste0(sample_data(phy)$timepoint, sample_data(phy)$treatment)
phym = merge_samples(phy, "dummy")

# repair the variables that you just destroyed
sample_data(phym)$timepoint <- levels(sample_data(phy)$timepoint)[get_variable(phym, "timepoint")]
sample_data(phym)$treatment <- levels(sample_data(phy)$treatment)[get_variable(phym, "treatment")]

# to do relative abundance (out of 100%), we need to also transform the counts to percents
phy.prop <- transform_sample_counts(phym, function(otu) 100 * otu/sum(otu))

# make the variables ordered factors so that they are in the correct order when plotting
sample_data(phy.prop)$treatment <- factor(sample_data(phy.prop)$treatment, levels = c("CS","CB","HS","HB"))

# convert phyloseq object to dataframe to make figures in ggplot2
phy.propdf<-psmelt(phy.prop)

#keep the same levels on the x axis
phy.propdf$family <- factor(phy.propdf$family,levels=rev(unique(phy.propdf$family)))

#you can define the colors
#colours <- c("#808080", "#0075DC","#993F00","#4C005C","#2BCE48","#FFCC99","#F0A3FF","#94FFB5","#8F7C00","#9DCC00","#C20088","#003380","#FFA405","#FFA8BB","#426600","#FF0010","#5EF1F2","#00998F","#740AFF","#990000","#FF8000")

# plot
m <- ggplot(phy.propdf, aes(x=treatment, fill = family, y=Abundance)) +
  facet_wrap(~timepoint) +
  geom_bar(aes(color=family), stat = 'identity', position = 'stack', colour = NA) +
  #scale_fill_manual(values = colours) +
  scale_y_continuous(expand = c(0,0), limits=c(0,100)) + theme(strip.background = element_rect(fill="black")) +
  theme_bw() +
  theme(axis.text.x=element_text(angle=45,hjust=1,vjust=1))
print(m)

## We can also make a barplot with just the top taxa at each level. Below shows the top ten order. The numbers can be changed.

# agglomerate at order level
phyg <- tax_glom(phy.prop, "family", NArm = FALSE)

# get just the top 10 order in the samples
top10otus = names(sort(taxa_sums(phyg), TRUE)[1:10])
taxtab10 = cbind(tax_table(phyg), phy10 = NA)
taxtab10[top10otus, "phy10"] <- as(tax_table(phyg)[top10otus, "family"], 
    "character")
tax_table(phyg) <- tax_table(taxtab10)

# alternatively get all the order that had abundances >1% in the samples
#x = taxa_sums(phyg) 
#keeptaxa = taxa_names(phyg)[which((x / sum(x)) > 0.01)]
#phyg1 = prune_taxa(keeptaxa, phyg)

# prune samples
phyg10 = prune_taxa(top10otus, phyg)

# make the variables ordered factors so that they are in the correct order when plotting
sample_data(phyg10)$treatment <- factor(sample_data(phyg10)$treatment, levels = c("CS", "CB","HS","HB","none"))

# convert phyloseq object to dataframe to make figures in ggplot2
phyg10df<-psmelt(phyg10)

#keep the same levels on the x axis
phyg10df$phy10 <- factor(phyg10df$phy10,levels=rev(unique(phyg10df$phy10)))

#you can define the colors
colours <- c("#FFA8BB","#426600","#FF0010","#5EF1F2","#00998F","#740AFF","#990000","#FF8000", "#993F00","#4C005C","#2BCE48")
# "#808080", "#0075DC","#FFCC99","#F0A3FF","#94FFB5","#8F7C00","#9DCC00","#C20088","#003380","#FFA405",

# plot
n <- ggplot(phyg10df, aes(x=treatment, fill = phy10, y=Abundance)) +
  facet_wrap(~timepoint) +
  geom_bar(aes(color=phy10), stat = 'identity', position = 'stack', colour = NA) +
  scale_fill_manual(values = colours) +
  scale_y_continuous(expand = c(0,0), limits=c(0,100)) + theme(strip.background = element_rect(fill="black")) +
  theme_bw() +
  theme(axis.text.x=element_text(angle=45,hjust=1,vjust=1))
print(n)

#------------------------------------------------------------------------------------------------

## Now we want to check the beta diversity by making a nMDS or PCoA, and run some tests.


##T0

# we can just use a subset of the samples so that we don't clog the plot
phy.t0 = subset_samples(phy, timepoint == "T0")
phy.t0.heat = subset_samples(phy.t0, temperature == "31C")
phy.t0.con = subset_samples(phy.t0, temperature == "26C")
phy.t0.bmc = subset_samples(phy.t0, inoculation == "BMC")
phy.t0.sal = subset_samples(phy.t0, inoculation == "Control")


# this is to normalize the abundances of the samples (here used as a percent)
phy.t0.norm = transform_sample_counts(phy.t0, function(x) 100 * x/sum(x))
phy.t0.norm.heat = transform_sample_counts(phy.t0.heat, function(x) 100 * x/sum(x))
phy.t0.norm.con = transform_sample_counts(phy.t0.con, function(x) 100 * x/sum(x))
phy.t0.norm.sal = transform_sample_counts(phy.t0.sal, function(x) 100 * x/sum(x))
phy.t0.norm.bmc = transform_sample_counts(phy.t0.bmc, function(x) 100 * x/sum(x))

# calculate the distance matrix
phy.dist.t0 = phyloseq::distance(phy.t0.norm, method="bray")
phy.dist.t0.heat = phyloseq::distance(phy.t0.norm.heat, method="bray")
phy.dist.t0.con = phyloseq::distance(phy.t0.norm.con, method="bray")
phy.dist.t0.sal = phyloseq::distance(phy.t0.norm.sal, method="bray")
phy.dist.t0.bmc = phyloseq::distance(phy.t0.norm.bmc, method="bray")

# calculate the ordination (using NMDS, PCoA, PCA) using distance (bray, euclidean)
ord <- ordinate(phy.t0.norm, method="PCoA", distance=phy.dist.t0)
ord <- ordinate(phy.t0.norm.heat, method="PCoA", distance=phy.dist.t0)
ord <- ordinate(phy.t0.norm.con, method="PCoA", distance=phy.dist.t0)

# see the eigen value of the axes
plot_scree(ord)

# plot
p <- plot_ordination(phy.t0.norm, ord, color = "treatment", shape = "timepoint")+
  geom_point(size=2, alpha=0.8)+
  scale_shape_manual(values=c(16, 17,18))+
  theme_bw()+
  scale_color_manual(values = c('#d7191c','#fdae61','#abdda4','#2b83ba'))+
  stat_ellipse()
plot(p)

# calculate the beta diversity
library(vegan)
adonis(phy.dist.t0 ~ sample_data(phy.t0.norm)$treatment, permutations = 9999)
adonis(phy.dist.t0.heat ~ sample_data(phy.t0.norm.heat)$treatment, permutations = 9999)
adonis(phy.dist.t0.con ~ sample_data(phy.t0.norm.con)$treatment, permutations = 9999)
adonis(phy.dist.t0.sal ~ sample_data(phy.t0.norm.sal)$treatment, permutations = 9999)
adonis(phy.dist.t0.bmc ~ sample_data(phy.t0.norm.bmc)$treatment, permutations = 9999)

##T1

# we can just use a subset of the samples so that we don't clog the plot
phy.t1 = subset_samples(phy, timepoint == "T1")
phy.t1.heat = subset_samples(phy.t1, temperature == "31C")
phy.t1.con = subset_samples(phy.t1, temperature == "26C")
phy.t1.bmc = subset_samples(phy.t1, inoculation == "BMC")
phy.t1.sal = subset_samples(phy.t1, inoculation == "Control")

# this is to normalize the abundances of the samples (here used as a percent)
phy.t1.norm = transform_sample_counts(phy.t1, function(x) 100 * x/sum(x))
phy.t1.norm.heat = transform_sample_counts(phy.t1.heat, function(x) 100 * x/sum(x))
phy.t1.norm.con = transform_sample_counts(phy.t1.con, function(x) 100 * x/sum(x))
phy.t1.norm.sal = transform_sample_counts(phy.t1.sal, function(x) 100 * x/sum(x))
phy.t1.norm.bmc = transform_sample_counts(phy.t1.bmc, function(x) 100 * x/sum(x))

# calculate the distance matrix
phy.dist.t1 = phyloseq::distance(phy.t1.norm, method="bray")
phy.dist.t1.heat = phyloseq::distance(phy.t1.norm.heat, method="bray")
phy.dist.t1.con = phyloseq::distance(phy.t1.norm.con, method="bray")
phy.dist.t1.sal = phyloseq::distance(phy.t1.norm.sal, method="bray")
phy.dist.t1.bmc = phyloseq::distance(phy.t1.norm.bmc, method="bray")

# calculate the ordination (using NMDS, PCoA, PCA) using distance (bray, euclidean)
ord <- ordinate(phy.t1.norm, method="PCoA", distance=phy.dist.t1)
ord <- ordinate(phy.t1.norm.heat, method="PCoA", distance=phy.dist.t1)
ord <- ordinate(phy.t1.norm.con, method="PCoA", distance=phy.dist.t1)

# see the eigen value of the axes
plot_scree(ord)

# plot
p <- plot_ordination(phy.t1.norm, ord, color = "treatment", shape = "timepoint")+
  geom_point(size=2, alpha=0.8)+
  scale_shape_manual(values=c(16, 17,18))+
  theme_bw()+
  scale_color_manual(values = c('#d7191c','#fdae61','#abdda4','#2b83ba'))+
  stat_ellipse()
plot(p)

# calculate the beta diversity
library(vegan)
adonis(phy.dist.t1 ~ sample_data(phy.t1.norm)$treatment, permutations = 9999)
adonis(phy.dist.t1.heat ~ sample_data(phy.t1.norm.heat)$treatment, permutations = 9999)
adonis(phy.dist.t1.con ~ sample_data(phy.t1.norm.con)$treatment, permutations = 9999)
adonis(phy.dist.t1.sal ~ sample_data(phy.t1.norm.sal)$treatment, permutations = 9999)
adonis(phy.dist.t1.bmc ~ sample_data(phy.t1.norm.bmc)$treatment, permutations = 9999)

##T2

# we can just use a subset of the samples so that we don't clog the plot
phy.t2 = subset_samples(phy, timepoint == "T2")
phy.t2.heat = subset_samples(phy.t2, temperature == "31C")
phy.t2.con = subset_samples(phy.t2, temperature == "26C")
phy.t2.bmc = subset_samples(phy.t2, inoculation == "BMC")
phy.t2.sal = subset_samples(phy.t2, inoculation == "Control")

# this is to normalize the abundances of the samples (here used as a percent)
phy.t2.norm = transform_sample_counts(phy.t2, function(x) 100 * x/sum(x))
phy.t2.norm.heat = transform_sample_counts(phy.t2.heat, function(x) 100 * x/sum(x))
phy.t2.norm.con = transform_sample_counts(phy.t2.con, function(x) 100 * x/sum(x))
phy.t2.norm.sal = transform_sample_counts(phy.t2.sal, function(x) 100 * x/sum(x))
phy.t2.norm.bmc = transform_sample_counts(phy.t2.bmc, function(x) 100 * x/sum(x))

# calculate the distance matrix
phy.dist.t2 = phyloseq::distance(phy.t2.norm, method="bray")
phy.dist.t2.heat = phyloseq::distance(phy.t2.norm.heat, method="bray")
phy.dist.t2.con = phyloseq::distance(phy.t2.norm.con, method="bray")
phy.dist.t2.sal = phyloseq::distance(phy.t2.norm.sal, method="bray")
phy.dist.t2.bmc = phyloseq::distance(phy.t2.norm.bmc, method="bray")

# calculate the ordination (using NMDS, PCoA, PCA) using distance (bray, euclidean)
ord <- ordinate(phy.t2.norm, method="PCoA", distance=phy.dist.t2)
ord <- ordinate(phy.t2.norm.heat, method="PCoA", distance=phy.dist.t2)
ord <- ordinate(phy.t2.norm.con, method="PCoA", distance=phy.dist.t2)

# see the eigen value of the axes
plot_scree(ord)

# plot
p <- plot_ordination(phy.t2.norm, ord, color = "treatment", shape = "timepoint")+
  geom_point(size=2, alpha=0.8)+
  scale_shape_manual(values=c(16, 17,18))+
  theme_bw()+
  scale_color_manual(values = c('#d7191c','#fdae61','#abdda4','#2b83ba'))+
  stat_ellipse()
plot(p)

# calculate the beta diversity
library(vegan)
adonis(phy.dist.t2 ~ sample_data(phy.t2.norm)$treatment, permutations = 9999)
adonis(phy.dist.t2.heat ~ sample_data(phy.t2.norm.heat)$treatment, permutations = 9999)
adonis(phy.dist.t2.con ~ sample_data(phy.t2.norm.con)$treatment, permutations = 9999)
adonis(phy.dist.t2.sal ~ sample_data(phy.t2.norm.sal)$treatment, permutations = 9999)
adonis(phy.dist.t2.bmc ~ sample_data(phy.t2.norm.bmc)$treatment, permutations = 9999)



###T3

# we can just use a subset of the samples so that we don't clog the plot
phy.t3 = subset_samples(phy, timepoint == "T3")
phy.t3.heat = subset_samples(phy.t3, temperature == "31C")
phy.t3.con = subset_samples(phy.t3, temperature == "26C")
phy.t3.bmc = subset_samples(phy.t3, inoculation == "BMC")
phy.t3.sal = subset_samples(phy.t3, inoculation == "Control")

# this is to normalize the abundances of the samples (here used as a percent)
phy.t3.norm = transform_sample_counts(phy.t3, function(x) 100 * x/sum(x))
phy.t3.norm.heat = transform_sample_counts(phy.t3.heat, function(x) 100 * x/sum(x))
phy.t3.norm.con = transform_sample_counts(phy.t3.con, function(x) 100 * x/sum(x))
phy.t3.norm.sal = transform_sample_counts(phy.t3.sal, function(x) 100 * x/sum(x))
phy.t3.norm.bmc = transform_sample_counts(phy.t3.bmc, function(x) 100 * x/sum(x))

# calculate the distance matrix
phy.dist.t3 = phyloseq::distance(phy.t3.norm, method="bray")
phy.dist.t3.heat = phyloseq::distance(phy.t3.norm.heat, method="bray")
phy.dist.t3.con = phyloseq::distance(phy.t3.norm.con, method="bray")
phy.dist.t3.sal = phyloseq::distance(phy.t3.norm.sal, method="bray")
phy.dist.t3.bmc = phyloseq::distance(phy.t3.norm.bmc, method="bray")

# calculate the ordination (using NMDS, PCoA, PCA) using distance (bray, euclidean)
ord <- ordinate(phy.t3.norm, method="PCoA", distance=phy.dist.t3)
ord <- ordinate(phy.t3.norm.heat, method="PCoA", distance=phy.dist.t3)
ord <- ordinate(phy.t3.norm.con, method="PCoA", distance=phy.dist.t3)

# see the eigen value of the axes
plot_scree(ord)

# plot
p <- plot_ordination(phy.t3.norm, ord, color = "treatment", shape = "inoculation")+
  geom_point(size=2, alpha=0.8)+
  scale_shape_manual(values=c(16, 17,18))+
  theme_bw()+
  scale_color_manual(values = c('#d7191c','#fdae61','#abdda4','#2b83ba'))+
  stat_ellipse()
plot(p)

# calculate the beta diversity
library(vegan)
adonis(phy.dist.t3 ~ sample_data(phy.t3.norm)$treatment, permutations = 9999)
adonis(phy.dist.t3.heat ~ sample_data(phy.t3.norm.heat)$treatment, permutations = 9999)
adonis(phy.dist.t3.con ~ sample_data(phy.t3.norm.con)$treatment, permutations = 9999)
adonis(phy.dist.t3.sal ~ sample_data(phy.t3.norm.sal)$treatment, permutations = 9999)
adonis(phy.dist.t3.bmc ~ sample_data(phy.t3.norm.bmc)$treatment, permutations = 9999)

## All timepoints

# this is to normalize the abundances of the samples (here used as a percent)
phy.norm = transform_sample_counts(phy, function(x) 100 * x/sum(x))

# calculate the distance matrix
phy.dist = phyloseq::distance(phy.norm, method="bray")

# calculate the ordination (using NMDS, PCoA, PCA) using distance (bray, euclidean)
ord <- ordinate(phy.norm, method="PCoA", distance=phy.dist)

# see the eigen value of the axes
plot_scree(ord)

# plot
p <- plot_ordination(phy.norm, ord, color = "treatment", shape = "timepoint")+
  geom_point(size=2, alpha=0.8)+
  scale_shape_manual(values=c(15,16, 17,18))+
  theme_bw()+
  scale_color_manual(values = c('#d7191c','#fdae61','#abdda4','#2b83ba'))
plot(p)

##Now for finding Indicator species:

rm(list=ls())

library (indicspecies)
library (stats)
library(tidyverse)

getwd()
setwd ("/Users/cardospm/Downloads/BMC experiment 2019/BMC_16s")
pc = read.csv("/Users/cardospm/Downloads/BMC experiment 2019/BMC_16s/pedro_data_no_euk_silva_t3_heat.csv", header=TRUE)

abund = pc[,6:ncol(pc)]
treat = pc$group

inv = multipatt(abund, treat, duleg = TRUE, func = "r.g", control = how(nperm=1000))
#inv = multipatt(abund, treat, restcomb = c(2,6,7), func = "r.g", control = how(nperm=1000))
#inv = multipatt(abund, treat, func = "r.g", control = how(nperm=1000))

p <- capture.output(summary(inv, indvalcomp=TRUE, alpha = 1))


write.csv(p,file='Indicsp_summary.csv',row.names=TRUE)

#comparing GROUPS of potential indicators instead of individual species
wetcomb = combinespecies(wetland, max.order = 2)$XC
dim(wetcomb)

indvalspcomb = multipatt(wetcomb, groups, duleg = TRUE, control = how(nperm=999))
summary(indvalspcomb, indvalcomp = TRUE)


p.adjust(c$p, method = "BH")


library(data.table)
indisp<- multipatt(abund, treat, duleg = TRUE, func = "r.g", control = how(nperm=1000))
#extract table of stats
indisp.sign<-as.data.table(indisp$sign, keep.rownames=TRUE)
#add adjusted p-value
indisp.sign[ ,p.value.bh:=p.adjust(p.value, method="BH")]
#now can select only the indicators with adjusted significant p-values
indisp.sign[p.value.bh<=0.05, ]


#######Now for indicator species analysis:


# Preparing the data
# Using the data file "pedro_data_silva_noeuk.txt" and one of the taxonomy files, we will go through phyloseq to look at the data and make figures. First, make sure you are in the right folder by running:
getwd()
setwd("/Users/cardospm/Downloads/BMC experiment 2019/BMC_16s")

library(tidyverse)
library (indicspecies)
library(dplyr)

# read OTU data, calling this 'dat' 
# use the 'stringsAsFactors' argument to read OTU names as character
data <- read.delim('pedro_data_silva_noeuk.txt', stringsAsFactors=F)

# read in the taxonomy, use 'stringsAsFactors=F' to keep text as 'character'
tax <- read.delim('pedro_taxonomy_SILVA_cleaned_2.txt', stringsAsFactors=F)

colnames(tax)

testdply <- data %>%
  pivot_longer(cols = starts_with("Zotu"), names_to = "id", values_to = "abd") %>%
  left_join(tax, by = join_by(id == XOTU)) %>%
  group_by(sample, timepoint, temperature, treatment, inoculation, id, genus) %>%
  summarize( abd = sum(abd)) %>%
  pivot_wider(names_from = id, values_from = abd)

dply <- inner_join(data, tax, 
                   by = join_by(sample == XOTU))

#testdply[, c(1:3, grep("Grand.*", colnames(testdply)))]



##--------------------------------------------------------------------------------

##Testing for T3 heat

ind.sp = subset(data, timepoint == 'T3')
ind.sp = subset(ind.sp, temperature == '31C')
ind.sp.abund = ind.sp[,6:ncol(ind.sp)]
ind.sp.rabund=t(apply(ind.sp.abund, 1, function(x)(x)*100/(sum(x))))
treat = ind.sp$treatment

indisp<- multipatt(ind.sp.rabund, treat, duleg = TRUE, func = "r.g", max.order = 1, control = how(nperm=9999))

#extract table of stats
indisp.sign<-data.table::as.data.table(indisp$sign, keep.rownames=TRUE)
#add adjusted p-value
indisp.sign[ ,p.value.bh:=p.adjust(p.value, method="BH")]
#now can select only the indicators with adjusted significant p-values
indtable <- indisp.sign[p.value<=0.05]
write.csv(indtable,"/Users/cardospm/Downloads/BMC experiment 2019/BMC_16s/Indic_species_T3_heat.csv", row.names = TRUE)

##--------------------------------------------------------------------------------

fill_isa_both <- tax$phylum

indplot <- ggplot(indtable, aes(x= rn, y=stat, fill = rn)) +
  geom_bar(stat = "identity") +
  coord_polar("x", start = 0, clip = "off") +
  geom_text(aes(label = rn, y = ifelse(stat>0, stat, -2))) +
  facet_wrap(~ index, ncol = 3) + 
  #scale_fill_manual(values = fill_isa_both, name = "Phylum") +
  guides(fill = guide_legend(ncol = 4)) +
  labs(x = NULL, y = "Number of indicator ASVs") +
  theme(axis.text.x = element_blank(), legend.position = "none")
indplot

##Testing for T3 26 °C

ind.sp = subset(data, timepoint == 'T3')
ind.sp = subset(ind.sp, temperature == '26C')
ind.sp.abund = ind.sp[,6:ncol(ind.sp)]
ind.sp.rabund=t(apply(ind.sp.abund, 1, function(x)(x)*100/(sum(x))))
treat = ind.sp$treatment

indisp<- multipatt(ind.sp.rabund, treat, duleg = TRUE, func = "r.g", max.order = 1, control = how(nperm=9999))

#extract table of stats
indisp.sign<-data.table::as.data.table(indisp$sign, keep.rownames=TRUE)
#add adjusted p-value
indisp.sign[ ,p.value.bh:=p.adjust(p.value, method="BH")]
#now can select only the indicators with adjusted significant p-values
indtable <- indisp.sign[p.value<=0.05]
write.csv(indtable,"/Users/cardospm/Downloads/BMC experiment 2019/BMC_16s/Indic_species_T3_26C.csv", row.names = TRUE)

##--------------------------------------------------------------------------------

indplot <- ggplot(indtable, aes(x= rn, y=stat, fill = rn)) +
  geom_bar(stat = "identity") +
  coord_polar("x", start = 0, clip = "off") +
  geom_text(aes(label = rn, y = ifelse(stat>0, stat, -2))) +
  facet_wrap(~ index, ncol = 3) + 
  #scale_fill_manual(values = fill_isa_both, name = "Phylum") +
  guides(fill = guide_legend(ncol = 4)) +
  labs(x = NULL, y = "Number of indicator ASVs") +
  theme(axis.text.x = element_blank(), legend.position = "none")
indplot

##Testing for T2 heat

ind.sp = subset(data, timepoint == 'T2')
ind.sp = subset(ind.sp, temperature == '31C')
ind.sp.abund = ind.sp[,6:ncol(ind.sp)]
ind.sp.rabund=t(apply(ind.sp.abund, 1, function(x)(x)*100/(sum(x))))
treat = ind.sp$treatment

indisp<- multipatt(ind.sp.rabund, treat, duleg = TRUE, func = "r.g", max.order = 1, control = how(nperm=9999))

#extract table of stats
indisp.sign<-data.table::as.data.table(indisp$sign, keep.rownames=TRUE)
#add adjusted p-value
indisp.sign[ ,p.value.bh:=p.adjust(p.value, method="BH")]
#now can select only the indicators with adjusted significant p-values
indtable <- indisp.sign[p.value<=0.05]
write.csv(indtable,"/Users/cardospm/Downloads/BMC experiment 2019/BMC_16s/Indic_species_T2_heat.csv", row.names = TRUE)

##--------------------------------------------------------------------------------

indplot <- ggplot(indtable, aes(x= rn, y=stat, fill = rn)) +
  geom_bar(stat = "identity") +
  coord_polar("x", start = 0, clip = "off") +
  geom_text(aes(label = rn, y = ifelse(stat>0, stat, -2))) +
  facet_wrap(~ index, ncol = 3) + 
  #scale_fill_manual(values = fill_isa_both, name = "Phylum") +
  guides(fill = guide_legend(ncol = 4)) +
  labs(x = NULL, y = "Number of indicator ASVs") +
  theme(axis.text.x = element_blank(), legend.position = "none")
indplot

##Testing for T2 26 °C

ind.sp = subset(data, timepoint == 'T2')
ind.sp = subset(ind.sp, temperature == '26C')
ind.sp.abund = ind.sp[,6:ncol(ind.sp)]
ind.sp.rabund=t(apply(ind.sp.abund, 1, function(x)(x)*100/(sum(x))))
treat = ind.sp$treatment

indisp<- multipatt(ind.sp.rabund, treat, duleg = TRUE, func = "r.g", max.order = 1, control = how(nperm=9999))

#extract table of stats
indisp.sign<-data.table::as.data.table(indisp$sign, keep.rownames=TRUE)
#add adjusted p-value
indisp.sign[ ,p.value.bh:=p.adjust(p.value, method="BH")]
#now can select only the indicators with adjusted significant p-values
indtable <- indisp.sign[p.value<=0.05]
write.csv(indtable,"/Users/cardospm/Downloads/BMC experiment 2019/BMC_16s/Indic_species_T2_26C.csv", row.names = TRUE)

##--------------------------------------------------------------------------------


indplot <- ggplot(indtable, aes(x= rn, y=stat, fill = rn)) +
  geom_bar(stat = "identity") +
  coord_polar("x", start = 0, clip = "off") +
  geom_text(aes(label = NA, y = ifelse(stat>0, stat, -2))) +
  facet_wrap(~ index, ncol = 3) + 
  #scale_fill_manual(values = fill_isa_both, name = "Phylum") +
  guides(fill = guide_legend(ncol = 4)) +
  labs(x = NULL, y = "Number of indicator ZOTUs") +
  theme(axis.text.x = element_blank(), legend.position = "none")
indplot

BETA diversity Swabs
#install.packages("raster")
library(vegan) 
library(ggplot2)
library(ggpubr)
library(data.table)
library(phyloseq)
library(qiime2R)
library(tidyr)
library(naniar)
library(raster)
library(plotly)
#transpose function

setwd("~/Desktop/eunice/PhD/IBDProb/") 
metadata <- read.csv("NG_METADATAIBDADLIBCorrect.csv", na.strings = c("","NA"), header=TRUE)
str(metadata)
metadata$TreatmentGroup <- as.factor(metadata$TreatmentGroup)
levels(metadata$TreatmentGroup) <- list("DSS+H2O"="DSS+H2O", "DSS+LbcVec"="DSS+LbcVec","DSS+LbcLapLin"="DSS+LbcLapLin")
rownames(metadata) = metadata$ID
metadata$dayCat = metadata$day
metadata$dayCat = as.factor(metadata$dayCat)
levels(metadata$TreatmentGroup) <- list("DSS+H2O"="DSS+H2O", "DSS+LbcVec"="DSS+LbcVec","DSS+LbcLapLin"="DSS+LbcLapLin")

ASVs <- read_qza("rarefied_table.qza") #548 ASVs, 34 samples
ASV_s <- as.data.frame(ASVs$data)
ASV_table <- as.data.frame(ASVs$data) #548 ASVs
ASV_table$ASVnos <- paste0("ASV", 1:nrow(ASV_table))
ASV_table$ASVstring <- rownames(ASV_table)
rownames(ASV_table) <- ASV_table$ASVnos ##We change the ASV name created in Qiime to ASVn
ASVkey <- ASV_table[, (ncol(ASV_table)-1):ncol(ASV_table)] #the key withe the names
ASV_table <- ASV_table[,-(ncol(ASV_table)-1):-ncol(ASV_table)]
ASV_table <- t(ASV_table)

#change names on the ASV table 
ASV_table <- as.data.frame(ASV_table)
ASV_table <- merge(metadata, ASV_table, by.x = "ID", by.y = 0)
row.names(ASV_table) <- ASV_table$name
ASV_table <- ASV_table[,-c(1:8)]
ASV_table <- as.matrix(ASV_table)
x = as.data.frame(t(ASV_table[c(1,2),]))
sum(x$F1-DAY0) #34412

rownames(metadata) = metadata$name

#Taxonomy of each OTU
##Adding taxonomy
#Taxonomy of each OTU
tax <- read_qza("taxonomy.qza")
tax <- as.data.frame(tax$data)
tax2 = separate(tax, Taxon, into = c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species"), sep=";")
#This warning means that some cells are empty and that R is replacing empty cells with NA. Now there are other cells that are unclassified that  say, for example `s__` or `g__`. All these need to be removed and replaced with `NA`. 
#All this is OK except that in future use of the taxonomy table, these ASVs will be ignored because they are not classified. Why are ASVs not classified? Its because there is not a close enough match in the database. Just because there is not a good match in the database does not mean they don’t exist, so I wanted to make sure this data was not lost. So in my new code, from lines 300 – 334 I make it so that ASVs that are unclassified at any level are classified as the lowest taxonomic level for which there is a classification.

#All the strings that need to be removed and replaced with NA
na_strings <- c(" s__", " g__", " f__", " o__", " c__")

tax3 = replace_with_na_all(tax2, condition = ~.x %in% na_strings)

#This code is great because an ASV that is unclassified at a certain level are all listed as `NA`.
#Unfortunately this command changed ou Feature.ID names

#Next, all these `NA` classifications with the last level that was classified
tax3[] <- t(apply(tax3, 1, zoo::na.locf))
tax3 <- as.data.frame(tax3)
row.names(tax3) <- tax3[,1]
tax3 = tax3[,-c(1:2)]
tax.clean <- as.data.frame(tax3)
tax.clean$OTUs <- rownames(tax.clean)
#Would be good to check here to make sure the order of the two data frames was the same. You should do this on your own.

###Remove all the OTUs that don't occur in our OTU.clean data set
tax.final = tax.clean[row.names(tax.clean) %in% row.names(ASV_s),]

##Remove unneccesary information from the taxonomy names
tax.final$Phylum <- sub("d__*", "", tax.final[,1])
tax.final$Phylum <- sub("p__*", "", tax.final[,1])
tax.final$Phylum <- sub("c__*", "", tax.final[,1])
tax.final$Class <- sub("d__*", "", tax.final[,2])
tax.final$Class <- sub("p__*", "", tax.final[,2])
tax.final$Class <- sub("c__*", "", tax.final[,2])
tax.final$Class <- sub("o__*", "", tax.final[,2])
tax.final$Order <- sub("d__*", "", tax.final[,3])
tax.final$Order <- sub("p__*", "", tax.final[,3])
tax.final$Order <- sub("c__*", "", tax.final[,3])
tax.final$Order <- sub("o__*", "", tax.final[,3])
tax.final$Order <- sub("f__*", "", tax.final[,3])
tax.final$Family <- sub("d__*", "", tax.final[,4])
tax.final$Family <- sub("p__*", "", tax.final[,4])
tax.final$Family <- sub("c__*", "", tax.final[,4])
tax.final$Family <- sub("o__*", "", tax.final[,4])
tax.final$Family <- sub("f__*", "", tax.final[,4])
tax.final$Family <- sub("g__*", "", tax.final[,4])
tax.final$Family <- sub("D_9__*", "", tax.final[,4])
tax.final$Genus <- sub("d__*", "", tax.final[,5])
tax.final$Genus <- sub("p__*", "", tax.final[,5])
tax.final$Genus <- sub("c__*", "", tax.final[,5])
tax.final$Genus <- sub("o__*", "", tax.final[,5])
tax.final$Genus <- sub("f__*", "", tax.final[,5])
tax.final$Genus <- sub("g__*", "", tax.final[,5])
tax.final$Genus <- sub("s__*", "", tax.final[,5])
tax.final$Species <- sub("d__*", "", tax.final[,6])
tax.final$Species <- sub("p__*", "", tax.final[,6])
tax.final$Species <- sub("c__*", "", tax.final[,6])
tax.final$Species <- sub("o__*", "", tax.final[,6])
tax.final$Species <- sub("f__*", "", tax.final[,6])
tax.final$Species <- sub("g__*", "", tax.final[,6])
tax.final$Species <- sub("s__*", "", tax.final[,6])
#write.table(tax.final,"taxonomyNasal.txt",sep=",", row.names = FALSE) 
TaxASV <- merge(tax.final, ASVkey, by.x = 0, by.y = "ASVstring")
row.names(TaxASV) <- TaxASV[,10]
TaxASV = TaxASV[,-c(1,10)]
#write.table(TaxASV,"TaxASV.txt",sep=",", row.names = FALSE)

### Creating the Phyloseq Object
OTU.physeq = otu_table(as.matrix(ASV_table), taxa_are_rows=FALSE)
tax.physeq = tax_table(as.matrix(TaxASV))
#meta.physeq = sample_data(meta)
meta.physeq = sample_data(metadata)

#We then merge these into an object of class phyloseq.
physeq_deseq = phyloseq(OTU.physeq, tax.physeq, meta.physeq)
physeq_deseq # [ 548 taxa and 34 samples ]

colnames(tax_table(physeq_deseq))
## Filter any non-baxteria, chloroplast and mitochondria
physeq_deseq %>%
  subset_taxa(Family != "Mitochondria" & 
                Genus != "Mitochondria" &
                Species != "Mitochondria" &
                Order != "Chloroplast" &
                Family != "Chloroplast" &
                Genus != "Chloroplast" &
                Species != "Chloroplast") -> physeq_deseq
physeq_deseq #[ 10759 taxa and 368 samples ] Dairy

prunetable<- phyloseq_to_df(physeq_deseq, addtax = T, addtot = F, addmaxrank = F,
                            sorting = "abundance")

## no mitochondria or chloroplast in the data
NewTax <- prunetable[,c(1:9)]
row.names(NewTax) <- NewTax[,1]
NewTax = NewTax[,-c(1)]
#write.table(NewTax,"NewTax.txt",sep=",", row.names = TRUE)
NewASVTable <- prunetable
NewASVTable <- NewASVTable[,-c(2:9)]
row.names(NewASVTable) <- NewASVTable[,1]
NewASVTable = NewASVTable[,-c(1)]

##Overall day and treatment effect
NewASVTable2 = t(NewASVTable)
NewASVTable2 <- NewASVTable2[rownames(NewASVTable2) %in% rownames(metadata),]

OTU.physeq = otu_table(as.matrix(NewASVTable2), taxa_are_rows=FALSE)
tax.physeq = tax_table(as.matrix(NewTax))
#meta.physeq = sample_data(meta)
meta.physeq = sample_data(metadata)

#We then merge these into an object of class phyloseq.
physeq_deseq = phyloseq(OTU.physeq, tax.physeq, meta.physeq)
physeq_deseq # [ 548 taxa and 34 samples ]

#Permanova test
dist.bray <- phyloseq::distance(physeq_deseq, method = "bray")
dist.bray <- as.dist(dist.bray)

set.seed(783587)
main <- adonis2(dist.bray~ metadata$TreatmentGroup + metadata$dayCat + metadata$TreatmentGroup*metadata$dayCat,  permutations = 999)
main

set.seed(254287)
treat <- adonis2(dist.bray~ metadata$TreatmentGroup,  permutations = 999)
treat

treat1 <- pairwise.adonis2(dist.bray ~ TreatmentGroup,  data=metadata)
treat1 

set.seed(358345)
day <- adonis2(dist.bray~ metadata$dayCat,  permutations = 999)
day

day1 <- pairwise.adonis2(dist.bray ~ dayCat,  data=metadata)
day1 

inte <- pairwise.adonis2(dist.bray ~ TreatmentGroup*dayCat,  data=metadata)
inte 


#plot
ordu.bc <- ordinate(physeq_deseq, "PCoA", "bray")
components <- ordu.bc[["vectors"]]
components <- data.frame(components)
components = components[,-c(4:25)]
components = cbind(components, metadata$TreatmentGroup)
components = cbind(components, metadata$dayCat)
str(components)
colnames(components) <- c("Axis.1", "Axis.2", "Axis.3", "TreatmentGroup","dayCat")
str(components)
legend_order <- c("0", "7", "17")
components$dayCat <- factor(components$dayCat, levels = legend_order)
legend_order <- c("DSS+H2O", "DSS+LbcVec", "DSS+LbcLapLin")
components$TreatmentGroup <- factor(components$TreatmentGroup, levels = legend_order)


fig <- plot_ly(components, x = ~Axis.1 , y = ~Axis.2, z = ~Axis.3, symbol = ~dayCat, symbols = symbols, color = ~as.factor(dayCat), colors = "black")
fig
colors <- c("black", "red", "green")
fig <- add_markers(fig, color = ~TreatmentGroup, colors = colors)
fig

fig <- plot_ly(components, x = ~Axis.1 , y = ~Axis.2, z = ~Axis.3, color = ~as.factor(TreatmentGroup), colors = colors, symbol = ~as.factor(TreatmentGroup), symbols = "square")
fig
symbols <- c("circle", "cross", "diamond")
fig <- add_markers(fig, symbol = ~as.factor(dayCat), symbols = symbols)
fig

#Weighted UniFrac
physeq = phyloseq(OTU.physeq, tax.physeq)
physeq

##Making the tree
library("ape")
random_tree = rtree(ntaxa(physeq), rooted=TRUE, tip.label=taxa_names(physeq))
##Merging the tree with the phyloseq object
physeq = merge_phyloseq(physeq, meta.physeq, random_tree)
physeq 

OTU_Unifrac <- UniFrac(physeq, weighted = TRUE) ## Calculating the weighted unifrac distances

#running permanova for WU
set.seed(386537)
MainWU <- adonis2(OTU_Unifrac~ metadata$TreatmentGroup + metadata$dayCat + metadata$TreatmentGroup*metadata$dayCat, permutations = 999)
MainWU

trt1 <- pairwise.adonis2(OTU_Unifrac ~ TreatmentGroup,  data=metadata)
trt1

day1 <- pairwise.adonis2(OTU_Unifrac ~ dayCat,  data=metadata)
day1

ordu.bc <- ordinate(physeq,"PCoA", "unifrac", weighted=T)
components <- ordu.bc[["vectors"]]
components <- data.frame(components)
components = components[,-c(4:23)]
components = cbind(components, metadata$TreatmentGroup)
components = cbind(components, metadata$dayCat)
str(components)
colnames(components) <- c("Axis.1", "Axis.2", "Axis.3", "TreatmentGroup","dayCat")
str(components)
legend_order <- c("0", "7", "17")
components$dayCat <- factor(components$dayCat, levels = legend_order)
legend_order <- c("DSS+H2O", "DSS+LbcVec", "DSS+LbcLapLin")
components$TreatmentGroup <- factor(components$TreatmentGroup, levels = legend_order)

fig <- plot_ly(components, x = ~Axis.1 , y = ~Axis.2, z = ~Axis.3, symbol = ~dayCat, symbols = symbols, color = ~as.factor(dayCat), colors = "black")
fig
colors <- c("black", "red", "green")
fig <- add_markers(fig, color = ~TreatmentGroup, colors = colors)
fig

fig <- plot_ly(components, x = ~Axis.1 , y = ~Axis.2, z = ~Axis.3, color = ~as.factor(TreatmentGroup), colors = colors, symbol = ~as.factor(TreatmentGroup), symbols = "square")
fig
symbols <- c("circle", "cross", "diamond")
fig <- add_markers(fig, symbol = ~as.factor(dayCat), symbols = symbols)
fig

#Unweighted UniFrac
OTU_UnifracUN <- UniFrac(physeq, weighted = FALSE) ## Calculating the weighted unifrac distances

#running permanova for WU
set.seed(836453)
MainUN <- adonis2(OTU_UnifracUN~ metadata$TreatmentGroup + metadata$dayCat + metadata$TreatmentGroup*metadata$dayCat, permutations = 999)
MainUN

trt2 <- pairwise.adonis2(OTU_UnifracUN ~ TreatmentGroup,  data=metadata)
trt2

day2 <- pairwise.adonis2(OTU_UnifracUN ~ dayCat,  data=metadata)
day2

ordu.UN <- ordinate(physeq,"PCoA", "unifrac", weighted=F)
components <- ordu.UN[["vectors"]]
components <- data.frame(components)
components = components[,-c(4:32)]
components = cbind(components, metadata$TreatmentGroup)
components = cbind(components, metadata$dayCat)
str(components)
colnames(components) <- c("Axis.1", "Axis.2", "Axis.3", "TreatmentGroup","dayCat")
str(components)
legend_order <- c("0", "7", "17")
components$dayCat <- factor(components$dayCat, levels = legend_order)
legend_order <- c("DSS+H2O", "DSS+LbcVec", "DSS+LbcLapLin")
components$TreatmentGroup <- factor(components$TreatmentGroup, levels = legend_order)

fig <- plot_ly(components, x = ~Axis.1 , y = ~Axis.2, z = ~Axis.3, symbol = ~dayCat, symbols = symbols, color = ~as.factor(dayCat), colors = "black")
fig
colors <- c("black", "red", "green")
fig <- add_markers(fig, color = ~TreatmentGroup, colors = colors)
fig

fig <- plot_ly(components, x = ~Axis.1 , y = ~Axis.2, z = ~Axis.3, color = ~as.factor(TreatmentGroup), colors = colors, symbol = ~as.factor(TreatmentGroup), symbols = "square")
fig
symbols <- c("circle", "cross", "diamond")
fig <- add_markers(fig, symbol = ~as.factor(dayCat), symbols = symbols)
fig


##### separation by dat
NewASVTable2 = t(NewASVTable)
day0 = subset(metadata, dayCat == "0")
d0_OTU <- NewASVTable2[rownames(NewASVTable2) %in% rownames(day0),]

day7 = subset(metadata, dayCat == "7")
d7_OTU <- NewASVTable2[rownames(NewASVTable2) %in% rownames(day7),]

day17 = subset(metadata, dayCat == "17")
d17_OTU <- NewASVTable2[rownames(NewASVTable2) %in% rownames(day17),]

## create phyloseq object
### Creating the Phyloseq Object
#day0
OTU.physeqd0 = otu_table(as.matrix(d0_OTU ), taxa_are_rows=FALSE)
tax.physeqd0 = tax_table(as.matrix(NewTax))
#meta.physeq = sample_data(meta)
meta.physeqd0 = sample_data(day0)

#We then merge these into an object of class phyloseq.
physeq_deseqd0 = phyloseq(OTU.physeqd0, tax.physeqd0, meta.physeqd0)
physeq_deseqd0 # [ 548 taxa and 11 samples ]

#day7
OTU.physeqd7 = otu_table(as.matrix(d7_OTU ), taxa_are_rows=FALSE)
tax.physeqd7 = tax_table(as.matrix(NewTax))
#meta.physeq = sample_data(meta)
meta.physeqd7 = sample_data(day7)

#We then merge these into an object of class phyloseq.
physeq_deseqd7 = phyloseq(OTU.physeqd7, tax.physeqd7, meta.physeqd7)
physeq_deseqd7 # [ 548 taxa and 11 samples ]

#day17
OTU.physeqd17 = otu_table(as.matrix(d17_OTU ), taxa_are_rows=FALSE)
tax.physeqd17 = tax_table(as.matrix(NewTax))
#meta.physeq = sample_data(meta)
meta.physeqd17 = sample_data(day17)

#We then merge these into an object of class phyloseq.
physeq_deseqd17 = phyloseq(OTU.physeqd17, tax.physeqd17, meta.physeqd17)
physeq_deseqd17 # [ 548 taxa and 12 samples ]

## Calculating distances based on Bray-curtis and Weighted Unifrac using the physeq object
dist.bray0 <- phyloseq::distance(physeq_deseqd0, method = "bray")
dist.bray0 <- as.dist(dist.bray0)

dist.bray7 <- phyloseq::distance(physeq_deseqd7, method = "bray")
dist.bray7 <- as.dist(dist.bray7)

dist.bray17 <- phyloseq::distance(physeq_deseqd17, method = "bray")
dist.bray17 <- as.dist(dist.bray17)

## permanova Bray Curtis
Md0 <- adonis2(dist.bray0~ day0$TreatmentGroup, permutations = 999)
Md0

PH_BCday0 <- pairwise.adonis2(dist.bray0 ~ TreatmentGroup,  data=day0)
PH_BCday0 

## permanova Bray Curtis
Md7 <- adonis2(dist.bray7~ day7$TreatmentGroup, permutations = 999)
Md7

PH_BCday7 <- pairwise.adonis2(dist.bray7 ~ TreatmentGroup,  data=day7)
PH_BCday7

## permanova Bray Curtis
Md17 <- adonis2(dist.bray17~ day17$TreatmentGroup, permutations = 999)
Md17

PH_BCday17 <- pairwise.adonis2(dist.bray17 ~ TreatmentGroup,  data=day17)
PH_BCday17


###################
## Weighted Unifrac
# We need to create first a tree using OTU and taxa table-- we do this by creating a phyloseq object 
#We then merge these into an object of class phyloseq.
physeqd0 = phyloseq(OTU.physeqd0, tax.physeq)
physeqd0

##Making the tree
library("ape")
random_treed0 = rtree(ntaxa(physeqd0), rooted=TRUE, tip.label=taxa_names(physeqd0))
##Merging the tree with the phyloseq object
physeqd0 = merge_phyloseq(physeqd0, meta.physeqd0, random_treed0)
physeqd0 

OTU_Unifracd0 <- UniFrac(physeqd0, weighted = TRUE) ## Calculating the weighted unifrac distances

#We then merge these into an object of class phyloseq.
physeqd7 = phyloseq(OTU.physeqd7, tax.physeq)
physeqd7

##Making the tree
random_treed7 = rtree(ntaxa(physeqd7), rooted=TRUE, tip.label=taxa_names(physeqd7))
##Merging the tree with the phyloseq object
physeqd7 = merge_phyloseq(physeqd7, meta.physeqd7, random_treed7)
physeqd7 

OTU_Unifracd7 <- UniFrac(physeqd7, weighted = TRUE) ## Calculating the weighted unifrac distances

#We then merge these into an object of class phyloseq.
physeqd17 = phyloseq(OTU.physeqd17, tax.physeq)
physeqd17

##Making the tree
random_treed17 = rtree(ntaxa(physeqd17), rooted=TRUE, tip.label=taxa_names(physeqd17))
##Merging the tree with the phyloseq object
physeqd17 = merge_phyloseq(physeqd17, meta.physeqd17, random_treed17)
physeqd17 

OTU_Unifracd17 <- UniFrac(physeqd17, weighted = TRUE) ## Calculating the weighted unifrac distances

#running permanova for WU
MUd0 <- adonis2(OTU_Unifracd0~ day0$TreatmentGroup, permutations = 999)
MUd0

PH_WUday0 <- pairwise.adonis2(OTU_Unifracd0 ~ TreatmentGroup,  data=day0)
PH_WUday0 

#running permanova for WU
MUd7 <- adonis2(OTU_Unifracd7~ day7$TreatmentGroup, permutations = 999)
MUd7

PH_WUday7 <- pairwise.adonis2(OTU_Unifracd7 ~ TreatmentGroup,  data=day7)
PH_WUday7

#running permanova for WU
MUd17 <- adonis2(OTU_Unifracd17~ day17$TreatmentGroup, permutations = 999)
MUd17

PH_WUday17 <- pairwise.adonis2(OTU_Unifracd17 ~ TreatmentGroup,  data=day17)
PH_WUday17

### Unweigted UniFrac
OTU_UnifracUNd0 <- UniFrac(physeqd0, weighted = FALSE) ## Calculating the weighted unifrac distances
OTU_UnifracUNd7 <- UniFrac(physeqd7, weighted = FALSE) ## Calculating the weighted unifrac distances
OTU_UnifracUNd17 <- UniFrac(physeqd17, weighted = FALSE) ## Calculating the weighted unifrac distances

#running permanova for WU UWEIGTHED
MUd0UN <- adonis2(OTU_UnifracUNd0~ day0$TreatmentGroup, permutations = 999)
MUd0UN

PH_WUday0UN <- pairwise.adonis2(OTU_UnifracUNd0 ~ TreatmentGroup,  data=day0)
PH_WUday0UN

#running permanova for WU
MUd7UN <- adonis2(OTU_UnifracUNd7~ day7$TreatmentGroup, permutations = 999)
MUd7UN

PH_WUday7UN <- pairwise.adonis2(OTU_UnifracUNd7 ~ TreatmentGroup,  data=day7)
PH_WUday7UN

#running permanova for WU
MUd17UN <- adonis2(OTU_UnifracUNd17~ day17$TreatmentGroup, permutations = 999)
MUd17UN

PH_WUday17UN <- pairwise.adonis2(OTU_UnifracUNd17 ~ TreatmentGroup,  data=day17)
PH_WUday17UN

######## creating plots
my_color1 <- c("black", "red", "green")

## Weighted unifrac
rm(wt.unifrac)
set.seed(123)
ordu.wt.uni <- ordinate(physeqd0 , "PCoA", "unifrac", weighted=T)
wt.unifrac <- plot_ordination(physeqd0, 
                              ordu.wt.uni, color="TreatmentGroup") 

set.seed(123)
wt.unifrac <- wt.unifrac + ggtitle("Weighted UniFrac Day 0") + geom_point(size = 2)
wt.unifrac <- wt.unifrac + theme_classic() 
print(wt.unifrac )
a = print(wt.unifrac + 
            scale_color_manual(values = my_color1) +
            theme(legend.text = element_text(size=12)) +
            theme(legend.title = element_text(size = 12, face= "bold")) +
            theme(axis.title.x = element_text(color="black", size=12, face="bold"), axis.title.y = element_text(color="black", size=12, face="bold")) + 
            theme(axis.text.x = element_text(color = "black", size = 12), axis.text.y = element_text(color = "black", size = 12)))
a

#day 7
rm(ordu.wt.unid7)
set.seed(456)
ordu.wt.unid7 <- ordinate(physeqd7 , "PCoA", "unifrac", weighted=T)
wt.unifracd7 <- plot_ordination(physeqd7, 
                              ordu.wt.unid7, color="TreatmentGroup") 

set.seed(456)
wt.unifracd7 <- wt.unifracd7 + ggtitle("Weighted UniFrac Day 7") + geom_point(size = 2)
wt.unifracd7 <- wt.unifracd7 + theme_classic() 
print(wt.unifracd7)
b = print(wt.unifracd7 + 
            scale_color_manual(values = my_color1) +
            theme(legend.text = element_text(size=12)) +
            theme(legend.title = element_text(size = 12, face= "bold")) +
            theme(axis.title.x = element_text(color="black", size=12, face="bold"), axis.title.y = element_text(color="black", size=12, face="bold")) + 
            theme(axis.text.x = element_text(color = "black", size = 12), axis.text.y = element_text(color = "black", size = 12)))
b

#day 17
rm(ordu.wt.unid17)
set.seed(789)
ordu.wt.unid17 <- ordinate(physeqd17 , "PCoA", "unifrac", weighted=T)
wt.unifracd17 <- plot_ordination(physeqd17, 
                                ordu.wt.unid17, color="TreatmentGroup") 

set.seed(789)
wt.unifracd17 <- wt.unifracd17 + ggtitle("Weighted UniFrac Day 17") + geom_point(size = 2)
wt.unifracd17 <- wt.unifracd17 + theme_classic() 
print(wt.unifracd17)
c = print(wt.unifracd17 + 
            scale_color_manual(values = my_color1) +
            theme(legend.text = element_text(size=12)) +
            theme(legend.title = element_text(size = 12, face= "bold")) +
            theme(axis.title.x = element_text(color="black", size=12, face="bold"), axis.title.y = element_text(color="black", size=12, face="bold")) + 
            theme(axis.text.x = element_text(color = "black", size = 12), axis.text.y = element_text(color = "black", size = 12)))
c

## UNWeighted unifrac
rm(wt.unifrac)
set.seed(123)
ordu.wt.uniUN <- ordinate(physeqd0 , "PCoA", "unifrac", weighted=F)
wt.unifracUN <- plot_ordination(physeqd0, 
                                ordu.wt.uniUN , color="TreatmentGroup") 

set.seed(123)
wt.unifracUN <- wt.unifracUN + ggtitle("UnWeighted UniFrac Day 0") + geom_point(size = 2)
wt.unifracUN <- wt.unifracUN + theme_classic() 
print(wt.unifracUN)
d = print(wt.unifracUN + 
            scale_color_manual(values = my_color1) +
            theme(legend.text = element_text(size=12)) +
            theme(legend.title = element_text(size = 12, face= "bold")) +
            theme(axis.title.x = element_text(color="black", size=12, face="bold"), axis.title.y = element_text(color="black", size=12, face="bold")) + 
            theme(axis.text.x = element_text(color = "black", size = 12), axis.text.y = element_text(color = "black", size = 12)))
d

#day 7
rm(wt.unifracdUN7)
set.seed(456)
ordu.wt.uniUNd7 <- ordinate(physeqd7 , "PCoA", "unifrac", weighted=F)
wt.unifracdUN7 <- plot_ordination(physeqd7, 
                                ordu.wt.uniUNd7, color="TreatmentGroup") 

set.seed(456)
wt.unifracdUN7 <- wt.unifracdUN7 + ggtitle("UnWeighted UniFrac Day 7") + geom_point(size = 2)
wt.unifracdUN7 <- wt.unifracdUN7 + theme_classic() 
print(wt.unifracdUN7)
e = print(wt.unifracdUN7 + 
            scale_color_manual(values = my_color1) +
            theme(legend.text = element_text(size=12)) +
            theme(legend.title = element_text(size = 12, face= "bold")) +
            theme(axis.title.x = element_text(color="black", size=12, face="bold"), axis.title.y = element_text(color="black", size=12, face="bold")) + 
            theme(axis.text.x = element_text(color = "black", size = 12), axis.text.y = element_text(color = "black", size = 12)))
b

#day 17
rm(wt.unifracdUN17)
set.seed(789)
ordu.wt.uniUNd17 <- ordinate(physeqd17 , "PCoA", "unifrac", weighted=F)
wt.unifracdUN17 <- plot_ordination(physeqd17, 
                                  ordu.wt.uniUNd17, color="TreatmentGroup") 

set.seed(789)
wt.unifracdUN17 <- wt.unifracdUN17 + ggtitle("UnWeighted UniFrac Day 17") + geom_point(size = 2)
wt.unifracdUN17 <- wt.unifracdUN17 + theme_classic() 
print(wt.unifracdUN17)
f = print(wt.unifracdUN17 + 
            scale_color_manual(values = my_color1) +
            theme(legend.text = element_text(size=12)) +
            theme(legend.title = element_text(size = 12, face= "bold")) +
            theme(axis.title.x = element_text(color="black", size=12, face="bold"), axis.title.y = element_text(color="black", size=12, face="bold")) + 
            theme(axis.text.x = element_text(color = "black", size = 12), axis.text.y = element_text(color = "black", size = 12)))
f

ggarrange(a,b,c,d,e,f, ncol =3, nrow=2, labels = c("A", "B", "C", "D", "E", "F"),
          font.label = list(size = 11))

## Bray_curtis
ordu.bc <- ordinate(physeqd0, "PCoA", "bray")
Bray <- plot_ordination(physeqd0, 
                        ordu.bc, color="TreatmentGroup") 
Bray <- Bray + ggtitle("Bray Curtis Day 0") + geom_point(size = 2)
Bray <- Bray + theme_classic()
x = print(Bray + 
            theme(legend.text = element_text(size=12)) +
            scale_color_manual(values = my_color1) +
            theme(legend.title = element_text(size = 12, face= "bold")) +
            theme(axis.title.x = element_text(color="black", size=12, face="bold"), axis.title.y = element_text(color="black", size=12, face="bold")) + 
            theme(axis.text.x = element_text(color = "black", size = 12), axis.text.y = element_text(color = "black", size = 12)))
x

#day 7
ordu.bcd7 <- ordinate(physeqd7, "PCoA", "bray")
Brayd7 <- plot_ordination(physeqd7, 
                        ordu.bcd7, color="TreatmentGroup") 
Brayd7 <- Brayd7 + ggtitle("Bray Curtis Day 7") + geom_point(size = 2)
Brayd7 <- Brayd7 + theme_classic()
y = print(Brayd7 + 
            theme(legend.text = element_text(size=12)) +
            scale_color_manual(values = my_color1) +
            theme(legend.title = element_text(size = 12, face= "bold")) +
            theme(axis.title.x = element_text(color="black", size=12, face="bold"), axis.title.y = element_text(color="black", size=12, face="bold")) + 
            theme(axis.text.x = element_text(color = "black", size = 12), axis.text.y = element_text(color = "black", size = 12)))

#day 17
ordu.bcd17 <- ordinate(physeqd17, "PCoA", "bray")
Brayd17 <- plot_ordination(physeqd17, 
                          ordu.bcd17, color="TreatmentGroup") 
Brayd17 <- Brayd17 + ggtitle("Bray Curtis Day 17") + geom_point(size = 2)
Brayd17 <- Brayd17 + theme_classic()
z = print(Brayd17 + 
            theme(legend.text = element_text(size=12)) +
            scale_color_manual(values = my_color1) +
            theme(legend.title = element_text(size = 12, face= "bold")) +
            theme(axis.title.x = element_text(color="black", size=12, face="bold"), axis.title.y = element_text(color="black", size=12, face="bold")) + 
            theme(axis.text.x = element_text(color = "black", size = 12), axis.text.y = element_text(color = "black", size = 12)))

ggarrange(a,b,c,d,e,f,x,y,z, ncol =3, nrow=3, labels = c("A", "B", "C", "D", "E", "F", "G", "H", "I"),
          font.label = list(size = 11))





BRD_Wu <- betadisper(OTU_Unifrac, type = c("centroid"), group = metadata$dayCat)
BRD_Wu
boxplot(BRD_Wu)
TukeyHSD(BRD_Wu) ## calculate the distance from the centroids to one group to another
boxplot(BRD_Wu) #dispersion weighted UniFrac significant

#install.packages("usedist")
library(usedist)

pBRD_WU<- permutest(BRD_Wu, permutations = 999)
pBRD_WU# Significant

#Bray_curtis dispersion
BRD_BC <- betadisper(dist.bray, type = c("centroid"), group = metadata$Status)
BRD_BC
boxplot(BRD_BC)

pBRD_BC<- permutest(BRD_BC, permutations = 999)
pBRD_BC# significant 

##O Make the plots with ellipses
## Calculating the centroids of the data based on Bray-curtis and weighted unifrac
#D33 centroids for state
my_color <- c(
  "lightpink", "#56B4E9","orange", "#66CC99","#5E738F","darkseagreen", "olivedrab", "palevioletred",
  "skyblue", "#CBD588","#D14385", "#653936", "#CD9BCD", 
  "#8569D5", "#5E738F","#D1A33D", "#8A7C64", "#599861", "black","lightblue"
)

## Weighted unifrac 
ordu.wt.uni <- ordinate(physeq1 , "PCoA", "unifrac", weighted=T)
wt.unifrac <- plot_ordination(physeq1, 
                              ordu.wt.uni, color="State") 
wuaxis1 <- wt.unifrac[["data"]][["Axis.1"]]
wuaxis1 <- as.data.frame(wuaxis1)
wuaxis1$number <- rownames(wuaxis1)
wuaxis2 <- wt.unifrac[["data"]][["Axis.2"]]
wuaxis2 <- as.data.frame(wuaxis2)
wuaxis2$number <- rownames(wuaxis2)
wuaxis3 <- wt.unifrac[["data"]][["ID"]]
wuaxis3 <- as.data.frame(wuaxis3)
wuaxis3$number <- rownames(wuaxis3)
wuaxis4 <- wt.unifrac[["data"]][["State"]]
wuaxis4 <- as.data.frame(wuaxis4)
wuaxis4$number <- rownames(wuaxis4)

wuaxis <- merge(wuaxis3, wuaxis4, by.x = "number", by.y = "number")
wuaxis <- merge(wuaxis, wuaxis1, by.x = "number", by.y = "number")
wuaxis <- merge(wuaxis, wuaxis2, by.x = "number", by.y = "number")

wt.unifrac <- wt.unifrac + ggtitle("Dairy Weighted UniFrac") + geom_point(size = 2)
wt.unifrac <- wt.unifrac + theme_classic() 
print(wt.unifrac + stat_ellipse())
a = print(wt.unifrac + stat_ellipse() +
            scale_color_manual(values = my_color) +
            theme(legend.text = element_text(size=12)) +
            theme(legend.title = element_text(size = 12, face= "bold")) +
            theme(axis.title.x = element_text(color="black", size=12, face="bold"), axis.title.y = element_text(color="black", size=12, face="bold")) + 
            theme(axis.text.x = element_text(color = "black", size = 12), axis.text.y = element_text(color = "black", size = 12)))
a
##other plot
str(wuaxis)

centroids <- aggregate(wuaxis[,4:5], list(Group=wuaxis$wuaxis4), mean)
colnames(centroids) <- c('BRD','groupX', 'groupY')


#distance betweeen centroids
pointDistance(c(0.03168905, -0.0010942445), c(-0.02609687, 0.0009011425), lonlat=FALSE)

## Bray_curtis
ordu.bc <- ordinate(physeq1, "PCoA", "bray")
Bray <- plot_ordination(physeq1, 
                        ordu.bc, color="State") 
Bray1 <- Bray[["data"]][["Axis.1"]]
Bray1 <- as.data.frame(Bray1)
Bray1$number <- rownames(Bray1)
Bray2 <- Bray[["data"]][["Axis.2"]]
Bray2 <- as.data.frame(Bray2)
Bray2$number <- rownames(Bray2)
Bray3 <- Bray[["data"]][["ID"]]
Bray3 <- as.data.frame(Bray3)
Bray3$number <- rownames(Bray3)
Bray4 <- Bray[["data"]][["State"]]
Bray4 <- as.data.frame(Bray4)
Bray4$number <- rownames(Bray4)

bray <- merge(Bray3, Bray4, by.x = "number", by.y = "number")
bray <- merge(bray, Bray1, by.x = "number", by.y = "number")
bray <- merge(bray, Bray2, by.x = "number", by.y = "number")

Bray <- Bray + ggtitle("Dairy Bray Curtis") + geom_point(size = 2)
Bray <- Bray + theme_classic()
b = print(Bray + stat_ellipse() +
            theme(legend.text = element_text(size=12)) +
            scale_color_manual(values = my_color) +
            theme(legend.title = element_text(size = 12, face= "bold")) +
            theme(axis.title.x = element_text(color="black", size=12, face="bold"), axis.title.y = element_text(color="black", size=12, face="bold")) + 
            theme(axis.text.x = element_text(color = "black", size = 12), axis.text.y = element_text(color = "black", size = 12)))
b

str(bray)
ggarrange(b,a,ncol =2, labels = c("A", "B"),
          font.label = list(size = 11))


### separate by state 
BrayB <- plot_ordination(physeq1, 
                         ordu.bc, color="Status") 
Bray1 <- BrayB[["data"]][["Axis.1"]]
Bray1 <- as.data.frame(Bray1)
Bray1$number <- rownames(Bray1)
Bray2 <- BrayB[["data"]][["Axis.2"]]
Bray2 <- as.data.frame(Bray2)
Bray2$number <- rownames(Bray2)
Bray3 <- BrayB[["data"]][["ID"]]
Bray3 <- as.data.frame(Bray3)
Bray3$number <- rownames(Bray3)
Bray4 <- BrayB[["data"]][["Status"]]
Bray4 <- as.data.frame(Bray4)
Bray4$number <- rownames(Bray4)

bray <- merge(Bray3, Bray4, by.x = "number", by.y = "number")
bray <- merge(bray, Bray1, by.x = "number", by.y = "number")
bray <- merge(bray, Bray2, by.x = "number", by.y = "number")

my_colors <- c("dodgerblue3","goldenrod3")
BrayB <- BrayB + ggtitle("Bray Curtis") + geom_point(size = 2)
BrayB <- BrayB + theme_classic()
print(BrayB + stat_ellipse() + facet_wrap(State~.)+
        theme_bw() + scale_color_manual(values = my_colors) +
        theme(legend.text = element_text(size=12)) +
        theme(legend.title = element_text(size = 12, face= "bold")) +
        theme(axis.title.x = element_text(color="black", size=12, face="bold"), axis.title.y = element_text(color="black", size=12, face="bold")) + 
        theme(axis.text.x = element_text(color = "black", size = 12), axis.text.y = element_text(color = "black", size = 12)) +
        ggtitle("Dairy Farms: Bray-Curtis dissimilarity")
)

#Separate by State- Dairy WU
WUB <- plot_ordination(physeq1, 
                       ordu.wt.uni, color="Status") 
WU1 <- WUB[["data"]][["Axis.1"]]
WU1 <- as.data.frame(WU1)
WU1$number <- rownames(WU1)
WU2 <- WUB[["data"]][["Axis.2"]]
WU2 <- as.data.frame(WU2)
WU2$number <- rownames(WU2)
WU3 <- WUB[["data"]][["ID"]]
WU3 <- as.data.frame(WU3)
WU3$number <- rownames(WU3)
WU4 <- WUB[["data"]][["Status"]]
WU4 <- as.data.frame(WU4)
WU4$number <- rownames(WU4)

WU <- merge(WU3, WU4, by.x = "number", by.y = "number")
WU <- merge(WU, WU1, by.x = "number", by.y = "number")
WU <- merge(WU, WU2, by.x = "number", by.y = "number")

my_colors <- c("dodgerblue3","goldenrod3")
WUB <- WUB + ggtitle("Weigthed UniFrac") + geom_point(size = 2)
WUB <- WUB + theme_classic()
print(WUB + stat_ellipse() + facet_wrap(State~.)+
        theme_bw() + scale_color_manual(values = my_colors) +
        theme(legend.text = element_text(size=12)) +
        theme(legend.title = element_text(size = 12, face= "bold")) +
        theme(axis.title.x = element_text(color="black", size=12, face="bold"), axis.title.y = element_text(color="black", size=12, face="bold")) + 
        theme(axis.text.x = element_text(color = "black", size = 12), axis.text.y = element_text(color = "black", size = 12)) +
        ggtitle("Dairy Farms: Weighted UniFrac")
)

#Function
phyloseq_to_df <- function(physeq, addtax = T, addtot = F, addmaxrank = F, sorting = "abundance"){
  
  # require(phyloseq)
  
  ## Data validation
  if(any(addtax == TRUE || sorting == "taxonomy")){
    if(is.null(phyloseq::tax_table(physeq, errorIfNULL = F))){
      stop("Error: taxonomy table slot is empty in the input data.\n")
    }
  }
  
  ## Prepare data frame
  if(taxa_are_rows(physeq) == TRUE){
    res <- data.frame(OTU = phyloseq::taxa_names(physeq), phyloseq::otu_table(physeq), stringsAsFactors = F)
  } else {
    res <- data.frame(OTU = phyloseq::taxa_names(physeq), t(phyloseq::otu_table(physeq)), stringsAsFactors = F)
  }
  
  ## Check if the sample names were silently corrected in the data.frame
  if(any(!phyloseq::sample_names(physeq) %in% colnames(res)[-1])){
    if(addtax == FALSE){
      warning("Warning: Sample names were converted to the syntactically valid column names in data.frame. See 'make.names'.\n")
    }
    
    if(addtax == TRUE){
      stop("Error: Sample names in 'physeq' could not be automatically converted to the syntactically valid column names in data.frame (see 'make.names'). Consider renaming with 'sample_names'.\n")
    }
  }
  
  ## Add taxonomy
  if(addtax == TRUE){
    
    ## Extract taxonomy table
    taxx <- as.data.frame(phyloseq::tax_table(physeq), stringsAsFactors = F)
    
    ## Reorder taxonomy table
    taxx <- taxx[match(x = res$OTU, table = rownames(taxx)), ]
    
    ## Add taxonomy table to the data
    res <- cbind(res, taxx)
    
    ## Add max tax rank column
    if(addmaxrank == TRUE){
      
      ## Determine the lowest level of taxonomic classification
      res$LowestTaxRank <- get_max_taxonomic_rank(taxx, return_rank_only = TRUE)
      
      ## Reorder columns (OTU name - Taxonomy - Max Rank - Sample Abundance)
      res <- res[, c("OTU", phyloseq::rank_names(physeq), "LowestTaxRank", phyloseq::sample_names(physeq))]
      
    } else {
      ## Reorder columns (OTU name - Taxonomy - Sample Abundance)
      res <- res[, c("OTU", phyloseq::rank_names(physeq), phyloseq::sample_names(physeq))]
      
    } # end of addmaxrank
  }   # end of addtax
  
  ## Reorder OTUs
  if(!is.null(sorting)){
    
    ## Sort by OTU abundance
    if(sorting == "abundance"){
      otus <- res[, which(colnames(res) %in% phyloseq::sample_names(physeq))]
      res <- res[order(rowSums(otus, na.rm = T), decreasing = T), ]
    }
    
    ## Sort by OTU taxonomy
    if(sorting == "taxonomy"){
      taxtbl <- as.data.frame( phyloseq::tax_table(physeq), stringsAsFactors = F )
      
      ## Reorder by all columns
      taxtbl <- taxtbl[do.call(order, taxtbl), ]
      # taxtbl <- data.table::setorderv(taxtbl, cols = colnames(taxtbl), na.last = T)
      res <- res[match(x = rownames(taxtbl), table = res$OTU), ]
    }
  }
  
  ## Add OTU total abundance
  if(addtot == TRUE){
    res$Total <- rowSums(res[, which(colnames(res) %in% phyloseq::sample_names(physeq))])
  }
  
  rownames(res) <- NULL
  return(res)
}

pairwise.adonis2 <- function(x, data, strata = NULL, nperm=999, ... ) {
  
  #describe parent call function 
  ststri <- ifelse(is.null(strata),'Null',strata)
  fostri <- as.character(x)
  #list to store results
  
  #copy model formula
  x1 <- x
  # extract left hand side of formula
  lhs <- x1[[2]]
  # extract factors on right hand side of formula 
  rhs <- x1[[3]]
  # create model.frame matrix  
  x1[[2]] <- NULL   
  rhs.frame <- model.frame(x1, data, drop.unused.levels = TRUE) 
  
  # create unique pairwise combination of factors 
  co <- combn(unique(as.character(rhs.frame[,1])),2)
  
  # create names vector   
  nameres <- c('parent_call')
  for (elem in 1:ncol(co)){
    nameres <- c(nameres,paste(co[1,elem],co[2,elem],sep='_vs_'))
  }
  #create results list  
  res <- vector(mode="list", length=length(nameres))
  names(res) <- nameres
  
  #add parent call to res 
  res['parent_call'] <- list(paste(fostri[2],fostri[1],fostri[3],', strata =',ststri, ', permutations',nperm ))
  
  
  #start iteration trough pairwise combination of factors  
  for(elem in 1:ncol(co)){
    
    #reduce model elements  
    if(inherits(eval(lhs),'dist')){	
      xred <- as.dist(as.matrix(eval(lhs))[rhs.frame[,1] %in% c(co[1,elem],co[2,elem]),
                                           rhs.frame[,1] %in% c(co[1,elem],co[2,elem])])
    }else{
      xred <- eval(lhs)[rhs.frame[,1] %in% c(co[1,elem],co[2,elem]),]
    }
    
    mdat1 <-  data[rhs.frame[,1] %in% c(co[1,elem],co[2,elem]),] 
    
    # redefine formula
    if(length(rhs) == 1){
      xnew <- as.formula(paste('xred',as.character(rhs),sep='~'))	
    }else{
      xnew <- as.formula(paste('xred' , 
                               paste(rhs[-1],collapse= as.character(rhs[1])),
                               sep='~'))}
    
    #pass new formula to adonis
    if(is.null(strata)){
      ad <- adonis2(xnew,data=mdat1, ... )
    }else{
      perm <- how(nperm = nperm)
      setBlocks(perm) <- with(mdat1, mdat1[,ststri])
      ad <- adonis2(xnew,data=mdat1,permutations = perm, ... )}
    
    res[nameres[elem+1]] <- list(ad[1:5])
  }
  #names(res) <- names  
  class(res) <- c("pwadstrata", "list")
  return(res)
} 





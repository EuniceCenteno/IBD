# for help installing phyloseq, see this website
# https://bioconductor.org/packages/release/bioc/html/phyloseq.html

# to install phyloseq:
# if (!requireNamespace("BiocManager", quietly = TRUE))
# install.packages("BiocManager")
#BiocManager::install("phyloseq")

library(qiime2R)
library(phyloseq)
library(zoo)
library(tidyverse)
library(tidyr) #for separate function
library(naniar)# Ffor replace all function
library(ggplot2)
library(tidyr) #separate function
library(reshape2) #melt function
library(dplyr)
library(data.table)
library(ggpubr)

setwd("~/Desktop/eunice/PhD/IBDProb/") 
metadata <- read.csv("NG_METADATAIBDADLIBCorrect.csv", na.strings = c("","NA"), header=TRUE)
str(metadata)
metadata$TreatmentGroup <- as.factor(metadata$TreatmentGroup)
levels(metadata$TreatmentGroup) <- list("DSS+H2O"="DSS+H2O", "DSS+LbcVec"="DSS+LbcVec","DSS+LbcLapLin"="DSS+LbcLapLin")
rownames(metadata) = metadata$ID
metadata$dayCat = metadata$day
metadata$dayCat = as.factor(metadata$dayCat)

ASVs <- read_qza("table.qza") #551 samples, but it will be 400 after removing Dairy 175, 176
ASV_s <- as.data.frame(ASVs$data)
ASV_table <- as.data.frame(ASVs$data)
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
sum(x$F1DAY0) #34412

rownames(metadata) = metadata$name
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

TaxASV <- merge(tax.final, ASVkey, by.x = 0, by.y = "ASVstring")
row.names(TaxASV) <- TaxASV[,10]
TaxASV = TaxASV[,-c(1,10)]
#write.csv(TaxASV,"taxonomy.csv", row.names = TRUE)

### Creating the Phyloseq Object
OTU.physeq = otu_table(as.matrix(ASV_table), taxa_are_rows=FALSE)
tax.physeq = tax_table(as.matrix(TaxASV))
#meta.physeq = sample_data(meta)
meta.physeq = sample_data(metadata)

#We then merge these into an object of class phyloseq.
physeq_deseq = phyloseq(OTU.physeq, tax.physeq, meta.physeq)
physeq_deseq #  [ 551 taxa and 34 samples ]

colnames(tax_table(physeq_deseq))
## Filter any non-baxteria, chloroplast and mitochondria
physeq_deseq %>%
  subset_taxa(Family != " Mitochondria" & 
                Genus != " Mitochondria" &
                Species != " Mitochondria" &
                Order != " Chloroplast" &
                Family != " Chloroplast" &
                Genus != " Chloroplast" &
                Species != " Chloroplast") -> physeq_deseq
physeq_deseq #[ 549 taxa and 34 samples ] 

#You need to run the phyloseq_to_df function
prunetable<- phyloseq_to_df(physeq_deseq, addtax = T, addtot = F, addmaxrank = F,
                            sorting = "abundance")

## no mitochondria or chloroplast in the data
NewTax <- prunetable[,c(1:9)]
row.names(NewTax) <- NewTax[,1]
NewTax = NewTax[,-c(1)]
NewTax$ASVs = rownames(NewTax)

NewASVTable <- prunetable[-c(2:9)]
row.names(NewASVTable) <- NewASVTable[,1]
NewASVTable = NewASVTable[,-c(1)]
NewASVTable = t(NewASVTable)

### separate samples per day
day0 = subset(metadata, dayCat == "0")
d0_OTU <- NewASVTable[rownames(NewASVTable) %in% rownames(day0),]

day7 = subset(metadata, dayCat == "7")
d7_OTU <- NewASVTable[rownames(NewASVTable) %in% rownames(day7),]

day17 = subset(metadata, dayCat == "17")
d17_OTU <- NewASVTable[rownames(NewASVTable) %in% rownames(day17),]

## day 0
##start with the healthy samples
OTU.physeqd0 = otu_table(as.matrix(d0_OTU), taxa_are_rows=FALSE)
tax.physeqd0 = tax_table(as.matrix(NewTax))
#meta.physeq = sample_data(meta)
meta.physeqd0 = sample_data(day0)
physeq_d0 = phyloseq(OTU.physeqd0, tax.physeqd0, meta.physeqd0)

##Making the tree, we need the phyloseq object 
library("ape") #to create the tree
random_treed0 = rtree(ntaxa(physeq_d0), rooted=TRUE, tip.label=taxa_names(physeq_d0))
##Merging the tree with the phyloseq object
physeqd0 = merge_phyloseq(physeq_d0, meta.physeqd0, random_treed0)
physeqd0 #This command should show the otu table, sample data, tax table and the phy tree
physeq_bar_plot = physeqd0

str(day0)
my_level <- c("Genus")
my_column <- c("TreatmentGroup")  #this is the metadata column that we will use in the taxa barplot

rm(taxa.summary)

abund_filter <- 0.01  # Our abundance threshold
#ml ="Genus"

for(ml in my_level){
  print(ml)
  
  taxa.summary <- physeq_bar_plot %>%
    tax_glom(taxrank = ml, NArm = FALSE) %>%  # agglomerate at `ml` level
    transform_sample_counts(function(x) {x/sum(x)} ) %>% # Transform to rel. abundance
    psmelt()  %>%                               # Melt to long format
    group_by(get(my_column), get(ml)) %>%
    summarise(Abundance.average=mean(Abundance)) 
  taxa.summary <- as.data.frame(taxa.summary)
  colnames(taxa.summary)[1] <- my_column
  colnames(taxa.summary)[2] <- ml
  
  physeq.taxa.max <- taxa.summary %>% 
    group_by(get(ml)) %>%
    summarise(overall.max=max(Abundance.average))
  
  physeq.taxa.max <- as.data.frame(physeq.taxa.max)
  colnames(physeq.taxa.max)[1] <- ml
  
  # merging the phyla means with the metadata #
  physeq_meta <- merge(taxa.summary, physeq.taxa.max)
  
  
  physeq_meta_filtered <- filter(physeq_meta, overall.max>abund_filter)
}
write.csv(physeq_meta_filtered, "physeq_meta_filteredGenusDay0.csv")

## day 7
##start with the healthy samples
OTU.physeqd7 = otu_table(as.matrix(d7_OTU), taxa_are_rows=FALSE)
tax.physeqd7 = tax_table(as.matrix(NewTax))
#meta.physeq = sample_data(meta)
meta.physeqd7 = sample_data(day7)
physeq_d7 = phyloseq(OTU.physeqd7, tax.physeqd7, meta.physeqd7)

##Making the tree, we need the phyloseq object 
random_treed7 = rtree(ntaxa(physeq_d7), rooted=TRUE, tip.label=taxa_names(physeq_d7))
##Merging the tree with the phyloseq object
physeqd7 = merge_phyloseq(physeq_d7, meta.physeqd7, random_treed7)
physeqd7 #This command should show the otu table, sample data, tax table and the phy tree
physeq_bar_plot7 = physeqd7

str(day7)
my_level <- c("Genus")
my_column <- c("TreatmentGroup")  #this is the metadata column that we will use in the taxa barplot

rm(taxa.summary)

abund_filter <- 0.01  # Our abundance threshold
#ml ="Genus"

for(ml in my_level){
  print(ml)
  
  taxa.summary <- physeq_bar_plot7 %>%
    tax_glom(taxrank = ml, NArm = FALSE) %>%  # agglomerate at `ml` level
    transform_sample_counts(function(x) {x/sum(x)} ) %>% # Transform to rel. abundance
    psmelt()  %>%                               # Melt to long format
    group_by(get(my_column), get(ml)) %>%
    summarise(Abundance.average=mean(Abundance)) 
  taxa.summary <- as.data.frame(taxa.summary)
  colnames(taxa.summary)[1] <- my_column
  colnames(taxa.summary)[2] <- ml
  
  physeq.taxa.max <- taxa.summary %>% 
    group_by(get(ml)) %>%
    summarise(overall.max=max(Abundance.average))
  
  physeq.taxa.max <- as.data.frame(physeq.taxa.max)
  colnames(physeq.taxa.max)[1] <- ml
  
  # merging the phyla means with the metadata #
  physeq_meta <- merge(taxa.summary, physeq.taxa.max)
  
  
  physeq_meta_filtered <- filter(physeq_meta, overall.max>abund_filter)
}
write.csv(physeq_meta_filtered, "physeq_meta_filteredGenusDay7.csv")

## day 17
##start with the healthy samples
OTU.physeqd17 = otu_table(as.matrix(d17_OTU), taxa_are_rows=FALSE)
tax.physeqd17 = tax_table(as.matrix(NewTax))
#meta.physeq = sample_data(meta)
meta.physeqd17 = sample_data(day17)
physeq_d17 = phyloseq(OTU.physeqd17, tax.physeqd17, meta.physeqd17)

##Making the tree, we need the phyloseq object 
random_treed17 = rtree(ntaxa(physeq_d17), rooted=TRUE, tip.label=taxa_names(physeq_d17))
##Merging the tree with the phyloseq object
physeqd17 = merge_phyloseq(physeq_d17, meta.physeqd17, random_treed17)
physeqd17 #This command should show the otu table, sample data, tax table and the phy tree
physeq_bar_plot17 = physeqd17

str(day17)
my_level <- c("Genus")
my_column <- c("TreatmentGroup")  #this is the metadata column that we will use in the taxa barplot

rm(taxa.summary)

abund_filter <- 0.01  # Our abundance threshold
#ml ="Genus"

for(ml in my_level){
  print(ml)
  
  taxa.summary <- physeq_bar_plot17 %>%
    tax_glom(taxrank = ml, NArm = FALSE) %>%  # agglomerate at `ml` level
    transform_sample_counts(function(x) {x/sum(x)} ) %>% # Transform to rel. abundance
    psmelt()  %>%                               # Melt to long format
    group_by(get(my_column), get(ml)) %>%
    summarise(Abundance.average=mean(Abundance)) 
  taxa.summary <- as.data.frame(taxa.summary)
  colnames(taxa.summary)[1] <- my_column
  colnames(taxa.summary)[2] <- ml
  
  physeq.taxa.max <- taxa.summary %>% 
    group_by(get(ml)) %>%
    summarise(overall.max=max(Abundance.average))
  
  physeq.taxa.max <- as.data.frame(physeq.taxa.max)
  colnames(physeq.taxa.max)[1] <- ml
  
  # merging the phyla means with the metadata #
  physeq_meta <- merge(taxa.summary, physeq.taxa.max)
  
  
  physeq_meta_filtered <- filter(physeq_meta, overall.max>abund_filter)
}
write.csv(physeq_meta_filtered, "physeq_meta_filteredGenusDay17.csv")

genus <- read.csv("GenusDayTrtIBDProb.csv", na.strings = c("","NA"), header=TRUE) #505
str(genus)
genus$TreatmentGroup = as.factor(genus$TreatmentGroup)
levels(genus$TreatmentGroup) <- list("DSS+H2O"="DSS+H2O", "DSS+LbcVec"="DSS+LbcVec","DSS+LbcLapLin"="DSS+LbcLapLin")
metadata$day = metadata$day

my_colors <- c(
  '#a6cee3','#1f78b4','#b2df8a','#33a02c','#fb9a99','#e31a1c',
  '#fdbf6f','#ff7f00','#cab2d6','#6a3d9a','#ffff99','#b15928', 
  "#CBD588", "#92C5DE","#DA5724", "#508578", "#CD9BCD",
  "#AD6F3B", "#673770","#652926", "#C84248", 
  "#8569D5", "#5E738F","#D1A33D", "#8A7C64", "#599861", "gray", "black"
) 

ggplot(genus, aes(x = TreatmentGroup, y = Abundance.average, fill =Genus)) + 
  geom_bar(stat = "identity") +
  theme_bw()+
  facet_wrap(day~.) +
  scale_fill_manual(values = my_colors) +
  #ylim(c(0,0.18)) +
  #guides(fill="none") +
  guides(fill = guide_legend(reverse = F, keywidth = 1.5, keyheight = .9, ncol = 1)) +
  theme(legend.text=element_text(size=10, face="italic")) +
  theme(strip.text = element_text(size = 11, face = "bold")) +
  theme(legend.text = element_text(size=12)) +
  theme(legend.title = element_text(size = 13, face= "bold")) +
  theme(legend.key.size = unit(10, "point")) +
  theme(axis.title.x = element_text(color="black", size=11, face="bold"), axis.title.y = element_text(color="black", size=11, face="bold")) + 
  theme(axis.text.x = element_text(color = "black", size = 10), axis.text.y = element_text(color = "black", size = 11)) +
  ggtitle("Genus >1% per sample") +
  ylab(paste0("Relative Abundance")) +  labs(x='Treatment')

family <- read.csv("FamilyDayTRTIBDProb.csv", na.strings = c("","NA"), header=TRUE) #505
str(family)
family$TreatmentGroup = as.factor(family$TreatmentGroup)
levels(family$TreatmentGroup) <- list("DSS+H2O"="DSS+H2O", "DSS+LbcVec"="DSS+LbcVec","DSS+LbcLapLin"="DSS+LbcLapLin")
family$day = family$day

ggplot(family, aes(x = TreatmentGroup, y = Abundance.average, fill =Family)) + 
  geom_bar(stat = "identity") +
  theme_bw()+
  facet_wrap(day~.) +
  scale_fill_manual(values = my_colors) +
  #ylim(c(0,0.18)) +
  #guides(fill="none") +
  guides(fill = guide_legend(reverse = F, keywidth = 1.5, keyheight = .9, ncol = 1)) +
  theme(legend.text=element_text(size=10, face="italic")) +
  theme(strip.text = element_text(size = 11, face = "bold")) +
  theme(legend.text = element_text(size=12)) +
  theme(legend.title = element_text(size = 13, face= "bold")) +
  theme(legend.key.size = unit(10, "point")) +
  theme(axis.title.x = element_text(color="black", size=11, face="bold"), axis.title.y = element_text(color="black", size=11, face="bold")) + 
  theme(axis.text.x = element_text(color = "black", size = 10), axis.text.y = element_text(color = "black", size = 11)) +
  ggtitle("Family >1% per sample") +
  ylab(paste0("Relative Abundance")) +  labs(x='Treatment')


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

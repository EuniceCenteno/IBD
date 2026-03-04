###Ancom
#combine the output of Qiime with the ASV metadata
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

##### subsetting data
##############################################
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
sum(x$F1DAY0) #34412

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
NewTax$ASV = rownames(NewTax)
write.csv(NewTax, "TaxonomyIBDProb.csv")
#write.table(NewTax,"NewTax.txt",sep=",", row.names = TRUE)
NewASVTable <- prunetable
NewASVTable <- NewASVTable[,-c(2:9)]
row.names(NewASVTable) <- NewASVTable[,1]
NewASVTable = NewASVTable[,-c(1)]

##### separation by dat
NewASVTable2 = t(NewASVTable)
day0 = subset(metadata, dayCat == "0")
d0_OTU <- NewASVTable2[rownames(NewASVTable2) %in% rownames(day0),]
d0_OTU = t(d0_OTU)

day7 = subset(metadata, dayCat == "7")
d7_OTU <- NewASVTable2[rownames(NewASVTable2) %in% rownames(day7),]
d7_OTU = t(d7_OTU)

day17 = subset(metadata, dayCat == "17")
d17_OTU <- NewASVTable2[rownames(NewASVTable2) %in% rownames(day17),]
d17_OTU = t(d17_OTU)

write.csv(day0 , "IBDProbday0.csv")
write.csv(day7 , "IBDProbday7.csv")
write.csv(day17 , "IBDProbday17.csv")

write.table(d0_OTU, file='d0_OTUTable.tsv', quote=FALSE, sep='\t')
write.table(d7_OTU, file='d7_OTUTable.tsv', quote=FALSE, sep='\t')
write.table(d17_OTU, file='d17_OTUTable.tsv', quote=FALSE, sep='\t')



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

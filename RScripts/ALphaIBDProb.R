library(afex)
library(lme4)
library(emmeans)
library(lubridate)
library(ggplot2)
library("cowplot")
theme_set(theme_grey())
#install.packages("sjstats")
library(jtools)
library(ggpubr)
library(sjstats)
library(ggdist)
library(tidyquant)
library(dplyr)
library(tidyverse)
#install.packages("wesanderson")
library(wesanderson)
library(qiime2R)

setwd("~/Desktop/eunice/PhD/IBDProb/") 
metadata <- read.csv("NG_METADATAIBDADLIBCorrect.csv", na.strings = c("","NA"), header=TRUE)
str(metadata)
metadata$TreatmentGroup <- as.factor(metadata$TreatmentGroup)
otu_table <- read.table("./AlphaIBD/exported/observed_asvs.tsv", header=TRUE, row.names=1, sep="\t")
shannon <-read.table("./AlphaIBD/exported/shannon.tsv", header=TRUE, row.names=1, sep="\t")
faith <- read.table("./AlphaIBD/exported/faith_pd.tsv", header=TRUE, row.names=1, sep="\t")
evenness <- read.table("./AlphaIBD/exported/evenness.tsv", header=TRUE, row.names=1, sep="\t")
chao1 <- read.table("./AlphaIBD/exported/chao1.tsv", header=TRUE, row.names=1, sep="\t")
alpha_diversity <- merge(otu_table, shannon, by.x = 0, by.y = 0)
alpha_diversity <- merge(alpha_diversity, faith, by.x = "Row.names", by.y = 0)
alpha_diversity <- merge(alpha_diversity, evenness, by.x = "Row.names", by.y = 0)
alpha_diversity <- merge(alpha_diversity, chao1, by.x = "Row.names", by.y = 0)
colnames(alpha_diversity) <- c("Row.names", "observed_otus", "shannon", "faith_pd", "pielou_e", "chao1")
metadata <- merge(metadata, alpha_diversity, by.x = "ID", by.y = "Row.names")
row.names(metadata) <- metadata$ID

## 
## check if the data is normal distributed or not
set_sum_contrasts() # important for afex

# full model
str(metadata)
summary(metadata)
metadata$dayCat = metadata$day
metadata$dayCat = as.factor(metadata$dayCat)

M1 <- mixed(pielou_e ~ dayCat +TreatmentGroup + dayCat*TreatmentGroup+  (1|animal), data = metadata,method = "KR", 
            control = lmerControl(optCtrl = list(maxfun = 1e6)), expand_re = TRUE)
anova(M1)
plot(M1$full_model)
shapiro.test(metadata$pielou_e) 
qqnorm(residuals(M1$full_model))
qqline(residuals(M1$full_model))

M2 <- mixed(observed_otus ~ dayCat +TreatmentGroup+ dayCat*TreatmentGroup +  (1|animal), data = metadata,method = "KR", 
            control = lmerControl(optCtrl = list(maxfun = 1e6)), expand_re = TRUE)
anova(M2)
plot(M2$full_model)
shapiro.test(metadata$observed_otus)
qqnorm(residuals(M2$full_model))
qqline(residuals(M2$full_model))

M3 <- mixed(chao1 ~ dayCat +TreatmentGroup +dayCat*TreatmentGroup+  (1|animal), data = metadata,method = "KR", 
            control = lmerControl(optCtrl = list(maxfun = 1e6)), expand_re = TRUE)
anova(M3)
plot(M3$full_model)
shapiro.test(metadata$chao1)
qqnorm(residuals(M3$full_model))
qqline(residuals(M3$full_model))

M4 <- mixed(faith_pd ~ dayCat +TreatmentGroup+dayCat*TreatmentGroup +  (1|animal), data = metadata,method = "KR", 
            control = lmerControl(optCtrl = list(maxfun = 1e6)), expand_re = TRUE)
anova(M4)
plot(M4$full_model)
shapiro.test(metadata$faith_pd)
qqnorm(residuals(M4$full_model))
qqline(residuals(M4$full_model))

M5 <- mixed(shannon ~ dayCat +TreatmentGroup+dayCat*TreatmentGroup +  (1|animal), data = metadata,method = "KR", 
            control = lmerControl(optCtrl = list(maxfun = 1e6)), expand_re = TRUE)
anova(M5)
plot(M5$full_model)
shapiro.test(metadata$shannon)
qqnorm(residuals(M5$full_model))
qqline(residuals(M5$full_model))

## test post-hoc
str(metadata)
emm_options(lmer.df = "asymptotic") # also possible: 'satterthwaite', 'kenward-roger'
emm_fsize <- emmeans(M1, "dayCat")
emm_fsize
pie_order <- data.frame(emm_fsize)
pairs(emm_fsize)
update(pairs(emm_fsize), by = NULL, adjust = "holm")

str(metadata)
emm_options(lmer.df = "asymptotic") # also possible: 'satterthwaite', 'kenward-roger'
emm_Obs <- emmeans(M2, "dayCat")
emm_Obs
obs_order <- data.frame(emm_Obs)
pairs(emm_Obs)
update(pairs(emm_Obs), by = NULL, adjust = "holm")

str(metadata)
my_color1 <- c("black", "red", "green")
levels(metadata$TreatmentGroup) <- list("DSS+H2O"="DSS+H2O", "DSS+LbcVec"="DSS+LbcVec","DSS+LbcLapLin"="DSS+LbcLapLin")

a = ggplot(metadata, aes(x=dayCat, y=pielou_e)) +
  geom_jitter(aes(x=dayCat, y=pielou_e, color=TreatmentGroup), width = 0.25, alpha=0.7) +
  geom_errorbar(data=pie_order, aes(x=dayCat, ymin=emmean - SE,
                                ymax=emmean + SE, y=NULL), color="black", width=0.2) +
  theme_bw()+ scale_color_manual(values = my_color1) +
  geom_point(data= pie_order, aes(y = emmean, x = dayCat)) +
  ylab("Evenness") +xlab ("Day") + #guides(color="none") +
  #scale_color_manual(values = my_color1) + ylim(0,3750) +
  #stat_compare_means(comparisons = my_comparisons, label = "p.signif")+
  theme(strip.text = element_text(size = 12, face = "bold")) +
  theme(legend.text = element_text(size=12)) +
  theme(legend.title = element_text(size = 12, face= "bold")) +
  theme(legend.key.size = unit(12, "point")) +
  theme(axis.title.x = element_text(color="black", size=12, face="bold"),
        axis.title.y = element_text(color="black", size=12, face="bold")) +
  theme(axis.text.x = element_text(color = "black", size = 12), axis.text.y =
          element_text(color = "black", size = 12))

b = ggplot(metadata, aes(x=dayCat, y=observed_otus)) +
  geom_jitter(aes(x=dayCat, y=observed_otus, color=TreatmentGroup), width = 0.25, alpha=0.5) +
  geom_errorbar(data=obs_order, aes(x=dayCat, ymin=emmean - SE,
                                    ymax=emmean + SE, y=NULL), color="black", width=0.2) +
  theme_bw()+ scale_color_manual(values = my_color1) +
  geom_point(data= obs_order, aes(y = emmean, x = dayCat)) +
  ylab("Observed ASVs") +xlab ("Day") + #guides(color="none") +
  #scale_color_manual(values = my_color1) + ylim(0,3750) +
  #stat_compare_means(comparisons = my_comparisons, label = "p.signif")+
  theme(strip.text = element_text(size = 12, face = "bold")) +
  theme(legend.text = element_text(size=12)) +
  theme(legend.title = element_text(size = 12, face= "bold")) +
  theme(legend.key.size = unit(12, "point")) +
  theme(axis.title.x = element_text(color="black", size=12, face="bold"),
        axis.title.y = element_text(color="black", size=12, face="bold")) +
  theme(axis.text.x = element_text(color = "black", size = 12), axis.text.y =
          element_text(color = "black", size = 12))

str(metadata)
emm_options(lmer.df = "asymptotic") # also possible: 'satterthwaite', 'kenward-roger'
emm_chao <- emmeans(M3, "dayCat")
emm_chao
chao_order <- data.frame(emm_chao)
pairs(emm_chao)
update(pairs(emm_chao), by = NULL, adjust = "holm")

c = ggplot(metadata, aes(x=dayCat, y=chao1)) +
  geom_jitter(aes(x=dayCat, y=chao1, color=TreatmentGroup), width = 0.25, alpha=0.5) +
  geom_errorbar(data=chao_order, aes(x=dayCat, ymin=emmean - SE,
                                    ymax=emmean + SE, y=NULL), color="black", width=0.2) +
  theme_bw()+ scale_color_manual(values = my_color1) +
  geom_point(data= chao_order, aes(y = emmean, x = dayCat)) +
  ylab("Chao1") +xlab ("Day") + #guides(color="none") +
  #scale_color_manual(values = my_color1) + ylim(0,3750) +
  #stat_compare_means(comparisons = my_comparisons, label = "p.signif")+
  theme(strip.text = element_text(size = 12, face = "bold")) +
  theme(legend.text = element_text(size=12)) +
  theme(legend.title = element_text(size = 12, face= "bold")) +
  theme(legend.key.size = unit(12, "point")) +
  theme(axis.title.x = element_text(color="black", size=12, face="bold"),
        axis.title.y = element_text(color="black", size=12, face="bold")) +
  theme(axis.text.x = element_text(color = "black", size = 12), axis.text.y =
          element_text(color = "black", size = 12))

emm_fai <- emmeans(M4, "dayCat")
emm_fai
fai_order <- data.frame(emm_fai)
pairs(emm_fai)
update(pairs(emm_fai), by = NULL, adjust = "holm")

d = ggplot(metadata, aes(x=dayCat, y=faith_pd)) +
  geom_jitter(aes(x=dayCat, y=faith_pd, color=TreatmentGroup), width = 0.25, alpha=0.5) +
  geom_errorbar(data=fai_order, aes(x=dayCat, ymin=emmean - SE,
                                     ymax=emmean + SE, y=NULL), color="black", width=0.2) +
  theme_bw()+ scale_color_manual(values = my_color1) +
  geom_point(data= fai_order, aes(y = emmean, x = dayCat)) +
  ylab("Phylogenetic Diversity") +xlab ("Day") + #guides(color="none") +
  #scale_color_manual(values = my_color1) + ylim(0,3750) +
  #stat_compare_means(comparisons = my_comparisons, label = "p.signif")+
  theme(strip.text = element_text(size = 12, face = "bold")) +
  theme(legend.text = element_text(size=12)) +
  theme(legend.title = element_text(size = 12, face= "bold")) +
  theme(legend.key.size = unit(12, "point")) +
  theme(axis.title.x = element_text(color="black", size=12, face="bold"),
        axis.title.y = element_text(color="black", size=12, face="bold")) +
  theme(axis.text.x = element_text(color = "black", size = 12), axis.text.y =
          element_text(color = "black", size = 12))

emm_sha <- emmeans(M5, "dayCat")
emm_sha
sha_order <- data.frame(emm_sha)
pairs(emm_sha)
update(pairs(emm_sha), by = NULL, adjust = "holm")

e = ggplot(metadata, aes(x=dayCat, y=shannon)) +
  geom_jitter(aes(x=dayCat, y=shannon, color=TreatmentGroup), width = 0.25, alpha=0.5) +
  geom_errorbar(data=sha_order, aes(x=dayCat, ymin=emmean - SE,
                                    ymax=emmean + SE, y=NULL), color="black", width=0.2) +
  theme_bw()+ scale_color_manual(values = my_color1) +
  geom_point(data= sha_order, aes(y = emmean, x = dayCat)) +
  ylab("Shannon") +xlab ("Day") + #guides(color="none") +
  #scale_color_manual(values = my_color1) + ylim(0,3750) +
  #stat_compare_means(comparisons = my_comparisons, label = "p.signif")+
  theme(strip.text = element_text(size = 12, face = "bold")) +
  theme(legend.text = element_text(size=12)) +
  theme(legend.title = element_text(size = 12, face= "bold")) +
  theme(legend.key.size = unit(12, "point")) +
  theme(axis.title.x = element_text(color="black", size=12, face="bold"),
        axis.title.y = element_text(color="black", size=12, face="bold")) +
  theme(axis.text.x = element_text(color = "black", size = 12), axis.text.y =
          element_text(color = "black", size = 12))

ggarrange(b,c,a,d,e, nrow=3, ncol = 2, common.legend = FALSE)

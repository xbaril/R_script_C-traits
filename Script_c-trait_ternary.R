

###################################################### 16S ######################################################
remove(list=ls())

########### ANCOM-BC
library(phyloseq)
library(ANCOMBC)


setwd ("...") # Set Working Directory


###### 16S ########

otumat1 = read.csv("ASV_table_16S_1percent.csv", row.names=1) #### ASV table for 1 percent carbon amendment
head(otumat1)
str(otumat1)

taxmat1 = read.csv("Taxonomy_table_16S.csv", row.names=1) #### Taxonomy table including a "Species" column 
head(taxmat1)
str(taxmat1)

otumat = data.matrix(otumat1, rownames.force = TRUE)
class(otumat)
#"matrix" "array"#

taxmat = as.matrix(taxmat1, rownames.force = TRUE)
class(taxmat)
#"matrix" "array"#

OTU = otu_table(otumat, taxa_are_rows = TRUE)
TAX = tax_table(taxmat)

# convert data to matrix
physeq = phyloseq(OTU, TAX)

meta = data.frame(
  Treatment = c("Starch", "Starch", "Sucrose", "Sucrose", "Sucrose", "Cellulose", "Cellulose", "Cellulose", "Starch", "A_Control", "A_Control", "A_Control"),
  Bloc =c("B",	"C",	"A",	"B",	"C",	"A",	"B",	"C",	"A",	"A",	"B",	"C"),
  row.names=sample_names(OTU),
  stringsAsFactors=FALSE)
str(meta)


meta1 = sample_data(meta)

sample_names(OTU)
sample_names(meta1)
sample_names(TAX)


physeq <- phyloseq(OTU, TAX, meta1)
physeq

# https://www.bioconductor.org/packages/release/bioc/html/ANCOMBC.html#



library(microbiome)

# The taxonomy table
tax_mat = as(tax_table(physeq), "matrix")
tax_mat
# Run ancombc function
#max.print=99999999


out <- ancombc2(
  physeq,
  assay_name = "counts",
  tax_level = "Species", ### Species refer to ASV
  fix_formula = "Treatment",
  rand_formula = NULL,
  p_adj_method = "holm",
  pseudo = 0,
  pseudo_sens = TRUE,
  prv_cut = 0.1,
  lib_cut = 0,
  s0_perc = 0.05,
  group = NULL,
  struc_zero = FALSE,
  neg_lb = FALSE,
  alpha = 0.05,
  n_cl = 1,
  verbose = FALSE,
  global = FALSE,
  pairwise = FALSE,
  dunnet = FALSE,
  trend = FALSE,
  iter_control = list(tol = 0.01, max_iter = 20, verbose = FALSE),
  em_control = list(tol = 1e-05, max_iter = 100),
  mdfdr_control = list(fwer_ctrl_method = "holm", B = 100),
  trend_control = list(contrast = NULL, node = NULL, solver = "ECOS", B = 100)
)


#### Show results

head(out$feature_table)
res <-out$res
head(res)
write.csv(res, file="res_16S_ANCOMBC2.csv")


###################################### 


#1 remove lines (ASV) with less than 3 non-zero observations out of 9 (Su, St & C treatments)  --------

ASV_16S <- read.csv("ASV_table_16S_1percent.csv")
head(ASV_16S)
nrow(ASV_16S)
#expression df > 0 produces a matrix with values TRUE where element > 0
#Function rowSums returns number of nonzero elements for every line (when summing it treats values TRUE as 1 and FALSE as 0)
#finally [selects only lines where the number of non-zero elements >= 3
ASV_16S_2 <- ASV_16S[,-c(11,12,13)]
ASV_16S_no <- ASV_16S_2[ rowSums(ASV_16S_2 > 0) >= 3, ]
head(ASV_16S_no)
nrow(ASV_16S_no)


#2 Kruskal Wallis test between treatments to ensure that library size is not significantly different  --------
## build data fram with library size

size_Su <- c(sum(ASV_16S_no$X16S_1_Sucrose_A),sum(ASV_16S_no$X16S_1_Sucrose_B), sum(ASV_16S_no$X16S_1_Sucrose_C))
size_C <- c(sum(ASV_16S_no$X16S_1_Cellulose_A),sum(ASV_16S_no$X16S_1_Cellulose_B), sum(ASV_16S_no$X16S_1_Cellulose_C))
size_St <- c(sum(ASV_16S_no$X16S_1_Starch_A),sum(ASV_16S_no$X16S_1_Starch_B), sum(ASV_16S_no$X16S_1_Starch_C))

df_size <- data.frame(size_Su, size_C, size_St)
df_size

# K-W rank sum test : kruskal.test {stats}
kruskal.test(df_size)
#data:  df_size




#3 Create the compositional matrix (log - ratio) ----------
### decostand
library(vegan)
t_deco_16S <- decostand(t(ASV_16S_no[,-1]), method = "clr", pseudocount = 1)
head(t_deco_16S) ## you can validate in excel that the function is performed on the samples
sum(t_deco_16S[1,])
sum(t_deco_16S[,1])

## resetting the matrix with samples in columns
deco_16S <- t(t_deco_16S)
head(deco_16S)

## matrice -> data frame
df_16S_clr <- as.data.frame(deco_16S)

## make the matrix positive: add the absolute of the minimum value of each column to all the cells in that column (sample)
summary(df_16S_clr)

starch_A_16S <- (df_16S_clr$X16S_1_Starch_A + abs(min(df_16S_clr$X16S_1_Starch_A)))
starch_B_16S <- (df_16S_clr$X16S_1_Starch_B + abs(min(df_16S_clr$X16S_1_Starch_B)))
starch_C_16S <- (df_16S_clr$X16S_1_Starch_C + abs(min(df_16S_clr$X16S_1_Starch_C)))
sucrose_A_16S <- (df_16S_clr$X16S_1_Sucrose_A + abs(min(df_16S_clr$X16S_1_Sucrose_A)))
sucrose_B_16S <- (df_16S_clr$X16S_1_Sucrose_B + abs(min(df_16S_clr$X16S_1_Sucrose_B)))
sucrose_C_16S <- (df_16S_clr$X16S_1_Sucrose_C + abs(min(df_16S_clr$X16S_1_Sucrose_C)))
cellulose_A_16S <- (df_16S_clr$X16S_1_Cellulose_A + abs(min(df_16S_clr$X16S_1_Cellulose_A)))
cellulose_B_16S <- (df_16S_clr$X16S_1_Cellulose_B + abs(min(df_16S_clr$X16S_1_Cellulose_B)))
cellulose_C_16S <- (df_16S_clr$X16S_1_Cellulose_C + abs(min(df_16S_clr$X16S_1_Cellulose_C)))
ASV <- ASV_16S_no$ASV

df_16S_pos <- data.frame(ASV, starch_A_16S, starch_B_16S, starch_C_16S, sucrose_A_16S, sucrose_B_16S, sucrose_C_16S, cellulose_A_16S, cellulose_B_16S, cellulose_C_16S)
head(df_16S_pos)

#4 data frame creation with ASV frequency averaging ---------
ASV_16S_no2 <- df_16S_pos
names(ASV_16S_no2)
# renaming columns
names(ASV_16S_no2) <- c("ASV","starch_A", "starch_B", "starch_C", "sucrose_A", "sucrose_B", "sucrose_C", "cellulose_A", "cellulose_B", "cellulose_C")
names(ASV_16S_no2)

m_Starch <- (ASV_16S_no2$starch_A +ASV_16S_no2$starch_B + ASV_16S_no2$starch_C)/3
m_Cellulose <- (ASV_16S_no2$cellulose_A +ASV_16S_no2$cellulose_B + ASV_16S_no2$cellulose_C)/3
m_Sucrose <- (ASV_16S_no2$sucrose_A +ASV_16S_no2$sucrose_B + ASV_16S_no2$sucrose_C)/3
ASV <- ASV_16S_no2$ASV

m_16S <- data.frame(ASV, m_Starch, m_Cellulose, m_Sucrose)
View(m_16S)
names(m_16S)

#### Corrected abundance ratio on 1. 
##    where Su_ness = m_Sucrose/(m_Sucrose + m_Cellulose + m_Starch)

Su_ness <- m_16S$m_Sucrose/(m_16S$m_Sucrose+m_16S$m_Cellulose+m_16S$m_Starch)
C_ness <- (m_16S$m_Cellulose/(m_16S$m_Sucrose+m_16S$m_Cellulose+m_16S$m_Starch))
St_ness <- (m_16S$m_Starch/(m_16S$m_Sucrose+m_16S$m_Cellulose+m_16S$m_Starch))



df_SCS <- data.frame(ASV, Su_ness, C_ness, St_ness)
View(df_SCS)


nrow(df_SCS)
write.csv(df_SCS, file = "data_ness_16s.csv")



#5 Add info ANCOM control-treatment for the graphic ------------------
# If p.val > 0.05 effect of treatment is noted as NA
# If p.val < 0.05 effect of treatment is noted as Su+, St+, C+, Su-, St-, C-,
ancom_16S<- read.csv("res_16S_ANCOMBC2.csv")
head(ancom_16S)

ancom_16S$effect_Su <- ifelse(ancom_16S$q_TreatmentSucrose <= 0.05,ancom_16S$lfc_TreatmentSucrose, 0) ### determine is the ASV respond significantly to the treatment vs control
ancom_16S$effect_Su_col <- ifelse(ancom_16S$effect_Su < 0, "Su-", ifelse(ancom_16S$effect_Su>0, "Su+", "NA")) ### determine is the ASV response is positive (+) or negative (-)
head(ancom_16S)
ancom_16S$effect_St <- ifelse(ancom_16S$q_TreatmentStarch <= 0.05,ancom_16S$lfc_TreatmentStarch, 0)
ancom_16S$effect_St_col <- ifelse(ancom_16S$effect_St < 0, "St-", ifelse(ancom_16S$effect_St>0, "St+", "NA"))
head(ancom_16S)
ancom_16S$effect_C <- ifelse(ancom_16S$q_TreatmentCellulose <= 0.05,ancom_16S$lfc_TreatmentCellulose, 0)
ancom_16S$effect_C_col <- ifelse(ancom_16S$effect_C < 0, "C-", ifelse(ancom_16S$effect_C>0, "C+", "NA"))
head(ancom_16S)


ancom_16S$Col  <- paste(ancom_16S$effect_Su_col, ancom_16S$effect_St_col, ancom_16S$effect_C_col,  sep="_")
head(ancom_16S)
### data frame with only ASV present in the SCS frame
ASV_list <- df_SCS$ASV
ancom_16S

ancom_16S_1 <- ancom_16S
colnames(ancom_16S_1)[1] <- "ASV"

row.names(ancom_16S_1) <- ancom_16S_1$ASV
head(ancom_16S_1)

rows_to_keep <- c(ASV_list)

ancom_16S_2 <- ancom_16S_1[rows_to_keep,]
head(ancom_16S_2)

df_SCS$effect_col <- ancom_16S_2$Col
head(df_SCS)

###### create data frame with trait and taxonomy
taxo1 <- read.csv("Taxonomy_table_16S.csv", row.names = 1)
head(taxo1)
summary(taxo1)
taxo2 <- taxo1[rows_to_keep,]
head(taxo2)
summary(taxo2)


df_trait_taxo <- taxo2
head(df_SCS)
df_trait_taxo$Su_ness <- df_SCS$Su_ness
df_trait_taxo$C_ness <- df_SCS$C_ness
df_trait_taxo$St_ness <- df_SCS$St_ness
df_trait_taxo$effect_col <- df_SCS$effect_col
head(df_trait_taxo)

write.csv(df_trait_taxo, file = "df_trait_taxo_16s.csv")


#6 Ternary plot ------------------
library(ggrepel)
library(ggplot2)
library(ggtern)


### Reorder Data Frame Rows to show significant ASV up front
library(tidyverse)

data_trait_16S <- read.csv("df_trait_taxo_16s.csv")
head(data_trait_16S)
data_trait_16S_2 <- data_trait_16S %>% arrange(effect_col)### reorder the "effer_col" column in alphabetical order
View(data_trait_16S_2)
data_trait_16S_3 <- data_trait_16S_2 %>% slice(5:1039,1:4) ### reorder manually



### PLot
plot <- ggtern(data = data_trait_16S_3, aes(C_ness, Su_ness, St_ness), label=ASV)+
  geom_point(aes(colour = effect_col, size = 0.5), show.legend = FALSE) +
  scale_color_manual(values = c("Su+_NA_NA" = "#c21717","Su-_NA_NA" = "#f77979",
                                "NA_St+_NA"="#259609", "NA_St-_NA"="#aee0a2", "NA_NA_C+"="#0a23a3", "NA_NA_C-"="#9aa6e3",
                                "Su+_St+_C+"="#000000","Su-_St-_C-"= "#ffffff","NA_NA_NA" = "#adadad99", "Su+_St+_NA"="#5c5a0c",
                                "Su+_NA_C+"="#a91db3", "NA_St+_C+"="#0b8185", "Su-_St-_NA"="#ebe646"))+ ### Colors can be added for each condition
  theme_rgbw() +
  guides(color = "none", fill = "none", alpha = "none")

plot

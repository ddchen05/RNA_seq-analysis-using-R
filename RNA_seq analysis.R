
# script to manipulate gene expression data
# install libraries if needed
# Install requested libraries 
#install.packages('BiocManager')
#BiocManager::install('limma')
#BiocManager::install('DESeq2')
#BiocManager::install('edgeR')


# load libraries
library(dplyr)
library(tidyverse)
library(GEOquery)
library(tidyverse)
library(ggplot2)

# read in the data ---------
data_1 <- read.csv(file = "data.csv")
dim(data_1)


# get metadata --------
metadata1 <- getGEO(GEO = 'GSE183947', GSEMatrix = TRUE)

metadata <- pData(phenoData(metadata1[[1]]))
head(metadata)

# select, mutate, rename the metadata
metadata.m <- metadata %>%
  select(1,10,11,17) %>%
  rename(tissue = characteristics_ch1) %>%
  rename(metastasis = characteristics_ch1.1) %>%
  mutate(tissue = gsub("tissue: ", "", tissue)) %>%
  mutate(metastasis = gsub("metastasis: ", "", metastasis))


# looking at gene expression data ---------
head(data_1)

# reshaping data - from wide to long--------
data.long <- data_1 %>%
  rename(gene = X) %>%
  gather(key = 'samples', value = 'FPKM', -gene)


# join dataframes = data.long + metadata.m
data.combined <- data.long %>%
  left_join(., metadata.m, by = c("samples" = "description")) 

# visulize specific gene expression
# 1. bar plots
data.combined %>%
  filter(gene == 'TP53') %>%
  ggplot(., aes(x = samples, y = FPKM, fill = tissue),fig(18,4)) +
  geom_col()+ggtitle(" TP53_expression_1")
  ggsave('TP53_expression_1.png', device = "png")

# 2. boxplot
data.combined %>%
  filter(gene == 'TP53') %>%
  ggplot(., aes(x = tissue, y = FPKM)) +
  geom_boxplot()+ggtitle(" TP53_expression_2")
  ggsave('TP53_expression_2.png', device = "png")
    
# 3. heatmap
genes.of.interest <- c('BRCA1', 'BRCA2','PTEN','PALB2','CDH1', 'TP53', 'ATM', 'CHEK2')
data.combined %>%
  filter(gene %in% genes.of.interest)%>%
  ggplot(., aes(x = tissue, y = gene, fill = FPKM),fig(20,4)) +
  geom_tile() +
  scale_fill_gradient(low = 'white', high = 'red')+ggtitle(" differiential gene expression")
  ggsave('interesting genes.png', device = "png")
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
# Supp. tables
# Author: Chao Cheng; Lamis Naddaf
# 04/03/2026
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


Dir="./Supp_Tables/"
myinf1 = "./data/Zhu_GO30140_IMbrave150_atezo_bev.rda"
myinf2 = "./data/Zhu_GO30140_IMbrave150_atezo_bev__GenomicEvent_iRAS.txt"

data <- read.table(myinf2, header=T, sep="\t", row.names=1, quote="")
cnum = ncol(data)/4
data = data[, 1:cnum] - data[, (cnum+1):(2*cnum)]
tmp = colnames(data)
tmp = gsub("\\.up\\.ES", "", tmp)
colnames(data) = tmp
se = grep("mul.adj__", colnames(data))
data = data[, se]
for(i in 1:ncol(data))
{
  data[,i] = data[,i]/sd(abs(data[,i]))
}
colnames(data) = gsub("mul.adj__", "", colnames(data))


load(file= myinf1)
info = info	##  samples

se = which(info$Visit =="Pre-treatment")
info = info[se,]
se = which(!info$Confirmed.Response_IRF %in% c("", "NE"))
info = info[se, ]

comxx = intersect(row.names(data), row.names(info))
data = data[comxx,]
info = info[comxx,]
rna<-rna[,comxx]

score = as.numeric(data[, "TP53__MUT"])
info$score = score

#-----------------------------
# TP53 high subset
#-----------------------------
info$TP53_group <- ifelse(info$score > median(info$score), "High", "Low")
info_high <- info[info$TP53_group=="High", ]
rna_high <- rna[, rownames(info_high)]

#-----------------------------
# Create CR/PR vs PD group
#-----------------------------
info_high$ResponseGroup <- ifelse(info_high$Confirmed.Response_IRF %in% c("CR","PR"), "CR_PR",
                                  ifelse(info_high$Confirmed.Response_IRF=="PD", "PD", NA))
info_high <- info_high[!is.na(info_high$ResponseGroup), ]
rna_high <- rna_high[, rownames(info_high)]
rna_high<-log2(rna_high+1)
#-----------------------------
# DEG analysis
#-----------------------------
table(info_high$ResponseGroup)
index1<-which(info_high$ResponseGroup=="CR_PR")
dat_CR_PR<-rna_high[,index1]

index1<-which(info_high$ResponseGroup=="PD")
dat_PD<-rna_high[,index1]

# myavg1 = apply(dat_CR_PR, 1, mean)
# myavg2 = apply(dat_PD, 1, mean)
# myvar1 = apply(dat_CR_PR, 1, var)
# myvar2 = apply(dat_PD, 1, var)
# n1 = ncol(dat_CR_PR)
# n2 = ncol(dat_PD)
# tscore = (myavg1-myavg2)/sqrt(myvar1/n1+myvar2/n2)
# df= (myvar1/n1+myvar2/n2)^2/((myvar1/n1)^2/(n1-1) + (myvar2/n2)^2/(n2-1))
# pval = pt(-abs(tscore), df)*2
# 


dge <- DGEList(counts = rna_high)
dge <- calcNormFactors(dge)

design <- model.matrix(~ ResponseGroup, data = info_high)
v <- voom(dge, design, plot = TRUE)

fit <- lmFit(v, design)
fit <- eBayes(fit)

# coef=2 corresponds to ResponseGroup PD vs CR_PR
res <- topTable(fit, coef=2, number = Inf, adjust.method = "BH")
head(res)

info_high$ResponseGroup <- factor(
  info_high$ResponseGroup,
  levels = c("CR_PR", "PD")
)

design <- model.matrix(~ 0 + ResponseGroup, data = info_high)
colnames(design) <- c("CR_PR", "PD")

contrast <- makeContrasts(
  CRPR_vs_PD = CR_PR - PD,
  levels = design
)

v <- voom(dge, design, plot = TRUE)
fit <- lmFit(v, design)
fit <- contrasts.fit(fit, contrast)
fit <- eBayes(fit)

res <- topTable(fit, number = Inf)


#-----------------------------
# Save results
#-----------------------------
saveRDS(res, paste0(Dir,"High_P53_CRPR_vs_PD_DE.rds"))




#+library(ggplot2)
library(clusterProfiler)
library(org.Hs.eg.db)
library(EnhancedVolcano)
library(dplyr)
myinf1 = "./data/Zhu_GO30140_IMbrave150_atezo_bev.rda"
load(file= myinf1)
##### enriched pathways

res<- readRDS(paste0(Dir,"High_P53_CRPR_vs_PD_DE.rds")) 
write.csv(res,paste0(Dir,"Table_S3_CR_PR_vs_PD_DEGs.csv"))

library(clusterProfiler)
library(org.Hs.eg.db)  # human gene annotation
library(dplyr)

# Assume your differential expression results are in `res` with columns: logFC, P.Value, gene
# Create a named vector of logFC ranked by effect size
gene_list <- res %>%
  arrange(desc(logFC)) %>%   # sort from most up to most down
  pull(logFC)                # extract logFC

names(gene_list) <- rownames(res)

# Make sure there are no duplicated gene names
gene_list <- sort(gene_list, decreasing = TRUE)
gene_list <- gene_list[!duplicated(names(gene_list))]


gsea_res <- gseGO(
  geneList = gene_list,
  OrgDb = org.Hs.eg.db,
  ont = "BP",            # Biological Process
  keyType = "SYMBOL",    # match your gene names
  nPerm = 1000,
  minGSSize = 10,
  maxGSSize = 500,
  pvalueCutoff = 0.05,
  verbose = TRUE
)

# Convert to data frame
gsea_df <- gsea_res@result

# Filter for immune-related pathways
# immune_keywords <- c("immune", "interferon", "cytokine", "T cell", "B cell", "inflammation", "lymphocyte", "leukocyte", "macrophage", "NK")
# immune_gsea <- gsea_df[grepl(paste(immune_keywords, collapse="|"), gsea_df$Description, ignore.case=TRUE), ]

# Filter for significance (FDR <= 0.25)
immune_gsea_sig <- immune_gsea[immune_gsea$p.adjust <= 0.25, ]
gsea_df_sig <- gsea_df[gsea_df$p.adjust <= 0.25, ]
write.csv(gsea_df_sig,paste0(Dir,"RespondersVsPD_pathways.csv"))

#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
# Fig.7
# Author: Chao Cheng; Lamis Naddaf
# 04/03/2026
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#+
#+#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
## how the HBV/HCV affect the TP53 score and TIS scores
rm(list=ls())
myinf1 = "./data/LIHC_RNAseqv2_Tumor_Symbol.rda"
myinf2 = "./data/LIHC_Clincial_info.txt"
myinf3 = "./data/Cao2016_SciRep_SupTab1_TCGA_Virus.txt"
myinf4 = "./data/TCGA_LIHC__GenomicEvent_iRAS.txt"
myinf5 = "./data/TCGA_LIHC_Freq_GenomicEvents.txt"
myinf6 = "./data/Thorsson_2018_TCGA_immunelandscape.csv"
figDir="./Figure7/"

#------------------
load(myinf1)
data = mydata
expr = log2(data+1)
tis_genes = c('CCL5', 'CD27', 'CD274', 'CD276', 'CD8A', 'CMKLR1', 'CXCL9', 'CXCR6',
              'HLA-DQA1', 'HLA-DRB1', 'IDO1', 'LAG3', 'NKG7', 'PDCD1LG2', 'PSMB10',
              'STAT1', 'TIGIT', 'TNFRSF9')  ## PMID: 28650338

se = which(row.names(expr) %in% tis_genes)
expr = expr[se,]
tis.score = apply(expr, 2, mean)

#--------------------------------------
info = read.table(myinf2, sep="\t", header=T, row.names=1, quote="")
se = c("vital_status", "days_to_death", "days_to_last_followup", "age_at_initial_pathologic_diagnosis", "gender", "stage_event.pathologic_stage")
info = info[,se]
xx = ifelse(!is.na(info$days_to_death), info$days_to_death, info$days_to_last_followup)
t.surv = as.numeric(xx)
e.surv = ifelse(info[, "vital_status"]=="dead", 1, 0)
info = cbind(t.surv, e.surv, info)
info = info[!is.na(info$t.surv), ]


data = read.table(myinf3, sep="\t", header=T, row.names=1, stringsAsFactors=F)
se = which(data$Cancer=="LIHC")
data = data[se,]
row.names(data) = substr(row.names(data), 1, 16)
data = data[, -1]
se = which(as.numeric(substr(row.names(data), 14, 15))==1)
data = data[se,]
row.names(data) = substr(row.names(data), 1, 12)
hcv.sam = row.names(data)[data$HCV>0]
hbv.sam = row.names(data)[data$HBV>0]
non.sam = row.names(data)[data$HCV==0 & data$HBV==0]
length(hcv.sam)
length(hbv.sam)
length(non.sam)

se = which(row.names(info) %in% row.names(data))
info = info[se,]
info$HCV = ifelse(row.names(info) %in% hcv.sam, 1, 0)
info$HBV = ifelse(row.names(info) %in% hbv.sam, 1, 0)
sum(row.names(info) %in% names(tis.score))
info$tis.score = tis.score[row.names(info)]

#--------------------------------------
data <- read.table(myinf4, header=T, sep="\t", row.names=1, quote="")
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

k = which(colnames(data) =="TP53__MUT")
p53.score = as.numeric(data[,k])
names(p53.score) = row.names(data)

info$p53.score = p53.score[row.names(info)]

#--------------------------------------
data <- read.table(myinf5, header=T, sep="\t", row.names=1, quote="")
P53.mut = data$TP53__MUT
names(P53.mut) = row.names(data)
P53.mut = P53.mut[!is.na(P53.mut)]
P53.mut = ifelse(P53.mut>=1, 1, 0)
names(P53.mut) = substr(names(P53.mut), 1, 12)
info$P53.mut = P53.mut[row.names(info)]

#--------------------------------------
mys1 = info$p53.score
mys2 = info$tis.score

mytag = rep("", nrow(info))
mytag[mys1>=median(mys1) & mys2>=median(mys2)] = "HH"
mytag[mys1>=median(mys1) & mys2<median(mys2)] = "HL"
mytag[mys1<median(mys1) & mys2>=median(mys2)] = "LH"
mytag[mys1<median(mys1) & mys2<median(mys2)] = "LL"
# mytag[mys1<median(mys1)] = "P53L"
table(mytag)
info$mytag = mytag




#------------------------------------------------------------
# Create viral group annotation
#------------------------------------------------------------
info$virus_group <- "Non-viral"
info$virus_group[info$HBV==1 & info$HCV==0] <- "HBV"
info$virus_group[info$HBV==0 & info$HCV==1] <- "HCV"
info$virus_group[info$HBV==1 & info$HCV==1] <- "HBV & HCV"
info$virus_group <- factor(info$virus_group,
                           levels=c("HBV & HCV","HCV","HBV","Non-viral"))

virus_cols <- c(
  "HBV & HCV"   = "brown3",
  "HCV"       = "mediumseagreen",
  "HBV"       = "dodgerblue2",
  "Non-viral" = "grey40"
)
library(ggplot2)
library(dplyr)

# Create summary dataframe
virus_df <- as.data.frame(table(info$virus_group))
colnames(virus_df) <- c("virus_group", "count")

# Calculate percentages and label text
virus_df <- virus_df %>%
  mutate(
    percent = count / sum(count) * 100,
    label = paste0(count, "\n(", sprintf("%.1f", percent), "%)")
  )

virus_df <- virus_df %>%
  mutate(
    percent = count / sum(count) * 100,
    label =count
  )

library(ggplot2)
library(dplyr)

virus_df <- virus_df %>%
  arrange(desc(virus_group)) %>%  # reverse if needed
  mutate(
    cum_count = cumsum(count),
    pos = cum_count - count / 2    # default midpoint
  ) %>%
  mutate(
    # Move "Non-viral" label to first quarter of its slice
    pos = ifelse(virus_group == "Non-viral",
                 cum_count - count * 0.75,   # 1st quarter: start + 1/4 of slice
                 pos)
  )


# Pie chart with counts inside and labels outside like x-axis
P <- ggplot(virus_df, aes(x = "x", y = count, fill = virus_group)) +
  geom_bar(stat = "identity", color = "white", width = 1) +
  coord_polar(theta = "y", direction = -1, clip = "off") +
  
  # Add counts inside slices
  geom_text(aes(y = pos, label = count),
            color = "white",
            size = 4,
            fontface = "bold") +
  
  # Set virus group labels outside the pie
  scale_y_continuous(
    breaks = virus_df$pos,
    labels = virus_df$virus_group
  ) +
  
  scale_fill_manual(values = virus_cols, guide = FALSE) +  # remove legend
  theme_void() +
  theme(
    axis.text.x = element_text(size = 12, face = "bold"),
    plot.title = element_text(hjust = 0.5, size = 15)
  ) 


graphics.off()
pdf(paste0(figDir,"Pie.pdf"),width = 2.5,height = 2.5)
print(P)
dev.off()

table(info$virus_group)
table(info$mytag[which(info$virus_group=="HBV & HCV")])
####### 


library(ggplot2)
library(dplyr)

# Contingency table
tab <- with(info, table(virus_group, mytag))

# Chi-square test
chi_res <- chisq.test(tab)

# Extract observed, expected, and standardized residuals
observed  <- chi_res$observed
expected  <- chi_res$expected
stdres    <- chi_res$stdres

# Create a data frame with OE ratio and significance
P2_df <- as.data.frame(as.table(observed)) %>%
  setNames(c("virus_group", "mytag", "observed")) %>%
  mutate(
    expected = as.vector(expected),
    OE_ratio = observed / expected,
    stdres   = as.vector(stdres),
    p_value  = 2 * pnorm(-abs(stdres)),   # two-sided
    p_adj    = p.adjust(p_value, method = "BH"),
    signif   = case_when(
      p_adj < 0.001 ~ "***",
      p_adj < 0.01  ~ "**",
      p_adj < 0.05  ~ "*",
      TRUE ~ ""
    ),
    virus_group = factor(virus_group,
                         levels = c("HBV & HCV","HCV","HBV","Non-viral"))
  )

# Plot O/E ratio with significance labels
P <- ggplot(P2_df, aes(x = mytag, y = OE_ratio, fill = virus_group)) +
  geom_col(width = 0.6) +
  geom_hline(yintercept = 1, linetype = "dashed") +
  
  # Add significance labels slightly above the bar
  geom_text(aes(label = signif, y = OE_ratio + 0.05),
            size = 5, fontface = "bold", vjust = 0) +
  
  facet_wrap(~ virus_group, ncol = 4) +
  scale_fill_manual(values = virus_cols) +
  labs(
    x = "TIS/TP53 category",
    y = "Observed / Expected ratio"
  ) +
  theme_classic(base_size = 12) +
  theme(
    legend.position = "none",
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12),
    axis.title.x = element_text(size = 14),
    axis.title.y = element_text(size = 14),
    strip.text = element_text(size = 12, face = "bold"),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8)
  )

P

graphics.off()
pdf(paste0(figDir,"5A.pdf"),width = 11,height = 3.5)
print(P)
dev.off()


#------------------------------------------------------------
# Figure 5.2A — TP53 score by viral status
#------------------------------------------------------------

# Define comparisons: Stage 1 vs 3 and Stage 1 vs 4
my_comparisons <- list(
  c("HBV", "Non-viral"),
  c("HCV", "HBV"),
  c("HBV & HCV", "HBV"))

my_comparisons <- list(
  c("HBV", "Non-viral"),
  c("HCV", "HBV"))

my_comparisons <- list(
  c("HBV", "Non-viral"),
  c("HCV", "Non-viral"))


P_viral_Score_box<-ggplot(info, aes(x=virus_group, y=p53.score, color=virus_group)) +
  geom_violin(fill = "grey90", color = NA, trim = FALSE, width = 1) +
  geom_boxplot(size=1,width = 0.6, outlier.shape = NA, fill = NA) +
  geom_jitter(width = 0.15, size = 0.5, alpha = 0.7) +
    scale_fill_manual(values=NA) +
  scale_color_manual(values=virus_cols) +
  theme(legend.position="none") +
  labs(
       x="", y="TP53 Score") +
  theme_classic() +
  theme(
    legend.position = "none",
    axis.text.x = element_text(size = 12,angle = 45, hjust = 1),     
    axis.text.y = element_text(size = 12),       
    axis.title.x = element_text(size = 14),      
    axis.title.y = element_text(size = 14)
    # ,      
    # plot.title = element_text(size = 16, face = "bold")
  ) +
    stat_compare_means(
    comparisons = my_comparisons,
    method = "wilcox.test",   
    position = c(4, 5),
    label = "p.signif",
    size = 5,label.y = c(4 ,4.7)  )+
  coord_cartesian(ylim = c(NA, 5))

graphics.off()
pdf(paste0(figDir,"5B.pdf"),width = 3.3,height = 3.5)
P_viral_Score_box
dev.off()

#------------------------------------------------------------
# Figure 5.2B — TIS score by viral status
#------------------------------------------------------------
my_comparisons <- list(
  c("HCV", "Non-viral"),
  c("HCV", "HBV"))

my_comparisons <- list(
  c("HBV", "Non-viral"),
  c("HCV", "Non-viral"))
P_viral_tis_box<-ggplot(info, aes(x=virus_group, y=tis.score, color=virus_group)) +
  geom_violin(fill = "grey90", color = NA, trim = FALSE, width = 1) +
  geom_boxplot(size=1,width = 0.6, outlier.shape = NA, fill = NA) +
  geom_jitter(width = 0.15, size = 0.5, alpha = 0.7) +
  scale_fill_manual(values=NA) +
  scale_color_manual(values=virus_cols) +
  theme(legend.position="none") +
  labs(
       x="", y="TIS Score") +
  theme_classic() +
  theme(
    legend.position = "none",
    axis.text.x = element_text(size = 12,angle = 45, hjust = 1),        
    axis.text.y = element_text(size = 12),       
    axis.title.x = element_text(size = 14),      
    axis.title.y = element_text(size = 14)
  ) +
  stat_compare_means(
    comparisons = my_comparisons,
    method = "wilcox.test",      
    label = "p.signif",
    size = 5  # Increase size of p-values
  )




  

graphics.off()
pdf(paste0(figDir,"5C.pdf"),width = 3.5,height = 3.5)
P_viral_tis_box
dev.off()

#------------------------------------------------------------
# Figure 5.2C — TP53 mutation frequencies by virus
#------------------------------------------------------------
library(dplyr)
library(ggplot2)

library(dplyr)

mutation_ratio <- info %>%
  filter(!is.na(P53.mut)) %>%
  group_by(virus_group) %>%
  summarise(
    total_n = n(),
    mutated_n = sum(P53.mut == 1),
    mutation_ratio = mean(P53.mut == 1),   # proportion
    mutation_percent = mutation_ratio * 100
  )

mutation_ratio

P_viral_mutaion_bar <- info %>%
  filter(!is.na(P53.mut)) %>%
  ggplot(aes(x = virus_group, fill = factor(P53.mut))) +
  geom_bar(width=0.7, position = "fill") +
  scale_fill_manual(
    name = "TP53",
    values = c("0" = "snow3",   # WT
               "1" = "indianred"),  # Mut
    labels = c("0" = "WT",
               "1" = "Mut")
  ) +
  labs(
    y = "Proportion",
    x = ""
  ) +
  theme_classic(base_size = 12) +
  theme(
    axis.text.x = element_text(size = 12, angle = 45, hjust = 1),
    axis.text.y = element_text(size = 12),
    axis.title.x = element_text(size = 14),
    axis.title.y = element_text(size = 14),
    legend.position = "top",
    legend.direction = "horizontal"
    # ,
    # legend.title = element_text(size = 10),  # smaller title
    # legend.text = element_text(size = 9),    # smaller labels
    # legend.key.size = unit(0.5, "lines")     # smaller key boxes
  )

# Fisher test between HBV and HCV for TP53 mutation

tab_hbv_hcv <- info %>%
  filter(virus_group %in% c("HBV", "HCV"),
         !is.na(P53.mut)) %>%
  with(table(virus_group, P53.mut))



fisher_res <- fisher.test(tab_hbv_hcv[2:3,])

fisher_res$p.value

tab_all <- info %>%
  filter(!is.na(P53.mut),
         !is.na(virus_group)) %>%
  with(table(virus_group, P53.mut))

fisher.test(tab_all)


graphics.off()
pdf(paste0(figDir,"5D.pdf"),width = 2.3,height = 3.5)
P_viral_mutaion_bar
dev.off()










########
info1 = read.table(myinf6, sep=",", header=T, row.names=1, quote="")
se = which(info1$TCGA.Study=="LIHC")
info1<-info1[se,]

commonID<-intersect(rownames(info),rownames(info1))
info1<-info1[commonID,]
info<-info[commonID,]


#------------------------------------------------------------
# Create viral group annotation
#------------------------------------------------------------

features <- c(    "Proliferation",
                  "IFN.gamma.Response",
                  "Lymphocyte.Infiltration.Signature.Score",
                  "Leukocyte.Fraction")

t <- c(    "Proliferation",
                  "IFN-gamma response",
                  "TIL level",
                  "Leukocyte fraction")


virus_cols <- c(
  "HBV & HCV"   = "brown3",
  "HCV"       = "mediumseagreen",
  "HBV"       = "dodgerblue2",
  "Non-viral" = "grey40"
)

my_comparisons <- list(
  c("HCV", "HBV"))

my_comparisons <- list(
  c( "HBV","Non-viral"),c("HCV", "Non-viral"))

table(info$virus_group)

####### 
library(ggplot2)


library(ggplot2)
library(ggpubr)

feat="Proliferation"
i=1
for (feat in features) {
  i=i+1
  
  # Add feature temporarily into info dataframe
  info$feature_value <- info1[rownames(info), feat]
  
  P <- ggplot(info, aes(x = virus_group,
                        y = feature_value,
                        color = virus_group)) +
    # geom_violin(fill = "grey90", color = NA, trim = FALSE, width = 1) +
    geom_boxplot(size=0.8, outlier.shape = NA, fill = NA) +
    # geom_jitter(width = 0.15, size = 0.5, alpha = 0.7) +
    scale_color_manual(values = virus_cols) +
    labs(
      title = t[i],
      x = "",
      y = t[i]
    ) +
    scale_y_continuous(
      limits = c(min(info$feature_value, na.rm = TRUE),
                 max(info$feature_value, na.rm = TRUE) * 2)
    )+
    theme_classic(base_size = 12) +
    theme(
      legend.position = "none",
       axis.text.x = element_text(size = 12,angle = 45, hjust = 1),
      axis.text.y = element_text(size = 12),
      axis.title.x = element_text(size = 14),
      axis.title.y = element_text(size = 14)
    ) +
    stat_compare_means(
      comparisons = my_comparisons,
      method = "wilcox.test",
      label = "p.signif",
      size = 5
      # ,
      # label.y = c(0.5 ,0.6)
    )
  
  # Save automatically
  pdf(paste0(figDir, "Feature_", feat, ".pdf"),
      width = 2.5, height = 3.5)
  print(P)
  dev.off()
}

graphics.off()


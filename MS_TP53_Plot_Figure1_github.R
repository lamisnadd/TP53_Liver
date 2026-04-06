#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
# Fig.1
# Author: Chao Cheng; Lamis Naddaf
# 04/03/2026
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

rm(list=ls())
FigDir<-"./Figure1/"
myinf1 = "./data/Fig_Data_4Datasets.rda"
## data.List, clin.List

load(myinf1)


library(survival)
library(survminer)
dataset<-"TCGA"
se = which(names(data.List)==dataset)
#se = which(names(data.List)=="NCI")ƒ
# se = which(names(data.List)=="UTSMC")

data = data.List[[se]]
info = clin.List[[se]]

commonID<-intersect(rownames(data),rownames(info))
data<-data[commonID,]
info<-info[commonID,]

# data<-data[-1,]
# info<-info[-1,]
surv = info$OS.time
event = info$OS.event

mycox = coxph(Surv(surv, event)~P53.mut, info) 
mycox = summary(mycox)
coxph.pval1_Mut = mycox$coefficients[5]
tmp = mycox$conf.int
hr1_Mut = tmp[1]
lb1_Mut = tmp[3]
ub1_Mut = tmp[4]


library(survival)

survreg.pval1 = survreg.pval2 = coxph.pval1 = coxph.pval2 =rep(0, ncol(data))
hr1 = lb1 = ub1 = hr2 =lb2 = ub2 = rep(0, ncol(data))
for(k in 1:ncol(data)){
  cat("\r", k)
  mytf = as.numeric(data[,k])
  xx = cbind(mytf, info)
  
  mycox = coxph(Surv(surv, event)~mytf, xx) 
  mycox = summary(mycox)
  coxph.pval1[k] = mycox$coefficients[5]
  tmp = mycox$conf.int
  hr1[k] = tmp[1]
  lb1[k] = tmp[3]
  ub1[k] = tmp[4]
}

coxph.pval1 = c(coxph.pval1_Mut,coxph.pval1)
hr1 = c(hr1_Mut,hr1)

coxph.qval1 = p.adjust(coxph.pval1, "BH")

name = rep("",15)
name[1]<-"TP53 Mutaion"
name[3]<-"TP53 Signature"
res = data.frame( coxph.pval1, hr1)


#---------------------------------------
data=res
data$category= ifelse(data$coxph.pval1>0.01, "NS", ifelse(data$hr1>1, "Hazardous", "Protective"))
data$label =  name

#---------------------------------------

my.color = c("red3", "grey", "skyblue")
my.label = c("Hazardous", "Not significant", "Protective")
my.xlab= "Log2HR"
my.ylab= "Log10(P-value)"
my.xlim= range(log2(data$hr1))
my.ylim= c(0, max(-log10(data$coxph.pval1)))
data$label =  rep("",15)
data$label[1]<-"TP53 Mutaion"
data$label[3]<-"TP53 Signature"
my.color = c( "grey", "black","red3")
my.label = c("TP53 Signature", "", "TP53 Mutaion")
data$name=name

library(tidyverse) # includes ggplot2, for data visualisation. dplyr, for data manipulation.
library(RColorBrewer) # for a colourful plot
library(ggrepel) # for nice annotations

mygg <- ggplot(data, aes(x = log2(hr1), y = -log10(coxph.pval1), col = label, label = label)) +
  geom_vline(xintercept = 0, col = "black", linetype = 'dashed') +
  geom_hline(yintercept = -log10(0.01), col = "black", linetype = 'dashed') + 
  geom_point(size = 2, show.legend = FALSE) + 
  scale_color_manual(values = my.color, 						## 
                     labels = my.label) + 		## 
  # geom_text_repel(
  #   data = subset(data, name == "TP53 Mutaion"),
  #   label = "TP53 mutation",
  #   color = "black",
  #   size = 4,
  #   vjust = -1,
  #   hjust = 0.5
  # )+

  guides(color="none") +	
  coord_cartesian(xlim = my.xlim, ylim = my.ylim) + 
  labs(color = '', #legend_title, 
       x = my.xlab, y = my.ylab) + 
  #  scale_x_continuous(breaks = seq(-10, 10, 2)) + 
  ggtitle(dataset) + 													   # Plot title 
  geom_text_repel(max.overlaps = Inf) + 											# To show all labels 
  theme_classic()+theme(
    axis.title.x = element_text(size = 14),		# face = "bold"
    axis.title.y = element_text(size = 14),
    axis.text.x  = element_text(size = 12),
    axis.text.y  = element_text(size = 12),
    plot.title   = element_text(size = 16, hjust = 0.5),  # center title
    legend.key.height = unit(0.05, "cm")
  )

graphics.off()
pdf(paste0(FigDir,"HR_OS_volcano.pdf"),width =2.3,height = 2.3)
mygg
dev.off()
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

myinf1 = "./data/Fig_Data_4Datasets.rda"

load(myinf1)

se = which(names(data.List)=="TCGA")
data = data.List[[se]]
score = data$TP53__MUT
names(score)<-rownames(data)
                       
hist(score)         
library(ggplot2)

df <- data.frame(score = score)

med <- median(df$score, na.rm = TRUE)

p <- ggplot(df, aes(x = score)) +
  geom_histogram(
    aes(y = after_stat(density)),
    bins = 30,
    fill = "grey90",
    color = "grey45") +
  geom_density(
    # color = "slateblue",
    color = "dodgerblue4",
    linewidth = 1
  ) +
  geom_vline(
    xintercept = med,
    # color = "tomato",
    color = "red3",
    
    linetype = "dashed",
    linewidth = 1
  ) +
  labs(
    x = "TP53 score",
    y = "Density",
    title = ""
  ) +
  theme_classic()+
  theme(
    axis.text.x = element_text(size = 11),       
    axis.text.y = element_text(size =11),       
    axis.title.x = element_text(size = 13),      
    axis.title.y = element_text(size = 13),      
    plot.title = element_text(size = 13,  hjust = 0.5))


p

graphics.off()
pdf(paste0(FigDir,"ScoreDistribution.pdf"), width = 2.8,height = 2)
print(p)
dev.off()
   
###############################

myinf5 = "./data/TCGA_LIHC_Freq_GenomicEvents.txt" 
myinf2 = "./data/LIHC_CNV_Symbol.rda"


#-------------------

## data.List, clin.List
myinf1 = "./data/Fig_Data_4Datasets.rda"
load(myinf1)

library(survival)
library(survminer)

# se = which(names(data.List)=="RIKEN")
se = which(names(data.List)=="TCGA")

data = data.List[[se]]
info = clin.List[[se]]

score = data$TP53__MUT
names(score)<-rownames(data)


raw.data <- cbind(score = score[rownames(info)], info)
raw.data <- raw.data[!is.na(raw.data$score), ]

data = raw.data

mycat = ifelse(data$score>median(data$score), "High", "Low")
data = cbind(mycat, data)
data$mycat <- factor(data$mycat, levels = c("High", "Low"))
data$time = data$OS.time
data$event = data$OS.event

#--------------------------------------
data1 <- read.table(myinf5, header=T, sep="\t", row.names=1, quote="")
P53.mut = data1$TP53__MUT
names(P53.mut) = row.names(data1)
P53.mut = P53.mut[!is.na(P53.mut)]
P53.mut = ifelse(P53.mut>=1, 1, 0)
names(P53.mut) = substr(names(P53.mut), 1, 12)

common_ids <- intersect(row.names(data), row.names(info))
data <- data[common_ids, ]
info <- info[common_ids, ]
P53.mut <- P53.mut[common_ids]
data$P53.mut <- P53.mut

library(ggplot2)
# Treat NA as a separate category
data$P53.mut.cat <- ifelse(is.na(data$P53.mut), "Unknown",
                           ifelse(data$P53.mut == 1, "MUT", "WT"))

# Convert to factor with desired order
data$P53.mut.cat <- factor(data$P53.mut.cat, 
                           levels = c("WT", "MUT", "Unknown"))


# -----------------------------
#   SETUP
# -----------------------------
library(ggplot2)
library(dplyr)
library(patchwork)
library(ggpubr)

# Factor order
data$P53.mut.cat <- factor(data$P53.mut.cat,
                           levels = c("WT",  "MUT","Unknown"))
cutoff <- median(data$score, na.rm = TRUE)

data$score_group <- ifelse(data$score >= cutoff, "high", "low")
data$score_group <- factor(data$score_group, levels = c("low", "high"))



cols <- c(
  "WT" = "royalblue3",
  "Unknown" = "grey40",
  "MUT" = "lightsalmon"
)
# -----------------------------
#   FREQUENCY PLOT (TOP)
# -----------------------------
freq_df <- data |> dplyr::count(P53.mut.cat)

P_mut_freq <- ggplot(freq_df, aes(x = P53.mut.cat, y = n, fill = P53.mut.cat)) +
  geom_col(alpha = 1, width = 0.3) +
  
  scale_fill_manual(values = cols) +
  
  labs(x = NULL,
       y = "Count",
       title = "TCGA SNPs") +
  
  # theme_bw(base_size = 14) +
  theme_classic() +
  theme(
    # axis.text.x = element_text(size = 12),       
    axis.text.y = element_text(size = 12),       
    # axis.title.x = element_text(size = 0),      
    axis.title.y = element_text(size = 14),      
    plot.title = element_text(size = 14,  hjust = 0.5)  ,
    axis.title.x = element_blank(),
    axis.text.x  = element_blank(),    # remove x labels on top
    axis.ticks.x = element_blank(),
    plot.margin = ggplot2::margin(0, 5, -7, 5, unit = "pt"),
    legend.position = "none"
  )


freq_df$ratio <- freq_df$n / sum(freq_df$n)
freq_df$label <- paste0(freq_df$n, " \n(", scales::percent(freq_df$ratio, accuracy = 0.1), ")")
freq_df$label <- freq_df$n
p<-ggplot(freq_df, aes(x = "", y = n, fill = P53.mut.cat)) +
  geom_col(width = 1, color = "white") +
  
  geom_text( aes(label = label),position = position_stack(vjust = 0.5), size = 5,color="white",fontface = "bold") +
  
  coord_polar(theta = "y") +
  
  scale_fill_manual(values = cols) +
  
  labs( title = "TCGA SNPs", fill = NULL) +
  
  theme_classic() +
  theme(
    axis.line = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.title = element_blank(),
    plot.title = element_text(size = 0, hjust = 0.5),
    legend.text = element_text(size = 10),   # <-- increase legend text
    legend.title = element_text(size = 10),  # optional
    
    # legend.position = "none",
    plot.margin = ggplot2::margin(0, 5, -7, 5, unit = "pt")
  )

graphics.off()
pdf(paste0(FigDir,"MutationDistribution.pdf"), width = 2.8,height = 2.2)
print(p)
dev.off()

freq_df2 <- data |>
  dplyr::count(P53.mut.cat, score_group, name = "n")



  
P_mut_freq2 <- ggplot(freq_df2,aes(x = P53.mut.cat, y = n, fill = score_group)) +
  # geom_col(width = 0.4) +
  geom_col(width = 0.4, position = position_dodge(width = 0.45)) +
  # scale_fill_manual(values = c( "royalblue3","lightsalmon")) +
  scale_fill_manual(values = c( "steelblue","coral1")) +
  scale_y_continuous(breaks = seq(0, max(freq_df2$n), by = 40))+
  labs( y = "Count",title = "TCGA SNPs",fill = "TP53 score") +
  theme_classic() +
  theme(
    axis.text.y  = element_text(size = 12),
    axis.title.y = element_text(size = 14),
    plot.title   = element_text(size = 14, hjust = 0.5),
    axis.title.x = element_blank(),
    axis.text.x  = element_blank(),
    axis.ticks.x = element_blank(),
    # legend.title  = element_text(size = 9),
    # legend.text   = element_text(size =9),
    legend.key.size = unit(0.4, "cm"),   # smaller boxes
    plot.margin  = ggplot2::margin(0, 5, -7, 5, unit = "pt")
  )




freq_df2$ratio <- freq_df2$n / sum(freq_df2$n)
freq_df2$label <- paste0(freq_df2$n, " \n(", scales::percent(freq_df2$ratio, accuracy = 0.1), ")")


freq_df3 <- data |>
  dplyr::count( score_group, name = "n")
freq_df3$ratio <- freq_df3$n / sum(freq_df3$n)
freq_df3$label <- paste0(freq_df3$n, " \n(", scales::percent(freq_df3$ratio, accuracy = 0.1), ")")



my_comparisons <- list(
  c( "Unknown","WT"),
  c("MUT","WT")
)

# Calculate fold changes
fold_df <- lapply(my_comparisons, function(x) {
  grp1 <- data$score[data$P53.mut.cat == x[1]]
  grp2 <- data$score[data$P53.mut.cat == x[2]]
  data.frame(
    group1 = x[1],
    group2 = x[2],
    fold_change = round(mean(grp1, na.rm = TRUE) / mean(grp2, na.rm = TRUE), 1)
  )
}) |> bind_rows()

# Positions for labels
fold_df$y_pos <- max(data$score, na.rm = TRUE) * c(1.15, 1.3)  # stagger labels

# Plot
P_mut_score_viobox<-ggplot(data, aes(x = P53.mut.cat, y = score, color = P53.mut.cat)) +
  
  geom_violin(fill = "grey85", color = NA, trim = FALSE, width = 1.2) +
  geom_boxplot(width = 0.6, outlier.shape = NA, fill = NA,linewidth=0.8) +
  geom_jitter(width = 0.15, size = 0.3, alpha = 0.7) +
  geom_hline(yintercept = median(data$score, na.rm = TRUE), linetype = "dashed", linewidth = 0.8,color = "black"  ) +
  scale_color_manual(values = cols) +
   ylim(NA, 8)+
  stat_compare_means(
    comparisons = my_comparisons,
    method = "wilcox.test",
    label = "p.signif",
    # label = "p.format",
    size = 6, 
    label.y = c(5.2, 6.5),
    method.args = list(alternative = "greater"),  ### this is not working, my value is significant if it worked
    tip.length = 0.03
  ) +
  
  # Add fold change labels
  geom_text(
    data = fold_df,
    aes(
      x = group1,
      y = y_pos,
      label = paste0("FC=", fold_change)
    ),
    inherit.aes = FALSE,
    color = "black",
    size = 4
  ) +
  
  labs(x = "TP53 mutation status",
       y = "TP53 Score",
       title="") +
  
  # theme_bw(base_size = 14) +
  theme_classic() +
  theme(
    axis.text.x = element_text(size = 12, angle = 45, hjust = 1),      
    axis.text.y = element_text(size = 12),       
    axis.title.x = element_text(size = 14),      
    axis.title.y = element_text(size = 14),      
    plot.title = element_text(size = 14,  hjust = 0.5)  ,
    plot.margin = ggplot2::margin(0, 5, -7, 5, unit = "pt"),
    legend.position = "none"
  )



# -----------------------------
#   COMBINE — bar plot on TOP
# -----------------------------
plot1_a <- P_mut_freq2 /
  P_mut_score_viobox +
  plot_layout(heights = c(1, 3))  # adjust relative sizes

plot1_a

graphics.off()
pdf(paste0(FigDir,"1A.pdf"), width = 3.5,height = 3.7)
    print(plot1_a)
dev.off()

########### deletion and score 
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# [5.4] TP53-del vs. wt
# rm(list=ls())
myinf1 = "./data/TCGA_LIHC__GenomicEvent_iRAS.txt"
myinf2 = "./data/LIHC_CNV_Symbol.rda"

#-------------------
data <- read.table(myinf1, header=T, sep="\t", row.names=1, quote="")
cnum = ncol(data)/4
means<-colMeans(data)
a=2*cnum  +1
plot(colMeans(data)[1:cnum],colMeans(data)[a: (3*cnum)])

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

score = data[, "TP53__MUT"]
names(score) = row.names(data)

#-------------------
load(myinf2)
cnv = mydata
se = which(substr(colnames(cnv), 14,15)=="01") 
cnv = cnv[,se]
colnames(cnv) = substr(colnames(cnv), 1, 12)

thr = log2(1.6/2)
sam.mu = colnames(cnv)[cnv["TP53", ]<thr]
sam.mu = sam.mu[!is.na(sam.mu)]
sam.wt = colnames(cnv)[cnv["TP53", ]>=thr]
sam.wt = sam.wt[!is.na(sam.wt)]

se = which(names(score) %in% sam.mu)
xx1 = score[se]
se = which(names(score) %in% sam.wt)
xx2 = score[se]
length(xx1)			## 20
length(xx2)			## 7
boxplot(list(Del=xx1, WT=xx2))
wilcox.test(xx1, xx2, alternative="g")

df <- data.frame(
  Score = c(xx1, xx2),
  Group = factor(
    c(rep("DEL", length(xx1)),
      rep("WT",  length(xx2))),
    levels = c("WT", "DEL")
  )
)


cols <- c(
  "WT" = "royalblue3",
  "NA" = "grey30",
  "DEL" = "firebrick3"
)


my_comparisons <- list(
  c("DEL","WT")
)


fc_mean <- mean(xx1, na.rm = TRUE) / mean(xx2, na.rm = TRUE)
log2_fc_mean <- log2(fc_mean)


# Positions for labels
y_pos <- max(df$Score, na.rm = TRUE) * c(1.05, 1.15)  # stagger labels
ymax <- max(df$Score, na.rm = TRUE)

# Plot
P_del_score_box<-ggplot(df, aes(x = Group, y = Score, color = Group)) +
   # geom_violin(fill = "lightgrey", color = NA, trim = FALSE, width = 1.4) +
  geom_boxplot(width = 0.6, outlier.shape = NA, fill = NA,linewidth=1) +
   # geom_jitter(width = 0.15, size = 0.7, alpha = 0.7) +
  scale_color_manual(values = cols) +
  stat_compare_means(
    comparisons = my_comparisons,
    method = "wilcox.test",
    label = "p.signif",
    #label = "p.format",
    label.y = 4 ,
    size=6,
    method.args =  list(alternative = "greater"),  ### this is not working, my value is significant if it worked
    tip.length = 0.03
  ) +
  labs(x = "TP53 Copy Number",
       y = "TP53 Score",
       title = "TCGA CNV") +
  
  theme_classic(base_size = 12) +
 theme(
  legend.position = "none",
  axis.text.x = element_text(size = 12),       
  axis.text.y = element_text(size = 12),       
  axis.title.x = element_text(size = 13),      
  axis.title.y = element_text(size = 14),      
  plot.title = element_text(size = 14,  hjust = 0.5)  
)+coord_cartesian(ylim = c(NA, 5))

fig1b_1<-P_del_score_box
## W = 20860, p-value = 9.285e-06
graphics.off()
pdf(paste0(FigDir,"1B.pdf"), width = 2,height = 3)
print(fig1b_1)
dev.off()

#### cell line data
## compare TP53 signature score between TP53-mut and wt liver cancer cell lines (CCLE data)


## TP53 mutation
# rm(list=ls())
myinf1 = "./data/CCLE_LIHC_CellLine_TP53_score.txt"
myinf2 = "./data/CCLE_summarized_non_synonymous_mutation_count.xls"

#------------------
score = read.table(myinf1, sep="\t", header=T, quote="", row.names=1)
tmp = score[,2]
names(tmp) = row.names(score)
score = tmp

#------------------
mut = read.table(myinf2, sep="\t", header=T, quote="", row.names=1)
se = which(mut["TP53", ]>0)
sam1 = colnames(mut)[se]
se = which(mut["TP53", ]==0)
sam2 = colnames(mut)[se]

se = which(names(score) %in% sam1)
xx1 = score[se]
se = which(names(score) %in% sam2)
xx2 = score[se]
length(xx1)			## 20
length(xx2)			## 7
boxplot(list(MUT=xx1, WT=xx2))
wilcox.test(xx1, xx2, alternative="g")
## W = 116, p-value = 0.004658

library(ggplot2)
library(ggpubr)  # for stat_compare_means

# Combine scores and mutation status into a single data frame
df <- data.frame(
  Score = c(xx1, xx2),
  Group = factor(c(rep("MUT", length(xx1)), rep("WT", length(xx2))),
                 levels = c("WT", "MUT"))
)


cols <- c(
  "WT" = "royalblue3",
  "Unknown" = "grey30",
  "MUT" = "lightsalmon"
)
my_comparisons <- list(
  c("MUT","WT")
)

fig1b<-ggplot(df, aes(x = Group, y = Score, color = Group)) +
  geom_boxplot(width = 0.6, outlier.shape = NA, fill = NA,linewidth=1) +
  geom_jitter(width = 0.15, size = 0.7, alpha = 0.7) +
  scale_color_manual(values = cols) +
  stat_compare_means(
    comparisons = my_comparisons,
    method = "wilcox.test",
    alternative = "greater",
    size=6,
    method.args =  list(alternative = "greater"),
    label = "p.signif"
  ) +
  labs(
     x = "TP53 mutation status",
    y = "TP53 score",
    title = "LIHC cell lines"
  ) +
  theme_classic(base_size = 12) +
  theme(
    legend.position = "none",
    axis.text.x = element_text(size = 12),       
    axis.text.y = element_text(size = 12),       
    axis.title.x = element_text(size = 13),      
    axis.title.y = element_text(size = 14),      
    plot.title = element_text(size = 14,  hjust = 0.5)) +
  coord_cartesian(ylim = c(NA, 29))


graphics.off()
pdf(paste0(FigDir,"1B_2.pdf"), width = 2.3,height = 3.5)
print(fig1b)
dev.off()


######## expression vs protein
# rm(list=ls())
myinf1 = "./data/TCGA_LIHC__GenomicEvent_iRAS.txt"
myinf2  = "./data/LIHC_RPPA.txt"
myinf3 = "./data/LIHC_RNAseqv2_Tumor_Symbol.rda"
load(file= myinf3)
data = mydata	##  samples

expr = data["TP53", ]
names(expr) = colnames(data)
expr = log2(expr+1)
#-------------------
data <- read.table(myinf1, header=T, sep="\t", row.names=1, quote="")
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

score = data[, "TP53__MUT"]
names(score) = row.names(data)

#-------------------
data <- read.table(myinf2, header=T, sep="\t", row.names=1, quote="", check.names=F)
se = which(substr(colnames(data), 14,15)=="01")
data = data[,se]
colnames(data) = substr(colnames(data), 1, 12)


comxx = intersect(intersect(colnames(data), names(score)),names(expr))
data = data[, comxx]
score = score[comxx]
expr=expr[comxx]

xx = cor(t(data), score, method="s")
xx[order(xx),]

xx["p53",] 
-0.612268 

p53_protein <- as.numeric(data["p53", ])
tp53_score  <- as.numeric(score[colnames(data)])
tp53_expr   <- as.numeric(expr[colnames(data)])  # log2(TP53 + 1)

  df_p53 <- data.frame(
         Sample = colnames(data),
         p53_protein = p53_protein,
        TP53_score  = tp53_score,
        TP53_expr   = tp53_expr
     )
  
 
  
  ## remove NAs
  df_p53 <-   df_p53[complete.cases(  df_p53), ]
  library(ggplot2)
  
  ct <- cor.test(
    df_p53$TP53_expr,
    df_p53$TP53_score,
    method = "spearman"
  )
  
  rho  <- ct$estimate
  pval <- ct$p.value
  
  
  # Plot with correlation 
  fig1S2<-ggplot(df_p53, aes(x = TP53_score, y = TP53_expr)) +
    geom_point(size = 1.2, alpha = 0.8, color = "steelblue") +
    # geom_point(size = 1.2, alpha = 0.8, color = "seashell4") +
    geom_smooth(method = "lm", se = FALSE, color = "black", linetype = "dashed") +
    annotate(
      "text",
      x =  Inf, y = Inf,
      label = paste0(
        # "\u03C1 = ", round(rho, 2),
        "Rho= ", round(rho, 2),
         ",P = ", signif(pval, 1)
      ),
      hjust = 1.3, vjust = 1.9,
      size = 4.5
    ) +
    theme_classic(base_size = 12) +
    labs(
      x = "TP53 score",
      y = "TP53 Expression"
    )+
    theme(
      legend.position = "none",
      axis.text.x = element_text(size = 12),       
      axis.text.y = element_text(size = 12),       
      axis.title.x = element_text(size = 14),      
      axis.title.y = element_text(size = 14)
      # ,      
      # plot.title = element_text(size = 14,  hjust = 0.5)
      )    
  # plot.title = element_text(size = 14,  hjust = 0.5)) +
  # coord_cartesian(ylim = c(NA, 4.8))+coord_cartesian(ylim = c(NA, 29))
  
  graphics.off()
  pdf(paste0(FigDir,"1S11.pdf"), width = 3,height = 3.5)
  print(fig1S2)
  dev.off()
  
  
  

  
  # Compute Spearman correlation and p-value
  ct <- cor.test(df_p53$p53_protein, df_p53$TP53_score, method = "spearman")
  rho  <- ct$estimate
  pval <- ct$p.value
  
 # Plot with correlation 
  fig1c<-ggplot(df_p53, aes(x = p53_protein, y = TP53_score)) +
    geom_point(size = 1.2, alpha = 0.8, color = "steelblue") +
    # geom_point(size = 1.2, alpha = 0.8, color = "seashell4") +
    geom_smooth(method = "lm", se = FALSE, color = "black", linetype = "dashed") +
    annotate(
      "text",
      x =  Inf, y = Inf,
      label = paste0(
        # "\u03C1 = ", round(rho, 2),
         "Rho= ", round(rho, 2),
        # ",P = ", signif(pval, 3)
        ",P ", "< 2e-16"
      ),
      hjust = 1.3, vjust = 1.9,
      size = 4.5
    ) +
    theme_classic(base_size = 12) +
    labs(
      x = "P53 protein level (RPPA)",
      y = "TP53 score",
      title = " "
    )+
    theme(
      legend.position = "none",
      axis.text.x = element_text(size = 12),       
      axis.text.y = element_text(size = 12),       
      axis.title.x = element_text(size = 14),      
      axis.title.y = element_text(size = 14),      
      plot.title = element_text(size = 14,  hjust = 0.5))    
      # plot.title = element_text(size = 14,  hjust = 0.5)) +
    # coord_cartesian(ylim = c(NA, 4.8))+coord_cartesian(ylim = c(NA, 29))
  
  graphics.off()
  pdf(paste0(FigDir,"1C.pdf"), width = 3,height = 3.5)
  print(fig1c)
  dev.off()
## negatively correlated with P53 protein level

## use TP53 vs. p53 as control (lower than)
  
  #########################
  myinf5<-"./data/TCGA_LIHC_Freq_GenomicEvents.txt"
  myinf2  = "./data/LIHC_RPPA.txt"
  
  #-------------------
  #--------------------------------------
  data1 <- read.table(myinf5, header=T, sep="\t", row.names=1, quote="")
  P53.mut = data1$TP53__MUT
  names(P53.mut) = row.names(data1)
  P53.mut = P53.mut[!is.na(P53.mut)]
  P53.mut = ifelse(P53.mut>=1, 1, 0)
  names(P53.mut) = substr(names(P53.mut), 1, 12)
  
  
  data <- read.table(myinf2, header=T, sep="\t", row.names=1, quote="", check.names=F)
  se = which(substr(colnames(data), 14,15)=="01")
  data = data[,se]
  colnames(data) = substr(colnames(data), 1, 12)
  
  comxx = intersect(colnames(data), (names(P53.mut)))
  data = data[, comxx]
  P53.mut <- P53.mut[comxx]
  
  xx = cor(t(data), score, method="s")
  xx[order(xx),]
  data["p53",] 

  df <- data.frame(
    mutation = factor(
      ifelse(P53.mut == 1, "MUT", "WT"),
      levels = c("WT", "MUT")
    ),
    RPPA = as.numeric(data["p53", , drop = TRUE])
  )
 
  
  cols <- c(
    "WT" = "royalblue3",
    "NA" = "grey30",
    "MUT" = "darksalmon"
  )
  

  my_comparisons<-list(c("WT","MUT"))
  P_box_RPPA_mutaion<-ggplot(df, aes(x = mutation, y = RPPA, color = mutation)) +
    # geom_violin(fill = "lightgrey", color = NA, trim = FALSE, width = 1.4) +
    geom_boxplot(width = 0.6, outlier.shape = NA, fill = NA,linewidth=1) +
    # geom_jitter(width = 0.15, size = 0.7, alpha = 0.7) +
    scale_color_manual(values = cols) +
    stat_compare_means(
      comparisons = my_comparisons,
      method = "wilcox.test",
      label = "p.signif",
      #label = "p.format",
      label.y =0.5,
      size=6,
      method.args =  list(alternative = "greater"),  ### this is not working, my value is significant if it worked
      tip.length = 0.03
    ) +
    labs(x = "TP53 mutation status",
         y = "RPPA",
         title = "") +
    
    theme_classic(base_size = 12) +
    theme(
      legend.position = "none",
      axis.text.x = element_text(size = 12),       
      axis.text.y = element_text(size = 12),       
      axis.title.x = element_text(size = 13),      
      axis.title.y = element_text(size = 14),      
      plot.title = element_text(size = 14,  hjust = 0.5)  
    ) +coord_cartesian(ylim = c(NA, 0.7))
  
  
  
  
  
  graphics.off()
  pdf(paste0(FigDir,"P_box_RPPA_mutaion.pdf"), width = 2.4,height = 3)
  print(P_box_RPPA_mutaion)
  dev.off()
#--------------
# rm(list=ls())
myinf1 = "./data/LIHC_RNAseqv2_Tumor_Symbol.rda"
myinf2  = "./data/LIHC_RPPA.txt"

load(file= myinf1)
data = mydata	##  samples

score = data["TP53", ]
names(score) = colnames(data)
score = log2(score+1)

#-------------------
data <- read.table(myinf2, header=T, sep="\t", row.names=1, quote="", check.names=F)
se = which(substr(colnames(data), 14,15)=="01")
data = data[,se]
colnames(data) = substr(colnames(data), 1, 12)


comxx = intersect(colnames(data), names(score))
data = data[, comxx]
score = score[comxx]

xx = cor(t(data), as.numeric(score), method="s")
xx["p53", ]		## -0.2297533 


## extract paired vectors
p53_protein <- as.numeric(data["p53", ])
tp53_expr   <- as.numeric(score[colnames(data)])  # log2(TP53 + 1)

df_expr_prot <- data.frame(
  Sample = colnames(data),
  p53_protein = p53_protein,
  TP53_expr   = tp53_expr
)

## remove NAs
df_expr_prot <- df_expr_prot[complete.cases(df_expr_prot), ]
library(ggplot2)

ct <- cor.test(
  df_expr_prot$p53_protein,
  df_expr_prot$TP53_expr,
  method = "spearman"
)

rho  <- ct$estimate
pval <- ct$p.value


fig1d<-ggplot(df_expr_prot, aes(x = p53_protein, y = TP53_expr)) +
  geom_point(size = 1.2, alpha = 0.8, color = "steelblue") +
  geom_smooth(method = "lm", se = FALSE,
              color = "black", linetype = "dashed") +
  annotate(
    "text",
    x = Inf, y = Inf,
    label = paste0(
      # "\u03C1 = ", round(rho, 2),
      "Rho = ", round(rho, 2),
      ",P = ", signif(pval, 1)
    ),
    hjust = 1.0, vjust = 1.9,
    size = 5
  ) +
  theme_classic(base_size = 12) +
  labs(
    x = "P53 protein level (RPPA)",
    y = "TP53 expression"
  )+
  theme(
    legend.position = "none",
    axis.text.x = element_text(size = 12),       
    axis.text.y = element_text(size = 12),       
    axis.title.x = element_text(size = 14),      
    axis.title.y = element_text(size = 14))    

graphics.off()
pdf(paste0(FigDir,"1D.pdf"), width = 3,height = 3.5)
print(fig1d)
dev.off()

########
#########################
myinf5<-"./data/TCGA_LIHC_Freq_GenomicEvents.txt"
myinf2 = "./data/LIHC_RNAseqv2_Tumor_Symbol.rda"
#-------------------
#--------------------------------------
data1 <- read.table(myinf5, header=T, sep="\t", row.names=1, quote="")
P53.mut = data1$TP53__MUT
names(P53.mut) = row.names(data1)
P53.mut = P53.mut[!is.na(P53.mut)]
P53.mut = ifelse(P53.mut>=1, 1, 0)
names(P53.mut) = substr(names(P53.mut), 1, 12)


load(file= myinf2)
data = mydata	##  samples

score = data["TP53", ]
names(score) = colnames(data)
score = log2(score+1)

#-------------------

comxx = intersect(names(P53.mut), names(score))
P53.mut = P53.mut[comxx]
score = score[comxx]

df <- data.frame(
  mutation = factor(
    ifelse(P53.mut == 1, "MUT", "WT"),levels = c("WT", "MUT")), Expression = unlist(score))


cols <- c(
  "WT" = "royalblue3",
  "NA" = "grey30",
  "MUT" = "darksalmon"
)


my_comparisons<-list(c("WT","MUT"))
P_box_Expression_mutaion<-ggplot(df, aes(x = mutation, y = Expression, color = mutation)) +
  # geom_violin(fill = "lightgrey", color = NA, trim = FALSE, width = 1.4) +
  geom_boxplot(width = 0.6, outlier.shape = NA, fill = NA,linewidth=1) +
  # geom_jitter(width = 0.15, size = 0.7, alpha = 0.7) +
  scale_color_manual(values = cols) +
  stat_compare_means(
    comparisons = my_comparisons,
    method = "wilcox.test",
    label = "p.signif",
    #label = "p.format",
    label.y =12,
    size=6,
    method.args =  list(alternative = "greater"),  ### this is not working, my value is significant if it worked
    tip.length = 0.03
  ) +
  labs(x = "TP53 mutation status",
       y = "TP53 expression",
       title = "") +
  
  theme_classic(base_size = 12) +
  theme(
    legend.position = "none",
    axis.text.x = element_text(size = 12),       
    axis.text.y = element_text(size = 12),       
    axis.title.x = element_text(size = 13),      
    axis.title.y = element_text(size = 14),      
    plot.title = element_text(size = 14,  hjust = 0.5)  
  )  +coord_cartesian(ylim = c(NA, 13))



graphics.off()
pdf(paste0(FigDir,"P_box_Expression_mutaion.pdf"), width = 2.4,height = 3)
print(P_box_Expression_mutaion)
dev.off()



#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#----------------------------
myinf1 = "./data/LIHC_RNAseqv2_Tumor_Symbol.rda"
myinf2  = "/mount/ictr1/chenglab/cc59/PubDat/Dataset/Firehose/Methylation/Gene_Promoter_Methyl/LIHC_Methy450KAvg_Promoter_beta.rda"

load(myinf1)
data = mydata
myx = as.numeric(data["TP53", ])
names(myx) = colnames(data)

#-------------------
load(myinf2)
data = mydata
se = which(substr(colnames(data), 14,15)=="01")
data = data[,se]
colnames(data) = substr(colnames(data), 1, 12)
myy = data["TP53", ]

comxx = intersect(names(myx), names(myy))
myx = myx[comxx]
myy = myy[comxx]

plot(myx, myy)
hist(myy, br=30, xlim=c(0.45, 0.8))

df<-data.frame(myth=myy)
p_methyl_distribution <- ggplot(df, aes(x = myth)) +
  geom_histogram(fill = "grey90", color = "black",alpha = 1, bins = 30) +
  geom_histogram(
    data = subset(df, myth > 0.7),
    aes(y = ..count..),
    fill = "red3", alpha = 0.8, bins = 30)+
  scale_x_continuous(limits = c(0.45, 0.8)) +
  labs(
    x = "Methylation level",
    y = "Count",
    title="") +
  theme_classic()+
  theme(
    legend.position = "none",
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12),
    axis.title.x = element_text(size = 14),
    axis.title.y = element_text(size = 14),      
    plot.title = element_text(size = 14,  hjust = 0.5)
  )

graphics.off()
pdf(paste0(FigDir,"p_methyl_distribution.pdf"), width = 4,height = 3)
print(p_methyl_distribution)
dev.off()

###########
myinf1 = "./data/TCGA_LIHC__GenomicEvent_iRAS.txt"
myinf2 = "./data/LIHC_RNAseqv2_Tumor_Symbol.rda"
#-------------------
data <- read.table(myinf1, header=T, sep="\t", row.names=1, quote="")
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

score = data[, "TP53__MUT"]
names(score) = row.names(data)

#-------------------
load(myinf2)
data = mydata
myx = as.numeric(data["TP53", ])
names(myx) = colnames(data)

comxx = intersect(names(myx), names(score))
myx = myx[comxx]
score = score[comxx]

plot(score, myx)

df<-data.frame(p53_score=score, expression=myx)
index<-which(rownames(df)=="TCGA-CC-A9FU")
P_scatter_score_expression <- ggplot(df, aes(x = p53_score, y = expression)) +
  geom_point(size = 1, alpha = 0.8, color = "grey") +
  theme_classic(base_size = 12) +
  labs(x = "TP53 score", y = "Expression") +
  geom_vline(
    xintercept = median(df$p53_score, na.rm = TRUE),
    linetype = "dashed",
    color = "dodgerblue4",size=1.02 ) +
  geom_hline(
    yintercept = median(df$expression, na.rm = TRUE),
    linetype = "dashed",
    color = "dodgerblue4",size=1.02) +
  geom_point(data = df[index,],
             color = "red3", size = 2) +
  annotate("text", x = 3, y = 0.7,
           label = "TCGA-CC-A9FU", size = 3, vjust = -1) +
  theme(
    legend.position = "none",
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12),
    axis.title.x = element_text(size = 14),
    axis.title.y = element_text(size = 14))


graphics.off()
pdf(paste0(FigDir,"P_scatter_score_expression.pdf"), width = 2,height = 2)
print(P_scatter_score_expression)
dev.off()

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
## corr between MDM2 and P53-score
myinf1 = "./data/TCGA_LIHC__GenomicEvent_iRAS.txt"
myinf2  = "./data/LIHC_CNV_Symbol.rda"
myinf3 = "./data/TCGA_LIHC_Freq_GenomicEvents.txt"

#-------------------
data <- read.table(myinf3, header=T, sep="\t", row.names=1, quote="")
P53.mut = data$TP53__MUT
names(P53.mut) = row.names(data)
se = which(P53.mut==1)
mut.sam = names(P53.mut)[se]
se = which(P53.mut==0)
wt.sam = names(P53.mut)[se]
mut.sam = substr(mut.sam, 1, 12)
wt.sam = substr(wt.sam, 1, 12)


#-------------------
data <- read.table(myinf1, header=T, sep="\t", row.names=1, quote="")
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

score = data[, "TP53__MUT"]
names(score) = row.names(data)

#-------------------
load(myinf2)
data = mydata
se = which(substr(colnames(data), 14,15)=="01")
data = data[,se]
colnames(data) = substr(colnames(data), 1, 12)

comxx = intersect(colnames(data), names(score))
data = data[, comxx]
score = score[comxx]


#-------------------
myx = data["MDM2", ]
cor(score, myx, method="s")
se = which(names(myx) %in% mut.sam)
cor(score[se], myx[se], method="s")
se = which(names(myx) %in% wt.sam)
cor(score[se], myx[se], method="s")

plot(score[se], myx[se])


# Extract MDM2 CNV
mdm2 <- data["MDM2", ]
# Harmonize samples
common <- intersect(names(score), names(mdm2))
score <- score[common]
mdm2  <- mdm2[common]

# Mutation status
mutation <- ifelse(names(score) %in% mut.sam, "MUT", "WT")

# MDM2 CNV status
mdm2_status <- ifelse(mdm2 > median(mdm2), "Amplified", "Normal")
df <- data.frame(
  mutation = factor(mutation, levels = c("WT", "MUT")),
  MDM2_CNV = factor(mdm2_status, levels = c("Normal", "Amplified")),
  score = as.numeric(score))

library(dplyr)


# compute p-values per mutation group
stat_df <- df %>%
  group_by(mutation) %>%
  summarise(p = wilcox.test(score ~ MDM2_CNV)$p.value) %>%
  mutate(label = cut(
    p,
    breaks = c(-Inf, 0.001, 0.01, 0.05, Inf),
    labels = c("***", "**", "*", "ns")
  ))

ymax <- max(df$score, na.rm = TRUE)

box_MDM2<-ggplot(df, aes(x = mutation, y = score, color = MDM2_CNV)) +
  geom_boxplot(position = position_dodge(0.7), width = 0.5, outlier.shape = NA, linewidth = 1) +
  scale_color_manual(values = c("grey70","salmon3")) +
  annotate("segment",
           x = c(0.8, 1.8), xend = c(1.2, 2.2),
           y = ymax * 1.05, yend = ymax * 1.05) +
  annotate("text",
           x = c(1, 2),
           y = ymax * 1.1,
           label = stat_df$label) +
  labs(
    x = "TP53 mutation status",
    y = "TP53 Score",
    color = "MDM2 CNV"
  ) +
  theme_classic(base_size = 12)+
  theme(
    axis.text.x = element_text(size = 12),       
    axis.text.y = element_text(size = 12),       
    axis.title.x = element_text(size = 14),      
    axis.title.y = element_text(size = 14),      
    plot.title = element_text(size = 14,  hjust = 0.5)  
  )

graphics.off()
pdf(paste0(FigDir,"box_MDM2t.pdf"), width = 4,height = 3)
print(box_MDM2)
dev.off()

########

myinf1 = "./data/TCGA_LIHC__GenomicEvent_iRAS.txt"
myinf2  = "./data/LIHC_CNV_Symbol.rda"
myinf3 = "./data/TCGA_LIHC_Freq_GenomicEvents.txt"

#-------------------
data <- read.table(myinf3, header=T, sep="\t", row.names=1, quote="")
P53.mut = data$TP53__MUT
names(P53.mut) = row.names(data)
se = which(P53.mut==1)
mut.sam = names(P53.mut)[se]
se = which(P53.mut==0)
wt.sam = names(P53.mut)[se]
mut.sam = substr(mut.sam, 1, 12)
wt.sam = substr(wt.sam, 1, 12)


#-------------------
data <- read.table(myinf1, header=T, sep="\t", row.names=1, quote="")
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

score = data[, "TP53__MUT"]
names(score) = row.names(data)

#-------------------
load(myinf2)
data = mydata
se = which(substr(colnames(data), 14,15)=="01")
data = data[,se]
colnames(data) = substr(colnames(data), 1, 12)

comxx = intersect(colnames(data), names(score))
data = data[, comxx]
score = score[comxx]


#-------------------
myx = data["MDM2", ]
cor(score, myx, method="s")
se = which(names(myx) %in% mut.sam)
cor(score[se], myx[se], method="s")
se = which(names(myx) %in% wt.sam)
cor(score[se], myx[se], method="s")

plot(score[se], myx[se])



# Extract MDM2 CNV values
myx <- data["MDM2", ]

# Combine score and CNV into a single data frame
df <- data.frame(
  score = score,
  MDM2_CNV = myx,
  mutation = ifelse(names(score) %in% mut.sam, "MUT", "WT")
)


library(dplyr)

# Compute Spearman rho and p-value for each mutation group
cors <- df %>%
  group_by(mutation) %>%
  summarise(
    test = list(cor.test(score, MDM2_CNV, method = "spearman"))
  ) %>%
  mutate(
    rho = sapply(test, function(x) x$estimate),
    pval = sapply(test, function(x) x$p.value)
  )

# Prepare annotation text
cors_text <- paste0( "Rho = ", round(cors$rho, 2),
                     ", P = ", signif(cors$pval, 2))
# Plot with annotations
fig_MDM2<-ggplot(df, aes(x = score, y = MDM2_CNV, color = mutation)) +
  geom_point(size = 1, alpha = 0.7) +
  geom_smooth(method = "lm", se = FALSE, fullrange = TRUE) +
  scale_color_manual(values = c("lightsalmon","royalblue3")) +
  annotate("text", x = max(df$score)*0.6, y = max(df$MDM2_CNV)*0.9, 
           label = cors_text[1], color = "lightsalmon") +
  annotate("text", x = max(df$score)*0.6, y = max(df$MDM2_CNV)*0.75, 
           label = cors_text[2], color = "royalblue3") +
  labs(
    x = "TP53 Score",
    y = "MDM2 CNV",
    color = "TP53 Status",
    title = ""
  ) +
  theme_classic(base_size = 14) +
  theme(
    # legend.position = "none",
    axis.text.x = element_text(size = 12),       
    axis.text.y = element_text(size = 12),       
    axis.title.x = element_text(size = 14),      
    axis.title.y = element_text(size = 14)
  )


graphics.off()
pdf(paste0(FigDir,"fig_MDM2.pdf"), width = 4.5,height = 3.5)
print(fig_MDM2)
dev.off()


##########





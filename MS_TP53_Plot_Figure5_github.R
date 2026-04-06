#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
# Fig.5
# Author: Chao Cheng; Lamis Naddaf
# 04/03/2026
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

rm(list=ls())
FigDir<-"./Figure5/"
mycol <-c(
  "CR" = "steelblue4",
  "PR" = "steelblue2", 
  "SD" = "grey60", 
  "PD" = "indianred4")
######################
######################

myinf1 = "./data/Zhu_GO30140_IMbrave150_atezo_bev.rda"
myinf2 = "./data/Zhu_GO30140_IMbrave150_atezo_bev__GenomicEvent_iRAS.txt"


#---------------------------
data <- read.table(myinf2, header=T, sep="\t", row.names=1, quote="")
cnum = ncol(data)/2
data = data[, 1:cnum]
tmp = colnames(data)
tmp = gsub("\\.ES", "", tmp)
colnames(data) = tmp
cnum = ncol(data)/2
dat1 = data[,1:cnum]
dat2 = data[, (cnum+1):(2*cnum)]
xx = dat1-dat2
colnames(xx) = gsub("\\.up", "", colnames(dat1))
data = xx
for(k in 1:ncol(data))
{
  data[,k] = data[,k]/sd(abs(data[,k]))
}
data = data[grep("mul.adj__", colnames(data))]
colnames(data) = gsub("mul.adj__", "", colnames(data))


#---------------------------
load(file= myinf1)
info = info	##  samples

expr = log2(rna+1)
tis_genes = c('CCL5', 'CD27', 'CD274', 'CD276', 'CD8A', 'CMKLR1', 'CXCL9', 'CXCR6',
              'HLA-DQA1', 'HLA-DRB1', 'IDO1', 'LAG3', 'NKG7', 'PDCD1LG2', 'PSMB10',
              'STAT1', 'TIGIT', 'TNFRSF9')  ## PMID: 28650338

se = which(row.names(expr) %in% tis_genes)
expr = expr[se,]
tis.score = apply(expr, 2, mean)

info = cbind(tis.score, info)
se = which(info$Visit =="Pre-treatment")
info = info[se,]
se = which(!info$Confirmed.Response_IRF %in% c("", "NE"))
info = info[se, ]


#---------------------------
comxx = intersect(row.names(data), row.names(info))
data = data[comxx,]
info = info[comxx,]
info = info[, 1:9]


#--------------------------
#--------------------------
mys1 = data$TP53__MUT
mys2 = info$tis.score
plot(mys1, mys2)
abline(v=median(mys1))
abline(h=median(mys2))





library(ggplot2)

mytag = rep("", length(info$tis.score))
mytag[mys1>=median(mys1) & mys2>=median(mys2)] = "HH"
mytag[mys1>=median(mys1) & mys2<median(mys2)] = "HL"
mytag[mys1<median(mys1) & mys2>=median(mys2)] = "LH"
mytag[mys1<median(mys1) & mys2<median(mys2)] = "LL"

# mytag[mys1<median(mys1)] = "P53L"
table(mytag)
info$mytag = mytag

df_scatter <- data.frame(
  p53.score = data$TP53__MUT,
  tis.score = info$tis.score,
  cat=info$mytag
)


# Make sure the grouping variable is a factor
df_scatter$cat <- factor(df_scatter$cat, levels = c("LH","LL", "HL", "HH"))

# Compute medians
median_score <- median(df_scatter$p53.score, na.rm = TRUE)
median_tis   <- median(df_scatter$tis.score, na.rm = TRUE)

# Scatter plot
fig3i<-ggplot(df_scatter, aes(x = p53.score, y = tis.score, color = cat)) +
  geom_point(size = 1.5, alpha = 0.6) +
  
  # Add horizontal and vertical median lines
  geom_vline(xintercept = median_score, linetype = "dashed", color = "grey40") +
  geom_hline(yintercept = median_tis,   linetype = "dashed", color = "grey40") +
  
  # Labels and colors
  labs(
    x = "TP53 score",
    y = "TIS score",
    color = "Group",
    title = "Ate+Bev"
  ) +
  # scale_color_manual(values = c("LL" = "lightgoldenrod",
  #                               "LH" = "grey",
  #                               "HL"   = "darkseagreen",
  #                               "HH"   = "slateblue")) +
  # 
  # scale_color_manual(values = c(
  #   "LL" = "royalblue4",
  #   "LH" = "seagreen" ,
  #   "HL" = "yellowgreen",
  #   "HH" = "gold"
  # ))+


  scale_color_manual(values = c(
    "LL" = "#0D3B66",    # deeper, slightly muted blue for low/low
    "LH" = "#2E8B90",    # strong forest green for low/high
    "HL" = "palegreen3",    # bright yellow-green for high/low
    "HH" = "lightgoldenrod2"    # rich golden yellow for high/high
  ))+

  
  theme_classic()+
  theme(
    axis.text.x = element_text(size = 12),       
    axis.text.y = element_text(size = 12),       
    axis.title.x = element_text(size = 14),      
    axis.title.y = element_text(size = 14),      
    plot.title = element_text(size = 16, face = "bold")
    # ,
    # legend.position = "none"
  )


graphics.off()
pdf(paste0(FigDir,"3I.pdf"),width=3,height = 3.5)
print(fig3i)
dev.off()

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


myinf1 = "./data/ImmunoTargetedTherapy2/Zhu_GO30140_IMbrave150_atezo_bev.rda"
myinf2 = "./data/Zhu_GO30140_IMbrave150_atezo_bev__GenomicEvent_iRAS.txt"


#---------------------------
data <- read.table(myinf2, header=T, sep="\t", row.names=1, quote="")
cnum = ncol(data)/2
data = data[, 1:cnum]
tmp = colnames(data)
tmp = gsub("\\.ES", "", tmp)
colnames(data) = tmp
cnum = ncol(data)/2
dat1 = data[,1:cnum]
dat2 = data[, (cnum+1):(2*cnum)]
xx = dat1-dat2
colnames(xx) = gsub("\\.up", "", colnames(dat1))
data = xx
for(k in 1:ncol(data))
{
  data[,k] = data[,k]/sd(abs(data[,k]))
}
data = data[grep("mul.adj__", colnames(data))]
colnames(data) = gsub("mul.adj__", "", colnames(data))


#---------------------------
load(file= myinf1)
info = info	##  samples

expr = log2(rna+1)
tis_genes = c('CCL5', 'CD27', 'CD274', 'CD276', 'CD8A', 'CMKLR1', 'CXCL9', 'CXCR6',
              'HLA-DQA1', 'HLA-DRB1', 'IDO1', 'LAG3', 'NKG7', 'PDCD1LG2', 'PSMB10',
              'STAT1', 'TIGIT', 'TNFRSF9')  ## PMID: 28650338

se = which(row.names(expr) %in% tis_genes)
expr = expr[se,]
tis.score = apply(expr, 2, mean)

info = cbind(tis.score, info)
se = which(info$Visit =="Pre-treatment")
info = info[se,]
se = which(!info$Confirmed.Response_IRF %in% c("", "NE"))
info = info[se, ]


#---------------------------
comxx = intersect(row.names(data), row.names(info))
data = data[comxx,]
info = info[comxx,]
info = info[, 1:9]


#--------------------------
#--------------------------
mys1 = data$TP53__MUT
mys2 = info$tis.score
plot(mys1, mys2)
abline(v=median(mys1))
abline(h=median(mys2))

mytag = rep("", nrow(data))
mytag[mys1>=median(mys1) & mys2>=median(mys2)] = "HH"
mytag[mys1>=median(mys1) & mys2<median(mys2)] = "HL"
mytag[mys1<median(mys1) & mys2>=median(mys2)] = "LH"
mytag[mys1<median(mys1) & mys2<median(mys2)] = "LL"

# mytag[mys1<median(mys1)] = "P53L"
table(mytag)


xx = info
xx$mytag = mytag
tab=table(xx$mytag, xx$Confirmed.Response_IRF)
# fisher.test(tab)
raw.data = xx

data = raw.data
mycat = data$Confirmed.Response_IRF
mys = data$mytag
data = data.frame(mys, mycat)
se = which(data$mycat !="")
data = data[se,]



table(data$mycat)

xx= table(data$mys, data$mycat)

# raw.data = xx

data$mys <- factor(data$mys, levels = c("HH", "HL", "LH", "LL"))
data$mycat <- factor(data$mycat, levels = c("CR", "PR", "SD", "PD"))
# data<-data[-which(is.na(data$mycat)),]

mycol <-c(
  "CR" = "steelblue4",
  "PR" = "steelblue2", 
  "SD" = "grey60", 
  "PD" = "indianred4")


plot_data <- data %>%
  group_by(mys, mycat) %>%
  summarise(n = n(), .groups = 'drop') %>%
  group_by(mys) %>%
  mutate(freq = n / sum(n),
         label = paste0(round(freq * 100), "%"))


# Table of Response × Immune
tab <- table(data$mys, data$mycat)
tab
# a=fisher.test(tab, workspace = 2e8)


# label_df <- data.frame(
#   mys = 2,
#   y = 1.05,         
#   label = paste0("P= ",a$p.value)
# )

mygg4<-ggplot(plot_data,aes(y = freq, x = mys, fill = mycat)) +
  geom_flow(aes(alluvium = mycat), alpha= .5, color = "white",curve_type = "linear",   width = .7) +
  geom_col(width = .7, color = "white") +
  labs(title='Ate+Bev',x = "Patient group", y = "Proportion", fill = "Response") +
  scale_fill_manual(values = mycol) + 
  scale_y_continuous(NULL, expand = c(0,0)) +
  cowplot::theme_minimal_hgrid() +
  geom_text(
    aes(label = label),
    position = position_stack(vjust = 0.5),  # center inside each segment
    size = 3.5,
     # fontface = "bold",
    color = "black"
  ) +

  # geom_text(data = label_df, aes(x = mys, y = y, label = label),
  #           inherit.aes = FALSE,
  #           vjust = 0, hjust = 0.5,
  #           size = 4)  +
  theme_classic()+
  theme(
    axis.text.x = element_text(size = 12),       
    axis.text.y = element_text(size = 12),       
    axis.title.x = element_text(size = 14),      
    axis.title.y = element_text(size = 14),      
    plot.title = element_text(size = 16, face = "bold")
    # ,
    # legend.position = "none"
  )+ylim(c(NA,1.05))

# posx = 0.5
# posy = 1.1
# mytxt <- grobTree(textGrob("P= 1.2e-05", x=posx,  y=posy, hjust=0, vjust=1, gp=gpar(col="black", fontsize=14)))
# mygg4 <- mygg4 + annotation_custom(mytxt)
# mygg4


fig3j<-mygg4
graphics.off()
pdf(file =paste0(FigDir,"3J_22.pdf"), width = 4.5, height = 3.5) 
mygg4
dev.off()


data$mys <- factor(data$mys, levels = c("HH", "HL", "LH", "LL"))
data$mycat <- factor(data$mycat, levels = c("CR", "PR", "SD", "PD"))

data$mys 
data$mycat


tab <- table(data$mycat, data$mys )
tab


row_tot <- rowSums(tab)
col_tot <- colSums(tab)
N <- sum(tab)
exp <- outer(row_tot, col_tot) / N
oe <- (tab+1) / (exp+1)

pvals <- matrix(NA, nrow=nrow(tab), ncol=ncol(tab))

for(i in 1:nrow(tab)){
  for(j in 1:ncol(tab)){
    
    a <- tab[i,j]
    b <- sum(tab[i,-j])
    c <- sum(tab[-i,j])
    d <- sum(tab[-i,-j])
    
    mat <- matrix(c(a,b,c,d), nrow=2)
    
    pvals[i,j] <- fisher.test(mat)$p.value
  }
}

rownames(pvals) <- rownames(tab)
colnames(pvals) <- colnames(tab)
plot_df <- data.frame(
  Response = rep(rownames(tab), times=ncol(tab)),
  Group = rep(colnames(tab), each=nrow(tab)),
  Observed = as.vector(tab),
  Expected = as.vector(exp),
  OE = as.vector(oe),
  Pvalue = as.vector(pvals)
)

plot_df$logP <- -log10(plot_df$Pvalue)

plot_df$signedP <- plot_df$logP * sign(log2(plot_df$OE))

plot_df$Group <- factor(plot_df$Group,
                        levels=c("HH", "HL", "LH", "LL"))
plot_df$Response <- factor(plot_df$Response, levels = c("CR","PR","SD","PD"))

plot_df$color_val <- ifelse(plot_df$Pvalue > 0.2, NA, plot_df$signedP)
plot_df$color_val <- plot_df$signedP


plot_df$signedP <- plot_df$logP * sign(log2(plot_df$OE))


plot_df$logOE <- log2(plot_df$OE)
P<-ggplot(plot_df, aes(x = Response, y = Group)) +
  geom_point(aes(color = logOE, size = logP)) +
  labs(title='Ate+Bev',x = "Response", y = "TP53 X TIS") +
  scale_size_continuous(
    name = "-log10(P)",
    range = c(3, 8),      # min and max point size
    breaks = scales::pretty_breaks(n = 3)
  ) +
  scale_color_gradientn(
    # colors = c("steelblue4","lightcyan","white","rosybrown2","#B40426"),  
    colors = c("steelblue4","#8DB0FE","white","indianred2","#B40426"),
    values = scales::rescale(c(min(plot_df$logOE), 0, max(plot_df$logOE))),
    limits = c(min(plot_df$logOE), max(plot_df$logOE)),
    name = "log2(OE)"
  ) +
  
  theme_bw()+
  theme(
    axis.text.x = element_text(size = 12,angle = 90),       
    axis.text.y = element_text(size = 12),       
    axis.title.x = element_text(size = 14),      
    axis.title.y = element_text(size = 14),      
    plot.title = element_text(size = 16, face = "bold")  
  )


graphics.off()
pdf(file = paste0(FigDir,"dot_tis_score_4cat.pdf"), width = 3.3, height = 3.5) 
P
dev.off()

#######
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

data = raw.data
mycat = data$Confirmed.Response_IRF
mys = data$mytag
data = data.frame(mys, mycat)
se = which(data$mycat !="")
data = data[se,]



table(data$mycat)

xx= table(data$mys, data$mycat)
# myp = fisher.test(xx, alternative="g")$p.value
#myp = ifelse(myp<0.001, formatC(myp, format = "e", digits = 0), signif(myp, 1))
# txt.myp = paste("P=", myp, sep="")


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
library(ggplot2)
library(grid)
library(dplyr)

mytag = rep("", nrow(data))
mytag[mys1>=median(mys1) & mys2>=median(mys2)] = "HH"
mytag[mys1>=median(mys1) & mys2<median(mys2)] = "HL"
mytag[mys1<median(mys1) & mys2>=median(mys2)] = "LH"
mytag[mys1<median(mys1) & mys2<median(mys2)] = "LL"
# mytag[mys1<median(mys1)] = "P53L"
table(mytag)


xx = info
xx$mytag = mytag
table(xx$mytag, xx$Confirmed.Response_IRF)

raw.data = xx

data$mys <- factor(data$mys, levels = c("HH", "HL", "LH", "LL"))
data$mycat <- factor(data$mycat, levels = c("CR", "PR", "SD", "PD"))

plot_data <- data %>%
  group_by(mys, mycat) %>%
  summarise(n = n(), .groups = 'drop') %>%
  group_by(mys) %>%
  mutate(freq = n / sum(n),
         label = paste0(round(freq * 100), "%"))


# Table of Response × Immune
tab <- table(data$mys, data$mycat)
tab
# fisher.test(tab, workspace = 2e8)

# fisher.test(as.matrix(plot_data[,1:3]))

label_df <- data.frame(
  mys = 2,
  y = 1.05,         
  label = "P= 1.2e-05"
)

mygg4<-ggplot(plot_data,aes(y = freq, x = mys, fill = mycat)) +
  geom_flow(aes(alluvium = mycat), alpha= .5, color = "white",curve_type = "linear",   width = .7) +
  geom_col(width = .7, color = "white") +
  labs(title='Ate+Bev',x = "Patient group", y = "Proportion", fill = "Response") +
  scale_fill_manual(values = mycol) + 
  scale_y_continuous(NULL, expand = c(0,0)) +
  cowplot::theme_minimal_hgrid() +
  geom_text(
    aes(label = label),
    position = position_stack(vjust = 0.5),  # center inside each segment
    size = 3.5,
    # fontface = "bold",
    color = "black"
  ) +
# +  geom_text(data = label_df, aes(x = mys, y = y, label = label),
#                  inherit.aes = FALSE,
#                  vjust = 0, hjust = 0.5,
#                  size = 4) +
  # overall p-value label
  geom_text(data = label_df, aes(x = mys, y = y, label = label),
            inherit.aes = FALSE,
            vjust = 0, hjust = 0.5,
            size = 4)  +
theme_classic()+
  theme(
    axis.text.x = element_text(size = 12),       
    axis.text.y = element_text(size = 12),       
    axis.title.x = element_text(size = 14),      
    axis.title.y = element_text(size = 14),      
    plot.title = element_text(size = 16, face = "bold")
    # ,
    # legend.position = "none"
  )+ylim(c(NA,1.05))



posx = 0.5
posy = 1.1
mytxt <- grobTree(textGrob("P= 1.2e-05", x=posx,  y=posy, hjust=0, vjust=1, gp=gpar(col="black", fontsize=14)))
 mygg4 <- mygg4 + annotation_custom(mytxt)
mygg4

fig3j<-mygg4
graphics.off()
pdf(file =paste0(FigDir,"3J.pdf"), width = 4.5, height = 3.5) 
mygg4
dev.off()


#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

library(survival)
library(survminer)

mytag = rep("", nrow(data))
mytag[mys1>=median(mys1) & mys2>=median(mys2)] = "HH"
mytag[mys1>=median(mys1) & mys2<median(mys2)] = "HL"
mytag[mys1<median(mys1)] = "P53L"
table(mytag)


xx = info
xx$mytag = mytag
table(xx$mytag, xx$Confirmed.Response_IRF)

data = xx



mycat = data$mytag
data = cbind(mycat, data)
# data$mycat <- factor(data$mycat, levels = c("HH", "HL", "LH","LL"))
data$mycat <- factor(data$mycat, levels = c("HH", "HL", "P53L"))
data$time = data$PFS..in.days.IRF/30
# data$time = data$OS.in.days/30
data$event = 1-data$PFS.censoring..1.cens.0.evt..IRF
# data$event = 1-data$OS.censoring..1.cens.0.evt.
# data$time = data$OS.in.days
#   data$event = 1-data$OS.censoring..1.cens.0.evt.

data_hot=data[which(data$mycat%in%c("HH", "LH")),]
data_cold=data[which(data$mycat%in%c("HL", "LL")),]

fit <- survfit(Surv(time, event) ~ mycat, data = data)

fit_hot <- survfit(Surv(time, event) ~ mycat, data = data_hot)
fit_cold <- survfit(Surv(time, event) ~ mycat, data = data_cold)


mygg5 <- ggsurvplot(
  fit,
  data = data,
  risk.table = F,               # adds number at risk
  pval = TRUE,                     
  pval.coord = c(min(data$time), 0.05),
  # adds p-value
  conf.int = F,                 # confidence interval
  legend.title = "Group",
  legend.labs = c("HH", "HL", "P53L"),
   palette = c("lightgoldenrod1", "palegreen3",  "grey60"),
  xlab = "PFS(months)",
  ylab = "Survival Probability",
  break.time.by = 200,              # interval on x-axis
  risk.table.height = 0.25,       # adjust table height
  risk.table.y.text.col = TRUE,
  risk.table.y.text = FALSE,
  ggtheme = theme_classic() +
    theme(
      axis.title.x = element_text(size = 14),
      axis.title.y = element_text(size = 14),
      axis.text.x  = element_text(size = 12),
      axis.text.y  = element_text(size = 12),
      plot.title   = element_text(size = 16, hjust = 0.5)  # center title
    ),
  
)

fig1k<-mygg5
graphics.off()
pdf(paste0(FigDir,"3K.pdf"), width = 3, height = 3.5) 
print(mygg5$plot)
dev.off()

#--------------------------




#
#==============================
##heatmap of immune genes
#==============================


mm.inh = c("CTLA4", "PDCD1", "LAG3", "BTLA", "CD160", "IDO1", "IL10", "IL10RB", "TGFB1", "TGFBR1", "VTCN1", "CD244", "LGALS9", "HAVCR2", "ADORA2A", "TIGIT", "CSF1R", "KIR2DL1", "KIR2DL2", "KIR2DL3", "KDR", "CD96", "PVRL2", "C10orf54")
imm.sti = c("MICA", "MICB", "CD27", "CD274", "CD28", "CD40", "CD40LG", "CD70", "CD80", "CD86", "ICOS", "ICOSLG", "IL6", "IL6R", "PDCD1LG2", "TMEM173", "TNFRSF13B", "TNFRSF13C", "TNFRSF14", "TNFRSF17", "TNFRSF18", "TNFRSF4", "TNFRSF9", "TNFSF13", "TNFSF13B", "TNFSF18", "TNFSF4", "TNFSF9", "TNFSF15", "TNFRSF25", "HHLA2", "TMIGD2", "BTNL2", "CD276", "CD48", "TNFSF14", "TNFRSF8", "PVR", "LTA",  "IL2RA", "ENTPD1", "NT5E", "CXCR4", "CXCL12", "KLRK1", "NKG2A", "RAET1E", "ULBP1")
imm.oth = c("GZMA", "PRF1")
tis_genes = c('CCL5', 'CD27', 'CD274', 'CD276', 'CD8A', 'CMKLR1', 'CXCL9', 'CXCR6',
              'HLA-DQA1', 'HLA-DRB1', 'IDO1', 'LAG3', 'NKG7', 'PDCD1LG2', 'PSMB10',
              'STAT1', 'TIGIT', 'TNFRSF9')  ## PMID: 28650338
selected_Genes<-c("PDCD1", "CD274", "PDCD1LG2", "ICOS", "IDO1", "IFNG")
# Define response groups
info$ResponseGroup <- ifelse(
  info$Confirmed.Response_IRF %in% c("CR", "PR"),
  "CR/PR",
  ifelse(info$Confirmed.Response_IRF== "PD", "PD", NA)
)

# Keep only samples with defined groups
info_sub <- info[!is.na(info$ResponseGroup), ]



boxPlot <- function(genelist,Title,yLim,yPos){
  
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(ggpubr)
  
  # Keep only genes of interest
  genes_use <- intersect(genelist, rownames(rna))
  rna_sub <- rna[genes_use, rownames(info_sub)]
  
  # Transpose so samples are rows
  rna_t <- t(rna_sub)
  
  # Add group info
  rna_t <- cbind(rna_t, Group = info_sub$ResponseGroup)
  rna_t <- as.data.frame(rna_t)
  rna_t[genes_use] <- lapply(rna_t[genes_use], as.numeric)
  
  # Mean expression per gene per group
  gene_group_mean <- sapply(
    genes_use,
    function(g) {
      tapply(rna_t[[g]], rna_t$Group, mean, na.rm = TRUE)
    }
  )
  
  gene_group_mean <- as.data.frame(t(gene_group_mean))
  gene_group_mean$Gene <- rownames(gene_group_mean)
  rownames(gene_group_mean) <- NULL
  
  # Long format
  df_plot <- gene_group_mean %>%
    pivot_longer(
      cols = c("CR/PR", "PD"),
      names_to = "Group",
      values_to = "MeanExpression"
    )
  my_comparisons <- list(
    c("CR/PR", "PD"))
  
  # Plot with paired significance
  P <- ggplot(df_plot, aes(x = Group, y = MeanExpression, color = Group)) +
    geom_boxplot(size=1) +
     # geom_jitter(width = 0.15, size = 0.5, alpha = 0.7) +
    scale_color_manual(values = c("steelblue4","indianred3"))  +
    stat_compare_means(
      comparisons = my_comparisons,
      method = "wilcox.test",
      paired = T,
      label = "p.signif",
      size=4.5,
      label.y = yPos
    )+
    theme_classic() +
    labs(
      title = Title,
      # subtitle = "Each point = mean expression of one gene",
      x = "Response",
      y = "TIS"
    ) +
    theme(
      axis.text.x = element_text(size = 12),       
      axis.text.y = element_text(size = 12),       
      axis.title.x = element_text(size = 14),      
      axis.title.y = element_text(size = 14),      
      plot.title = element_text(size = 16,face="bold"),
      legend.position = "none"
    ) + coord_cartesian(ylim = yLim)
  
  return(P)
}



Inhibitors_Box<-boxPlot(mm.inh,"Immune Inhibitors",c(NA,55),50)
Stimulators_Box<-boxPlot(imm.sti,"Immune Stimulators",c(NA,30),22)
All_Box<-boxPlot(c(imm.sti,mm.inh,imm.oth),"All",c(NA,55),50)
Others_Box<-boxPlot(imm.oth,"Others",c(NA,55),50)
TIS_Box<-boxPlot(tis_genes,"IMbrave50",c(NA,22),17)
SelectedGenes_Box<-boxPlot(selected_Genes,"Selected Imm",c(NA,10),9)

graphics.off()
pdf(paste0(FigDir,"3F2.pdf"),width = 1.8, height = 3.5)
print(TIS_Box)
dev.off()



load(file= myinf1)

info$ResponseGroup <- ifelse(
  info$Confirmed.Response_IRF == "CR", "CR",
  ifelse(info$Confirmed.Response_IRF == "PR", "PR",
         ifelse(info$Confirmed.Response_IRF == "SD", "SD",
                ifelse(info$Confirmed.Response_IRF == "PD", "PD", NA)
         )
  )
)





##########################################
## Fig0_Boxplot_TP53_4Cat_Atex  -- Anti-PD-L1 only

myinf1 = "./data/Immunotherapy2/Zhu_GO30140_IMbrave150.rda"
myinf2 = "./data/Zhu_GO30140_IMbrave150_atezo__GenomicEvent_iRAS.txt"
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
se = which(info$Confirmed.Response_IRF %in% c("CR", "PR", "SD", "PD"))
info = info[se, ]

comxx = intersect(row.names(data), row.names(info))
data = data[comxx,]
info = info[comxx,]
raw.data = data
raw.info = info

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

data = raw.data
info = raw.info
mycat = info$Confirmed.Response_IRF
mys = data$TP53__MUT
data = data.frame(mys, mycat)
se = which(data$mycat !="")
data = data[se,]

myf = data$mycat
xx = table(myf)
myList = list(NULL)
for(k in 1:length(xx))
{
  se = which(myf==names(xx)[k])
  myList[[k]] = data$mys[se]
}
names(myList) = names(xx)

boxplot(myList)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
library(ggplot2)
library(grid)

data$mycat <- factor(data$mycat, levels = c("CR", "PR", "SD", "PD"))



my_comparisons <- list(
  # c("CR", "SD"),
  c("PR", "SD")
  ,c("PD", "PR")
)
mygg <- ggplot(data, aes(x=mycat, y=mys,  color=mycat)) + 
  # geom_boxplot(size=0.7,outlier.size = 0.3,linewidth=0.5) +
  geom_boxplot(size=1, outlier.shape = NA) +
  scale_fill_manual(values="white") +
  scale_color_manual(values=mycol) +
  # scale_fill_manual(values=mycol) +
  labs(title='Ate mono',x= 'Response', y = 'TP53 score')+
  theme_classic() + 
  theme(legend.position="none") +
  theme(
    axis.text.x = element_text(size = 12),       
    axis.text.y = element_text(size = 12),       
    axis.title.x = element_text(size = 14),      
    axis.title.y = element_text(size = 14),      
    plot.title = element_text(size = 16, face = "bold")  
  )+stat_compare_means(
    comparisons = my_comparisons,
    method = "wilcox.test",      
    label = "p.signif",
    size = 4.5   # Increase size of p-values
  )+coord_cartesian(ylim = c(NA, 3.5))

graphics.off()
pdf(paste0(FigDir,"5A.pdf"), width = 2.5, height = 3.5) 
mygg
dev.off()


#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

data = raw.data
info = raw.info
mycat = info$Confirmed.Response_IRF
mys = ifelse(data$TP53__MUT>median(data$TP53__MUT), "High", "Low")
data = data.frame(mys, mycat)
se = which(data$mycat !="")
data = data[se,]

xx= table(data$mys, data$mycat)
myp = fisher.test(xx, alternative="g")$p.value
myp = ifelse(myp<0.001, formatC(myp, format = "e", digits = 0), signif(myp, 1))
txt.myp = paste("P=", myp, sep="")


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
library(ggplot2)
library(grid)
library(dplyr)

data$mys <- factor(data$mys, levels = c( "High", "Low"))
data$mycat <- factor(data$mycat, levels = c("CR", "PR", "SD", "PD"))






# Table of Response × Immune
tab <- table(data$mys, data$mycat)
tab
fisher.test(tab, workspace = 2e8)

fisher.test(as.matrix(plot_data[,1:3]))

# Add a label above each bar group
label_df <- data.frame(
  mys = 1.5,
  y = 1.05,         
  label = "P = 0.008"
)


plot_data <- data %>%
  group_by(mys, mycat) %>%
  summarise(n = n(), .groups = 'drop') %>%
  group_by(mys) %>%
  mutate(freq = n / sum(n),
         label = paste0(round(freq * 100), "%"))


mygg2<-ggplot(plot_data,aes(y = freq, x = mys, fill = mycat)) +
  geom_flow(aes(alluvium = mycat), alpha= .5, color = "white",curve_type = "linear",   width = .7) +
  geom_col(width = .7, color = "white") +
  labs(title = "Ate mono",x = "TP53 score", y = "Proportion", fill = "Response") +
  scale_fill_manual(values = mycol) + 
  scale_y_continuous(NULL, expand = c(0,0)) +
  # theme_bw() +
  geom_text(
    aes(label = label),
    position = position_stack(vjust = 0.5),  # center inside each segment
    size = 4,
    # fontface = "bold",
    color = "black"
  ) +  geom_text(data = label_df, aes(x = mys, y = y, label = label),
                 inherit.aes = FALSE,
                 vjust = 0, hjust = 0.5,
                 size = 4)  +
  theme_classic()+
  theme(
    axis.text.x = element_text(size = 12),       
    axis.text.y = element_text(size = 12),       
    axis.title.x = element_text(size = 14),      
    axis.title.y = element_text(size = 14),      
    plot.title = element_text(size = 16, face = "bold")
    # ,
    # legend.position = "none"
  )+ylim(c(NA,1.05))






graphics.off()
pdf(file = paste0(FigDir,"5B.pdf"), width = 3.2, height = 3.5) 
mygg2
dev.off()

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
## -- Anti-PD-L1 combine TP53 with TLS
# rm(list=ls())
myinf1 = "./data/Immunotherapy2/Zhu_GO30140_IMbrave150.rda"
myinf2 = "./data/Zhu_GO30140_IMbrave150_atezo__GenomicEvent_iRAS.txt"

#---------------------------
data <- read.table(myinf2, header=T, sep="\t", row.names=1, quote="")
cnum = ncol(data)/2
data = data[, 1:cnum]
tmp = colnames(data)
tmp = gsub("\\.ES", "", tmp)
colnames(data) = tmp
cnum = ncol(data)/2
dat1 = data[,1:cnum]
dat2 = data[, (cnum+1):(2*cnum)]
xx = dat1-dat2
colnames(xx) = gsub("\\.up", "", colnames(dat1))
data = xx
for(k in 1:ncol(data))
{
  data[,k] = data[,k]/sd(abs(data[,k]))
}
data = data[grep("mul.adj__", colnames(data))]
colnames(data) = gsub("mul.adj__", "", colnames(data))


#---------------------------
load(file= myinf1)
info = info	##  samples

expr = log2(rna+1)
tis_genes = c('CCL5', 'CD27', 'CD274', 'CD276', 'CD8A', 'CMKLR1', 'CXCL9', 'CXCR6',
              'HLA-DQA1', 'HLA-DRB1', 'IDO1', 'LAG3', 'NKG7', 'PDCD1LG2', 'PSMB10',
              'STAT1', 'TIGIT', 'TNFRSF9')  ## PMID: 28650338

se = which(row.names(expr) %in% tis_genes)
expr = expr[se,]
tis.score = apply(expr, 2, mean)

info = cbind(tis.score, info)
se = which(info$Visit =="Pre-treatment")
info = info[se,]
se = which(info$Confirmed.Response_IRF %in% c("CR", "PR", "SD", "PD"))
info = info[se, ]


#---------------------------
comxx = intersect(row.names(data), row.names(info))
data = data[comxx,]
info = info[comxx,]
info = info[, 1:9]


#--------------------------
#--------------------------
mys1 = data$TP53__MUT
mys2 = info$tis.score
plot(mys1, mys2)
abline(v=median(mys1))
abline(h=median(mys2))

mytag = rep("", nrow(data))
mytag[mys1>=median(mys1) & mys2>=median(mys2)] = "HH"
mytag[mys1>=median(mys1) & mys2<median(mys2)] = "HL"
mytag[mys1<median(mys1)] = "P53L"

mytag[mys1<median(mys1) & mys2>=median(mys2)] = "LH"
mytag[mys1<median(mys1) & mys2<median(mys2)] = "LL"
table(mytag)


xx = info
xx$mytag = mytag
table(xx$mytag, xx$Confirmed.Response_IRF)

raw.data = xx


#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

data = raw.data
mycat = data$Confirmed.Response_IRF
mys = data$mytag
data = data.frame(mys, mycat)
se = which(data$mycat !="")
data = data[se,]



table(data$mycat)

xx= table(data$mys, data$mycat)
myp = fisher.test(xx, alternative="g")$p.value
myp = ifelse(myp<0.001, formatC(myp, format = "e", digits = 0), signif(myp, 1))
txt.myp = paste("P=", myp, sep="")


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
library(ggplot2)
library(grid)
library(dplyr)

# data$mys <- factor(data$mys, levels = c("HH", "HL", "P53L"))
data$mys <- factor(data$mys, levels = c("HH", "HL","LH", "LL"))
data$mycat <- factor(data$mycat, levels = c("CR", "PR", "SD", "PD"))


# Table of Response × Immune
tab <- table(data$mys, data$mycat)
tab
# fisher.test(tab, workspace = 2e8)


# # Add a label above each bar group
# label_df <- data.frame(
#   mys = 1.5,
#   y = 1.05,         
#   label = "P = 0.008"
# )


plot_data <- data %>%
  group_by(mys, mycat) %>%
  summarise(n = n(), .groups = 'drop') %>%
  group_by(mys) %>%
  mutate(freq = n / sum(n),
         label = paste0(round(freq * 100), "%"))


mygg2<-ggplot(plot_data,aes(y = freq, x = mys, fill = mycat)) +
  geom_flow(aes(alluvium = mycat), alpha= .5, color = "white",curve_type = "linear",   width = .7) +
  geom_col(width = .7, color = "white") +
  labs(title = "Ate mono",x = "TP53 score", y = "Proportion", fill = "Response") +
  scale_fill_manual(values = mycol) + 
  scale_y_continuous(NULL, expand = c(0,0)) +
  # theme_bw() +
  geom_text(
    aes(label = label),
    position = position_stack(vjust = 0.5),  # center inside each segment
    size = 4,
    # fontface = "bold",
    color = "black"
  )+
# +  geom_text(data = label_df, aes(x = mys, y = y, label = label),
#                  inherit.aes = FALSE,
#                  vjust = 0, hjust = 0.5,
#                  size = 4)  +
  theme_classic()+
  theme(
    axis.text.x = element_text(size = 12),       
    axis.text.y = element_text(size = 12),       
    axis.title.x = element_text(size = 14),      
    axis.title.y = element_text(size = 14),      
    plot.title = element_text(size = 16, face = "bold")
    # ,
    # legend.position = "none"
  )+ylim(c(NA,1.05))



# posx = 0.5
# posy = 1
# mytxt <- grobTree(textGrob(txt.myp, x=posx,  y=posy, hjust=0, vjust=1, gp=gpar(col="black", fontsize=14)))
# mygg2 <- mygg2 + annotation_custom(mytxt)
# mygg2


graphics.off()
pdf(file = paste0(FigDir,"5B_3.pdf"), width = 4.5, height = 3.5) 
mygg2
dev.off()





data$mys <- factor(data$mys, levels = c("HH", "HL", "LH", "LL"))
data$mycat <- factor(data$mycat, levels = c("CR", "PR", "SD", "PD"))

data$mys 
data$mycat


tab <- table(data$mycat, data$mys )
tab


row_tot <- rowSums(tab)
col_tot <- colSums(tab)
N <- sum(tab)
exp <- outer(row_tot, col_tot) / N
oe <- (tab+1) / (exp+1)

pvals <- matrix(NA, nrow=nrow(tab), ncol=ncol(tab))

for(i in 1:nrow(tab)){
  for(j in 1:ncol(tab)){
    
    a <- tab[i,j]
    b <- sum(tab[i,-j])
    c <- sum(tab[-i,j])
    d <- sum(tab[-i,-j])
    
    mat <- matrix(c(a,b,c,d), nrow=2)
    
    pvals[i,j] <- fisher.test(mat)$p.value
  }
}

rownames(pvals) <- rownames(tab)
colnames(pvals) <- colnames(tab)
plot_df <- data.frame(
  Response = rep(rownames(tab), times=ncol(tab)),
  Group = rep(colnames(tab), each=nrow(tab)),
  Observed = as.vector(tab),
  Expected = as.vector(exp),
  OE = as.vector(oe),
  Pvalue = as.vector(pvals)
)

plot_df$logP <- -log10(plot_df$Pvalue)

plot_df$signedP <- plot_df$logP * sign(log2(plot_df$OE))

plot_df$Group <- factor(plot_df$Group,
                        levels=c("HH", "HL", "LH", "LL"))
plot_df$Response <- factor(plot_df$Response, levels = c("CR","PR","SD","PD"))

plot_df$color_val <- ifelse(plot_df$Pvalue > 0.2, NA, plot_df$signedP)
 plot_df$color_val <- plot_df$signedP

 plot_df$logOE <- log2(plot_df$OE)
 
 
 P<-ggplot(plot_df, aes(x = Response, y = Group)) +
   geom_point(aes(color = logOE, size = logP)) +
   labs(title='Ate mono',x = "Response", y = "TP53 X TIS") +
   scale_size_continuous(
     name = "-log10(P)",
     range = c(3, 8),      # min and max point size
     breaks = scales::pretty_breaks(n = 3)
   ) +
   scale_color_gradientn(
     # colors = c("steelblue4","lightcyan","white","rosybrown2","#B40426"),  
     colors = c("steelblue4","#8DB0FE","white","indianred2","#B40426"),
     values = scales::rescale(c(min(plot_df$logOE), 0, max(plot_df$logOE))),
     limits = c(min(plot_df$logOE), max(plot_df$logOE)),
     name = "log2(OE)"
   ) +
   
   theme_bw()+
   theme(
     axis.text.x = element_text(size = 12,angle = 90),       
     axis.text.y = element_text(size = 12),       
     axis.title.x = element_text(size = 14),      
     axis.title.y = element_text(size = 14),      
     plot.title = element_text(size = 16, face = "bold")  
   )
 
 
 


graphics.off()
pdf(file = paste0(FigDir,"dot_tis_score_4cat_mono.pdf"), width = 3.3, height = 3.5) 
P
dev.off()




#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
## Fig0_Boxplot_TP53_4Cat_Sora  -- sorafenib only
# rm(list=ls())
myinf1 = "./data/Zhu_GO30140_IMbrave150_sorafenib.rda"
myinf2 = "./data/Zhu_GO30140_IMbrave150_sorafenib__GenomicEvent_iRAS.txt"

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
se = which(info$Confirmed.Response_IRF %in% c("CR", "PR", "SD", "PD"))
info = info[se, ]

comxx = intersect(row.names(data), row.names(info))
data = data[comxx,]
info = info[comxx,]
raw.data = data
raw.info = info

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

data = raw.data
info = raw.info
mycat = info$Confirmed.Response_IRF
mys = data$TP53__MUT
data = data.frame(mys, mycat)
se = which(data$mycat !="")
data = data[se,]

myf = data$mycat
xx = table(myf)
myList = list(NULL)
for(k in 1:length(xx))
{
  se = which(myf==names(xx)[k])
  myList[[k]] = data$mys[se]
}
names(myList) = names(xx)

boxplot(myList)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
library(ggplot2)
library(grid)

data$mycat <- factor(data$mycat, levels = c("CR", "PR", "SD", "PD"))


my_comparisons <- list(
  c("PR", "SD"),
  c("PD", "SD"),
  c("PD", "PR"))

mygg <- ggplot(data, aes(x=mycat, y=mys,  color=mycat)) + 
  # geom_boxplot(size=0.7,outlier.size = 0.3,linewidth=0.5) +
  geom_boxplot(size=1, outlier.shape = NA) +
  # geom_jitter(width = 0.15, size = 0.3, alpha = 0.7) +
  scale_fill_manual(values="white") +
  scale_color_manual(values=mycol) +
  labs(title='sorafenib',x= 'Response', y = 'TP53 score')+
  theme_classic() + 
  theme(legend.position="none") +
  theme(
    axis.text.x = element_text(size = 12),       
    axis.text.y = element_text(size = 12),       
    axis.title.x = element_text(size = 14),      
    axis.title.y = element_text(size = 14),      
    plot.title = element_text(size = 16, face = "bold")  
  )+
  stat_compare_means(
    comparisons = my_comparisons,
    method = "wilcox.test",      
    label = "p.signif",
    size = 4.5   # Increase size of p-values
  )+coord_cartesian(ylim = c(NA, 4.3))


graphics.off()
pdf(paste0(FigDir,"5D.pdf"), width = 2.5, height = 3.5) 
mygg
dev.off()

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

data = raw.data
info = raw.info
mycat = info$Confirmed.Response_IRF
mys = ifelse(data$TP53__MUT>median(data$TP53__MUT), "High", "Low")
data = data.frame(mys, mycat)
se = which(data$mycat !="")
data = data[se,]

xx= table(data$mys, data$mycat)
myp = fisher.test(xx, alternative="g")$p.value
myp = ifelse(myp<0.001, formatC(myp, format = "e", digits = 0), signif(myp, 1))
txt.myp = paste("P=", myp, sep="")


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
library(ggplot2)
library(grid)
library(dplyr)

data$mys <- factor(data$mys, levels = c( "High", "Low"))
data$mycat <- factor(data$mycat, levels = c("CR", "PR", "SD", "PD"))

plot_data <- data %>%
  group_by(mys, mycat) %>%
  summarise(n = n(), .groups = 'drop') %>%
  group_by(mys) %>%
  mutate(freq = n / sum(n),
         label = paste0(round(freq * 100), "%"))

# Table of Response × Immune
tab <- table(data$mys, data$mycat)
tab
fisher.test(tab, workspace = 2e8)


label_df <- data.frame(
  mys = 1.5,
  y = 1.05,         
  label = "P = 5e-04"
)

mygg1<-ggplot(plot_data,aes(y = freq, x = mys, fill = mycat)) +
  geom_flow(aes(alluvium = mycat), alpha= .5, color = "white",curve_type = "linear",   width = .7) +
  geom_col(width = .7, color = "white") +
  labs(title="sorafenib",x = "TP53 score", y = "Proportion", fill = "Response") +
  scale_fill_manual(values = mycol) + 
  scale_y_continuous(NULL, expand = c(0,0)) +
  geom_text(
    aes(label = label),
    position = position_stack(vjust = 0.5),  # center inside each segment
    size = 4,
    # fontface = "bold",
    color = "black"
  ) +  geom_text(data = label_df, aes(x = mys, y = y, label = label),
                 inherit.aes = FALSE,
                 vjust = 0, hjust = 0.5,
                 size = 4)  +
  theme_classic()+
  theme(
    axis.text.x = element_text(size = 12),       
    axis.text.y = element_text(size = 12),       
    axis.title.x = element_text(size = 14),      
    axis.title.y = element_text(size = 14),      
    plot.title = element_text(size = 16, face = "bold")
    # ,
    # legend.position = "none"
  )+ylim(c(NA,1.05))





graphics.off()
pdf(file =paste0(FigDir,"5E.pdf"), width = 3.2, height = 3.5) 
mygg1
dev.off()






################
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
## -- Sorafinib -- combine TP53 with TLS -- 3 cat
# rm(list=ls())
myinf1 = "./data/Zhu_GO30140_IMbrave150_sorafenib.rda"
myinf2 = "./data/Zhu_GO30140_IMbrave150_sorafenib__GenomicEvent_iRAS.txt"


#---------------------------
data <- read.table(myinf2, header=T, sep="\t", row.names=1, quote="")
cnum = ncol(data)/2
data = data[, 1:cnum]
tmp = colnames(data)
tmp = gsub("\\.ES", "", tmp)
colnames(data) = tmp
cnum = ncol(data)/2
dat1 = data[,1:cnum]
dat2 = data[, (cnum+1):(2*cnum)]
xx = dat1-dat2
colnames(xx) = gsub("\\.up", "", colnames(dat1))
data = xx
for(k in 1:ncol(data))
{
  data[,k] = data[,k]/sd(abs(data[,k]))
}
data = data[grep("mul.adj__", colnames(data))]
colnames(data) = gsub("mul.adj__", "", colnames(data))


#---------------------------
load(file= myinf1)
info = info	##  samples

expr = log2(rna+1)
tis_genes = c('CCL5', 'CD27', 'CD274', 'CD276', 'CD8A', 'CMKLR1', 'CXCL9', 'CXCR6',
              'HLA-DQA1', 'HLA-DRB1', 'IDO1', 'LAG3', 'NKG7', 'PDCD1LG2', 'PSMB10',
              'STAT1', 'TIGIT', 'TNFRSF9')  ## PMID: 28650338

se = which(row.names(expr) %in% tis_genes)
expr = expr[se,]
tis.score = apply(expr, 2, mean)

info = cbind(tis.score, info)
se = which(info$Visit =="Pre-treatment")
info = info[se,]
se = which(info$Confirmed.Response_IRF %in% c("CR", "PR", "SD", "PD"))
info = info[se, ]


#---------------------------
comxx = intersect(row.names(data), row.names(info))
data = data[comxx,]
info = info[comxx,]
info = info[, 1:9]


#--------------------------
#--------------------------
mys1 = data$TP53__MUT
mys2 = info$tis.score
plot(mys1, mys2)
abline(v=median(mys1))
abline(h=median(mys2))

mytag = rep("", nrow(data))
mytag[mys1>=median(mys1) & mys2>=median(mys2)] = "HH"
mytag[mys1>=median(mys1) & mys2<median(mys2)] = "HL"
mytag[mys1<median(mys1) & mys2>=median(mys2)] = "LH"
mytag[mys1<median(mys1) & mys2<median(mys2)] = "LL"

# mytag[mys1<median(mys1)] = "P53L"
table(mytag)


xx = info
xx$mytag = mytag
table(xx$mytag, xx$Confirmed.Response_IRF)

raw.data = xx


#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
data = raw.data
mycat = data$Confirmed.Response_IRF
mys = data$mytag
data = data.frame(mys, mycat)
se = which(data$mycat !="")
data = data[se,]



table(data$mycat)

xx= table(data$mys, data$mycat)
myp = fisher.test(xx, alternative="g")$p.value
myp = ifelse(myp<0.001, formatC(myp, format = "e", digits = 0), signif(myp, 1))
txt.myp = paste("P=", myp, sep="")


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
library(ggplot2)
library(grid)
library(dplyr)

data$mys <- factor(data$mys, levels = c("HH", "HL", "LH", "LL"))
data$mycat <- factor(data$mycat, levels = c("CR", "PR", "SD", "PD"))

plot_data <- data %>%
  group_by(mys, mycat) %>%
  summarise(n = n(), .groups = 'drop') %>%
  group_by(mys) %>%
  mutate(freq = n / sum(n),
         label = paste0(round(freq * 100), "%"))



# Table of Response × Immune
tab <- table(data$mys, data$mycat)
tab
fisher.test(tab, workspace = 2e8)


label_df <- data.frame(
  mys = 1.5,
  y = 1.05,
  label = "P = 0.003"
)

mygg2<-ggplot(plot_data,aes(y = freq, x = mys, fill = mycat)) +
  geom_flow(aes(alluvium = mycat), alpha= .5, color = "white",curve_type = "linear",   width = .7) +
  geom_col(width = .7, color = "white") +
  labs(title="sorafenib",x = "TP53 score", y = "Proportion", fill = "Response") +
  scale_fill_manual(values = mycol) + 
  scale_y_continuous(NULL, expand = c(0,0)) +
  geom_text(
    aes(label = label),
    position = position_stack(vjust = 0.5),  # center inside each segment
    size = 4,
    # fontface = "bold",
    color = "black"
  ) +  
  geom_text(data = label_df, aes(x = mys, y = y, label = label),
                 inherit.aes = FALSE,
                 vjust = 0, hjust = 0.5,
                 size = 4)  +
  theme_classic()+
  theme(
    axis.text.x = element_text(size = 12),       
    axis.text.y = element_text(size = 12),       
    axis.title.x = element_text(size = 14),      
    axis.title.y = element_text(size = 14),      
    plot.title = element_text(size = 16, face = "bold")
    # ,
    # legend.position = "none"
  )+ylim(c(NA,1.05))

graphics.off()
pdf(file =paste0(FigDir,"5F.pdf"), width = 4.5, height = 3.5) 
mygg2
dev.off()







data$mys <- factor(data$mys, levels = c("HH", "HL", "LH", "LL"))
data$mycat <- factor(data$mycat, levels = c("CR", "PR", "SD", "PD"))

data$mys 
data$mycat


tab <- table(data$mycat, data$mys )
tab


row_tot <- rowSums(tab)
col_tot <- colSums(tab)
N <- sum(tab)
exp <- outer(row_tot, col_tot) / N
oe <- (tab+1) / (exp+1)

pvals <- matrix(NA, nrow=nrow(tab), ncol=ncol(tab))

for(i in 1:nrow(tab)){
  for(j in 1:ncol(tab)){
    
    a <- tab[i,j]
    b <- sum(tab[i,-j])
    c <- sum(tab[-i,j])
    d <- sum(tab[-i,-j])
    
    mat <- matrix(c(a,b,c,d), nrow=2)
    
    pvals[i,j] <- fisher.test(mat)$p.value
  }
}

rownames(pvals) <- rownames(tab)
colnames(pvals) <- colnames(tab)
plot_df <- data.frame(
  Response = rep(rownames(tab), times=ncol(tab)),
  Group = rep(colnames(tab), each=nrow(tab)),
  Observed = as.vector(tab),
  Expected = as.vector(exp),
  OE = as.vector(oe),
  Pvalue = as.vector(pvals)
)

plot_df$logP <- -log10(plot_df$Pvalue)

plot_df$signedP <- plot_df$logP * sign(log2(plot_df$OE))

plot_df$Group <- factor(plot_df$Group,
                        levels=c("HH", "HL", "LH", "LL"))
plot_df$Response <- factor(plot_df$Response, levels = c("CR","PR","SD","PD"))

plot_df$color_val <- ifelse(plot_df$Pvalue > 0.2, NA, plot_df$signedP)
plot_df$color_val <- plot_df$signedP

plot_df<-plot_df[-which(plot_df$Response=="CR"),]

plot_df$logOE <- log2(plot_df$OE)

 
P<-ggplot(plot_df, aes(x = Response, y = Group)) +
  geom_point(aes(color = logOE, size = logP)) +
  labs(title='sorafenib (n=48)',x = "Response", y = "TP53 X TIS") +
  scale_size_continuous(
    name = "-log10(P)",
    range = c(3, 8),      # min and max point size
    breaks = scales::pretty_breaks(n = 3)
  ) +
  scale_color_gradientn(
    # colors = c("steelblue4","lightcyan","white","rosybrown2","#B40426"),  
    colors = c("steelblue4","#8DB0FE","white","indianred2","#B40426"),
    values = scales::rescale(c(min(plot_df$logOE), 0, max(plot_df$logOE))),
    limits = c(min(plot_df$logOE), max(plot_df$logOE)),
    name = "log2(OE)"
  ) +
  
  theme_bw()+
  theme(
    axis.text.x = element_text(size = 12,angle = 90),       
    axis.text.y = element_text(size = 12),       
    axis.title.x = element_text(size = 14),      
    axis.title.y = element_text(size = 14),      
    plot.title = element_text(size = 16, face = "bold")  
  )



graphics.off()
pdf(file = paste0(FigDir,"dot_tis_score_4cat_targeted.pdf"), width = 3.3, height = 3.5) 
P
dev.off()






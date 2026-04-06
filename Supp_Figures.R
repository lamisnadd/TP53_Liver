#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
# Fig.S
# Author: Chao Cheng; Lamis Naddaf
# 04/03/2026
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# rm(list=ls())
FigDir="./Supp_Figures/"
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

###################################
#####################################################
###############

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


####################################
myinf1 = "./data/Fig_Data_4Datasets.rda"
## data.List, clin.List

load(myinf1)

library(survival)
library(survminer)




######
RFsurvival.plot<-function(data.List,clin.List,dataset){
  library(survival)
  library(survminer)
  
  se = which(names(data.List)==dataset)
  #se = which(names(data.List)=="NCI")
  # se = which(names(data.List)=="UTSMC")
  
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
  data$time = data$RFS.time
  data$event = data$RFS.event
  
  fit <- survfit(Surv(time, event) ~ mycat, data = data)
  
  mygg <- ggsurvplot(
    fit,
    data = data,
    risk.table = F,               # adds number at risk
    pval = TRUE,                     
    pval.coord = c(min(data$time), 0.05),
    # adds p-value
    conf.int = FALSE,                 # confidence interval
    legend.title = "TP53-Mut score",
    legend.labs = c("High", "Low"),
    # palette = c("red3", "blue3"),
    palette = c("coral1", "steelblue"),
    xlab = "RFS(months)",
    ylab = "Survival Probability",
    # break.time.by = 10,              # interval on x-axis
    risk.table.height = 0.25,       # adjust table height
    risk.table.y.text.col = TRUE,
    risk.table.y.text = FALSE,
    ggtheme = theme_classic() +
      theme(
        axis.title.x = element_text(size = 14),		# face = "bold"
        axis.title.y = element_text(size = 14),
        axis.text.x  = element_text(size = 12),
        axis.text.y  = element_text(size = 12),
        plot.title   = element_text(size = 16, hjust = 0.5),  # center title
        legend.key.height = unit(0.05, "cm")
      ),
    title = dataset
  )
  
  return(mygg )}

P_UTSMC_RFS<-RFsurvival.plot(data.List,clin.List,"UTSMC")$plot

graphics.off()
pdf(paste0(FigDir,"1Supp.pdf"),width =4,height = 3.5)
print(P_UTSMC_RFS)
dev.off()




###################################################
####################################################
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
raw.data = data
raw.info = info

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

data = raw.data
info = raw.info



dat <- data.frame(
  P53 = data$TP53__MUT,
  Resp = info$Confirmed.Response_IRF,
  OS.time = info$PFS..in.days.IRF/30,
  OS.event = ifelse(info$PFS.censoring..1.cens.0.evt..IRF == 0, 1, 0)  
)


# filter response categories
dat = dat[dat$Resp %in% c("CR", "PR", "SD", "PD"), ]

# High–Low TP53 grouping
dat$P53.group = ifelse(dat$P53 > median(dat$P53, na.rm=TRUE), "High", "Low")


# ============================================================
# KM plotting function for one category
# ============================================================

km_plot <- function(cat, dat,Title) {
  
  sub = dat[dat$Resp %in% cat, ]
  if(nrow(sub) < 5) {
    return(ggplot() + ggtitle(paste(cat, "(N too small)")))
  }
  
  fit = survfit(Surv(OS.time, OS.event) ~ P53.group, data=sub)
  
  
  ggsurvplot(
    fit,
    data = sub,
    pval = TRUE,
    title = Title,
    legend.title = "TP53 score",
    legend.labs = c("High", "Low"),
    palette =c("coral1","steelblue"),
    risk.table = FALSE,     # no risk table in subplots
    xlab = "PFS (months)",
    ylab = "Survival probability",
    censor.shape = "|",
    censor.size = 2,
    ggtheme = theme_classic(base_size = 12)+
      theme(
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        axis.text.x  = element_text(size = 12),
        axis.text.y  = element_text(size = 12),
        plot.title = element_text(size = 16, face = "bold")    # center title
      ))$plot  # extract only the KM curve
  
}


# ============================================================
# Generate each subplot
# ============================================================


p_ALL_PFS = km_plot(c("CR","PR","SD","PD"), dat,"ALL response categories")
fig3g<-p_PD

graphics.off()
pdf(paste0(FigDir,"3D2_2.pdf"),height = 4,width = 3)
print(p_ALL_PFS)
dev.off()



############################################
###############################################################






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
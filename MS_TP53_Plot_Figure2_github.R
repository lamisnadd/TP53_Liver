#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
# Fig.2
# Author: Chao Cheng; Lamis Naddaf
# 04/03/2026
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# rm(list=ls())
FigDir<-"./Figure2/"
myinf1 = "./data/Fig_Data_4Datasets.rda"
## data.List, clin.List

load(myinf1)

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


library(survival)
library(survminer)


dataset="RIKEN"

survival.plot<-function(data.List,clin.List,dataset){
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
  data$time = data$OS.time
  data$event = data$OS.event
  
  
  
  fit <- survfit(Surv(time, event) ~ mycat, data = data)
  
  
  mygg <- ggsurvplot(
    fit,
    data = data,
    risk.table = F,               # adds number at risk
    pval = TRUE,                     
    pval.coord = c(min(data$time), 0.05),
    # adds p-value
    # conf.int = TRUE,                 # confidence interval
    conf.int = FALSE,    
    legend.title = "TP53 score",
    legend.labs = c("High", "Low"),
    # palette = c("red3", "blue3"),
    # palette = c("red3", "grey40"),
    # palette = c("firebrick4", "mediumseagreen"),
    # palette = c("firebrick3", "seagreen3"),
    palette = c("coral1", "steelblue"),
    xlab = "OS(months)",
    ylab = "Survival Probability",
    break.time.by = 10,              # interval on x-axis
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
  
  data$mycat <- relevel(data$mycat, ref = "Low")
  
  cox <- coxph(Surv(time, event) ~ mycat, data = data)
  cox_sum <- summary(cox)
  
  HR  <- round(cox_sum$coefficients[,"exp(coef)"], 2)
  hr_label <- paste0(
    "HR = ", HR )
  
  counts <- length(data$score)
  
  count_label <- paste0("n = ", counts,
                        collapse = "\n"
  )
  mygg$plot <- mygg$plot +
    annotate("text", x = min(data$time), y = 0.15,label = hr_label,size = 5, hjust = 0) +
    annotate( "text", x = min(data$time), y = 0.25,label = count_label,size = 5, hjust = 0)
  return(mygg )}


P_RIKEN<-survival.plot(data.List,clin.List,"RIKEN")
P_TCGA<-survival.plot(data.List,clin.List,"TCGA")
P_NCI<-survival.plot(data.List,clin.List,"NCI")
P_UTSMC<-survival.plot(data.List,clin.List,"UTSMC")

fig1e<-P_RIKEN

graphics.off()
pdf(paste0(FigDir,"1E.pdf"),width =4,height = 3.5)
print(fig1e$plot)
dev.off()



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

P_NCI_RFS<-RFsurvival.plot(data.List,clin.List,"NCI")
P_UTSMC_RFS<-RFsurvival.plot(data.List,clin.List,"UTSMC")$plot

graphics.off()
pdf(paste0(FigDir,"1Supp.pdf"),width =4,height = 3.5)
print(P_UTSMC_RFS)
dev.off()


# data$mycat <- factor(data$mycat, levels = c("High", "Low"))
data$time = data$OS.time
data$event = data$OS.event



##### Fig1c: Forest : score predicting hazards in differnt stages in Jabanese cohorot 

myinf1 = "./data/Fig_Data_4Datasets.rda"
## data.List, clin.List

load(myinf1)
se = which(names(data.List)=="RIKEN")
data = data.List[[se]]
info = clin.List[[se]]

xx = info$donor_age_at_diagnosis
info$age = xx
xx = info$donor_tumour_stage_at_diagnosis
info$stage = as.factor(xx)
xx = info$donor_sex
info$sex = xx
info$time = info$OS.time
info$event = info$OS.event

se = which(colnames(data)=="TP53__MUT")
mytf = as.numeric(data[,se])
mydat = cbind(mytf, info)

library(survival)

mydatHigh<-mydat[which(mydat$mytf>median(mydat$mytf)),]
mydatLow<-mydat[which(mydat$mytf<median(mydat$mytf)),]

mycox = coxph(Surv(time, event)~mytf + age + sex + stage, mydat) 
mycox1 = coxph(Surv(time, event)~stage, mydatHigh)
mycox2 = coxph(Surv(time, event)~stage, mydatLow)
mycox = summary(mycox)

PV1 = mycox$coefficient[,5]
HR1 = mycox$conf.int[,1]
HR.lo = mycox$conf.int[,3]
HR.hi = mycox$conf.int[,4]

raw.data = data.frame(PV1, HR1, HR.lo, HR.hi)

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

res = raw.data
colnames(res) = c("PV1", "HR1", "HR.lo", "HR.hi")
res

PV1        HR1     HR.lo       HR.hi
mytf    0.0000429106  1.4778358 1.2256433   1.7819202
age     0.9678234595  1.0007137 0.9666179   1.0360122
sexmale 0.0174090406  0.4512674 0.2342127   0.8694758
stage2  0.0524004284  7.3825891 0.9794637  55.6453688
stage3  0.0185723689 11.6834576 1.5090464  90.4565858
stage4  0.0054749163 19.2944737 2.3900937 155.7582109

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
library(tidyverse)
library(forestplot)

ids2 = row.names(res)
ids2 = c("TP53 score", "Age", "Sex (Male vs. Female)", "Stage (II vs. I)", "Stage (III vs. I)", "Stage (IV vs. I)")		

mhr = res$HR1
mup = res$HR.hi
mlw = res$HR.lo
myp = res$PV1
myp = ifelse(myp<0.001, formatC(myp, format = "e", digits = 0), signif(myp, 1))
hr = paste(formatC(mhr, 2, format = 'f'), ' (', formatC(mlw, 2, format = 'f') ,' to ',formatC(mup, 2, format = 'f'),')' , sep="")

data = tibble(Feature = ids2, mean  = mhr, lower = mlw, upper = mup,  hr=hr, coxphP = myp)
header <- tibble(Feature = c("Feature"), hr='HR (95%CI)',coxphP = c('P-value'))
data = bind_rows(header, data)
max(mhr)
max(mup)
min(mlw)
data

limits = c(0.2, 5)  
xticks <- c(0.5,  0, 2, 4)  
xtlab <- rep(c(TRUE), length.out = length(xticks))
attr(xticks, "labels") <- xtlab

data %>%  
  forestplot(labeltext = c(Feature, hr,coxphP),
             graph.pos = 2,
             zero = 1,
             is.summary=c(TRUE, rep(FALSE, nrow(data)-1)), 
             xlog = TRUE,
             # clip = limits,
             xlab = "Hazard ratio",
             ci.vertices = TRUE,
             hrzl_lines= TRUE,
             boxsize = 0.2, 
             line.margin = 0.1,
             mar = unit(rep(0.1, times = 4), "mm"),
             # xticks = xticks,
             graphwidth = "auto",    	## unit(5, "cm")
             graphhight = "auto",
             align = rep("l", 3),
             txt_gp=fpTxtGp(label = gpar(cex = 0.9),
                            title = gpar(cex = 0.9),
                            ticks = gpar(cex = 0.9),
                            xlab = gpar(cex = 0.9)),
             col = fpColors(box = "black",
                            line = "black")) |>
  fp_set_zebra_style("#EFEFEF") ->P2

fig1f<-P2
cairo_pdf(paste0(FigDir,"1F.pdf"), width = 5, height = 2.5) 
P2
dev.off()



######################## 
### fig1d: compare TP53 score within sages

# rm(list=ls())
myinf1 = "./data/Fig_Data_4Datasets.rda"
## data.List, clin.List

load(myinf1)
se = which(names(data.List)=="RIKEN")
data = data.List[[se]]
info = clin.List[[se]]

xx = info$donor_age_at_diagnosis
info$age = xx
xx = info$donor_tumour_stage_at_diagnosis
info$stage = as.factor(xx)
xx = info$donor_sex
info$sex = xx
info$time = info$OS.time
info$event = info$OS.event

se = which(colnames(data)=="TP53__MUT")
mytf = as.numeric(data[,se])
mydat = cbind(mytf, info)


library(ggplot2)
library(ggpubr)

# Ensure stage is a factor and ordered (optional but improves plot)
mydat$stage <- factor(mydat$stage, 
                      levels = sort(unique(mydat$stage)),
                      ordered = TRUE)


mycol<-c("burlywood2","darkgoldenrod1","sienna2","darkred")

library(ggpubr)
library(dplyr)

# Convert stage to factor with correct order
mydat <- mydat %>%
  mutate(stage = factor(stage, levels = c(1,2,3,4)))

# Define comparisons: Stage 1 vs 3 and Stage 1 vs 4
my_comparisons <- list(
  
  c("1", "4")
)


# Boxplot with two comparisons
p_box_score_stages <- ggplot(mydat, aes(x = stage, y = mytf, color = stage)) +
  # geom_violin(fill = "lightgrey", color = NA, trim = FALSE, width = 1) +
  geom_boxplot(width = 0.7, outlier.shape = NA, fill = NA) +
  geom_jitter(width = 0.15, size = 0.5, alpha = 0.7) +
  # scale_color_manual(values = c("burlywood2","sienna2","red3","darkred")) +
  scale_color_manual(values = mycol) +
  labs( title = "", x = "Tumor Stage",y = "TP53 Score") +
  theme_classic() +
  theme(
    axis.text.x = element_text(size = 12),       
    axis.text.y = element_text(size = 12),       
    axis.title.x = element_text(size = 14),      
    axis.title.y = element_text(size = 14),      
    plot.title = element_text(size = 14, face = "bold"),
    legend.position = "none") + ylim(NA,4)+
  stat_compare_means(
    comparisons = my_comparisons,
    method = "wilcox.test", 
    size = 6,
    label.y = 2.6,
    method.args =  list(alternative = "less"),
    label = "p.signif"
  )


fig1g<-p_box_score_stages


graphics.off()
pdf(paste0(FigDir,"1G.pdf"),width = 3,height = 2.8)
print(fig1g)
dev.off()

##################
#### fig1e: compare survival within stages based on P53
library(survival)
library(survminer)
library(dplyr)

# ---------------------------------------------
# 1. Prepare metadata
# ---------------------------------------------
mydat <- cbind(mytf = mytf, info)

mydat$TP53_group <- ifelse(mydat$mytf > median(mydat$mytf, na.rm = TRUE),
                           "High", "Low")
mydat$TP53_group <- factor(mydat$TP53_group, levels = c("Low", "High"))

mydat$stage <- as.factor(mydat$stage)

# Remove samples missing survival or stage
mydat <- mydat %>% 
  filter(!is.na(time), !is.na(event), !is.na(stage))

# ---------------------------------------------
# 2. Fit KM model: TP53 High vs Low within each stage
# ---------------------------------------------
fit_stage <- survfit(Surv(time, event) ~ TP53_group + stage, data = mydat)

# ---------------------------------------------
# 3. Plot: faceted KM curves per stage
# ---------------------------------------------
ann_df <- mydat %>%
  mutate(TP53 = factor(TP53_group, levels = c("Low", "High"))) %>%
  group_by(stage) %>%
  do({
    fit <- coxph(Surv(time, event) ~ TP53, data = .)
    s   <- summary(fit)
    
    data.frame(
      HR   = round(s$coefficients[,"exp(coef)"], 2),
      nLow  = sum(.$TP53 == "Low"),
      nHigh = sum(.$TP53 == "High")
    )
  }) %>%
  mutate(
    label = paste0(
      "HR = ", HR, "\n",
      "n(Low) = ", nLow, ", n(High) = ", nHigh
    )
  )

P_Survival_stages<-ggsurvplot_facet(
  fit_stage,
  data = mydat,
  facet.by = "stage",
  # palette = c( "seagreen3","firebrick3"),
  palette = c( "steelblue", "coral1"),
  legend.title = "TP53 level",
  legend.labs = c("Low", "High"),
  risk.table = TRUE,
  risk.table.height = 0.25,   # ← important
  pval = TRUE,
  # pval.coord = c(0, 0.95),
  pval.size = 9, 
  pval.coord = c(0, 0.15),  
  xlab = "Time (months)",
  ylab = "Survival probability",
  ggtheme = theme_classic(base_size = 15)
)
fig1h<-P_Survival_stages$plot



graphics.off()
pdf(paste0(FigDir,"1H.pdf"),height=5,width=5)
print(ggsurvplot_facet(
  fit_stage,
  data = mydat,
  facet.by = "stage",
  # palette = c( "seagreen3","firebrick3"),
  palette = c( "steelblue", "coral1"),
  legend.title = "TP53 level",
  legend.labs = c("Low", "High"),
  risk.table = TRUE,
  risk.table.height = 0.25,   # ← important
  pval = TRUE,
  # pval.coord = c(0, 0.95),
  pval.size = 9, 
  pval.coord = c(0, 0.15),  
  xlab = "Time (months)",
  ylab = "Survival probability",
  ggtheme = theme_classic(base_size = 15) +
    theme(
      panel.border = element_rect(
        color = "black",
        fill = NA,
        linewidth = 0.8)) ))
dev.off()

graphics.off()
pdf(paste0(FigDir,"1H2.pdf"),height=3.2,width=7.5)
print(print(
  ggsurvplot_facet(
    fit_stage,
    data = mydat,
    facet.by = "stage",
    nrow = 1,
    palette = c("steelblue", "coral1"),
    legend.title = "TP53 level",
    legend.labs = c("Low", "High"),
    risk.table = TRUE,
    risk.table.height = 0.25,
    pval = TRUE,
    pval.size = 9,
    pval.coord = c(0, 0.15),
    xlab = "Time (months)",
    ylab = "Survival probability",
    ggtheme = theme_classic(base_size = 15) +
      theme(
        panel.border = element_rect(
          color = "black",
          fill = NA,
          linewidth = 0.8)) ))$plot)
dev.off()


library(dplyr)
library(ggplot2)

count_df <- mydat %>%
  dplyr::count(stage, TP53_group)
fig1I<-ggplot(count_df, aes(x = stage, y = n, fill = TP53_group)) +
  geom_col(position = position_dodge(width = 0.7), width = 0.6) +
  scale_fill_manual(
    values = c("Low" = "steelblue", "High" = "coral1"),
    name = "TP53 score"
  ) +
  labs(
    x = "Stage",
    y = "Number"
  ) +
  theme_classic(base_size = 10) +
  theme(
    axis.text.x = element_text(size = 12),       
    axis.text.y = element_text(size = 12),       
    axis.title.x = element_text(size = 14),      
    axis.title.y = element_text(size = 14),      
    # plot.title = element_text(size = 14, face = "bold"),
    legend.position = "top"
  )

graphics.off()
pdf(paste0(FigDir,"1I.pdf"),width = 2.5,height = 2.2)
print(fig1I)
dev.off()


library(cowplot)

fig_combined <- plot_grid(
  fig1I + theme(
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank()
  ),
  fig1g,
  ncol = 1,
  align = "v",
  rel_heights = c(0.7, 1)   # adjust height ratio
)

library(cowplot)

fig_combined <- plot_grid(
  fig1I +
    theme(
      axis.title.x = element_blank(),
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      plot.margin = ggplot2::margin(t = 5, r = 5, b = 0, l = 5)
    ),
  fig1g +
    theme(
      plot.margin = ggplot2::margin(t = -5, r = 5, b = 5, l = 5)
    ),
  ncol = 1,
  align = "v",
  rel_heights = c(0.6, 1)
)



graphics.off()
pdf(paste0(FigDir,"1_stage_TP53_combined.pdf"), width = 3, height = 4.2)
print(fig_combined)
dev.off()


##### fig1f:  forest comparing mutaion vs TCGA
myinf1 = "./data/Fig_Data_TP53_OS_4Datasets.txt"

data = read.table(myinf1, sep="\t", header=T, row.names=1)
raw.data = data

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

data = raw.data
res = data[-1, 3:ncol(data)]
colnames(res) = c("PV1", "HR1", "HR.lo", "HR.hi")
res

PV1      HR1    HR.lo    HR.hi
TCGA  2.430038e-03 1.386014 1.122336 1.711639
NCI   2.065568e-03 1.333852 1.110482 1.602151
UTSMC 2.743108e-02 1.458788 1.042907 2.040509
RIKEN 2.692581e-06 1.569836 1.300344 1.895179

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
library(tidyverse)
library(forestplot)

ids2 = c("TCGA",  "NCI"  , "UTSMC" ,"RIKEN")


mhr = res$HR1
mup = res$HR.hi
mlw = res$HR.lo
myp = res$PV1
myp = ifelse(myp<0.001, formatC(myp, format = "e", digits = 0), signif(myp, 1))
hr = paste(formatC(mhr, 2, format = 'f'), ' (', formatC(mlw, 2, format = 'f') ,' to ',formatC(mup, 2, format = 'f'),')' , sep="")

data = tibble(Feature = ids2, mean  = mhr, lower = mlw, upper = mup,  hr=hr, coxphP = myp)
header <- tibble(Feature = c("Feature"), hr='HR (95%CI)',coxphP = c('P-value'))
data = bind_rows(header, data)
max(mhr)
max(mup)
min(mlw)
data

limits = c(0.5, 2.5)  
xticks <- c(0.5,  0, 1, 1.5,2)  
xtlab <- rep(c(TRUE), length.out = length(xticks))
attr(xticks, "labels") <- xtlab
# data<-data[c(1,6,5,4,3,2),]

data %>%  
  forestplot(data=data,labeltext = c(Feature, hr,coxphP),
             graph.pos = 2,
             zero = 1,
             is.summary=c(TRUE, rep(FALSE, nrow(data)-1)), 
             xlog = FALSE,
             clip = limits,
             xlab = "Hazard ratio",
             ci.vertices = TRUE,
             hrzl_lines= TRUE,
             boxsize = 0.2, 
             line.margin = 0.1,
             mar = unit(rep(0.1, times = 4), "mm"),
             xticks = xticks,
             graphwidth = "auto",    	## unit(5, "cm")
             graphhight = "auto",
             align = rep("l", 3),
             txt_gp=fpTxtGp(label = gpar(cex = 0.9),
                            title = gpar(cex = 0.9),
                            ticks = gpar(cex = 0.9),
                            xlab = gpar(cex = 0.9)),
             col = fpColors(box = "black",
                            line = "black")) |>
  fp_set_zebra_style("#EFEFEF")-> P1



fig1i<-P1
cairo_pdf(paste0(FigDir,"1J2.pdf"), width = 7, height = 2) 
P1
dev.off()



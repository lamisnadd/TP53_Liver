#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
# Fig.4
# Author: Chao Cheng; Lamis Naddaf
# 04/03/2026
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

rm(list=ls())
FigDir<-"./Figure4/"
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




mycol <-c(
  "CR" = "steelblue4",
  "PR" = "steelblue2", 
  "SD" = "grey60", 
  "PD" = "indianred4")

# 
#   mygg <- ggplot(data, aes(x=mycat, y=mys,  fill=mycat)) + 
#     geom_boxplot(size=1) +
#     scale_fill_manual(values=mycol) +
# Define comparisons: Stage 1 vs 3 and Stage 1 vs 4
my_comparisons <- list(c("CR", "SD"), c("PR", "SD"), c("PD", "SD"), c("PR", "PD"),c("CR", "PD"))

mygg1 <- ggplot(data, aes(x=mycat, y=mys, color=mycat)) + 
  geom_boxplot(size=1,outlier.size = 0.3,linewidth=1) +
  # geom_jitter(width = 0.15, size = 0.3, alpha = 0.7) +
  # scale_fill_manual(values="white") +
  scale_color_manual(values=mycol) +
  
  # scale_fill_manual(values=mycol) +
  labs(title='Ate+Bev', x= 'Response', y = 'TP53 score') +
  theme_classic() + 
  theme(legend.position="none") +
  theme(
    axis.text.x = element_text(size = 12),       
    axis.text.y = element_text(size = 12),       
    axis.title.x = element_text(size = 14),      
    axis.title.y = element_text(size = 14),      
    plot.title = element_text(size = 16, face = "bold")
    # ,
    # plot.margin = margin(10, 10, 10, 10)  # Add margins
  ) +
  stat_compare_means(
    comparisons = my_comparisons,
    method = "wilcox.test",      
    label = "p.signif",
    size = 4.5   # Increase size of p-values
  )+coord_cartesian(ylim = c(NA, 6.6))


fig3a<-mygg1

graphics.off()
pdf(paste0(FigDir,"3A.pdf"), width = 4, height = 4) 
mygg1
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

mycol <- c(
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

# Add a label above each bar group
label_df <- data.frame(
  mys = 0,
  y = 1.1,         
  label = "1e-05"
)

mygg3<-ggplot(plot_data,aes(y = freq, x = mys, fill = mycat)) +
  geom_flow(aes(alluvium = mycat), alpha= .5, color = "white",curve_type = "linear",   width = .7) +
  geom_col(width = .7, color = "white") +
  labs(title='Ate+Bev',x = "TP53 score", y = "Proportion", fill = "Response") +
  scale_fill_manual(values = mycol) + 
  # scale_y_continuous(NULL, expand = c(0,0)) +
  # theme_bw() +
  theme_classic() +   # classic theme
  # cowplot::theme_minimal_hgrid() +
  geom_text(
    aes(label = label),
    position = position_stack(vjust = 0.5),  # center inside each segment
    size = 4,
    fontface = "bold",
    color = "black"
  ) +  

  theme(
    axis.text.x = element_text(size = 12),       
    axis.text.y = element_text(size = 12),       
    axis.title.x = element_text(size = 14),      
    axis.title.y = element_text(size = 14),      
    plot.title = element_text(size = 16, face = "bold")
    # ,
    # legend.position = "none"
  )+ylim(NA, 1.05)
# +
#   theme(panel.grid.major = element_blank(), 
#         axis.text.y = element_blank(), 
#         axis.ticks.y = element_blank())

 





posx = 0.5
posy = 0.95
mytxt <- grobTree(textGrob(txt.myp, x=posx,  y=posy, hjust=0.5, vjust=0.5, gp=gpar(col="black", fontsize=16)))
mygg3 <- mygg3 + annotation_custom(mytxt)
mygg3

fig3b<-mygg3
graphics.off()
pdf(file = paste0(FigDir,"3B.pdf"), width = 3.4, height = 4) 
mygg3
dev.off()


#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

data = raw.data
info = raw.info
mycat = ifelse(info$Confirmed.Response_IRF=="SD", "SD", "Non-SD")
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

data$mycat[is.na(data$mycat)] <- "Non-SD"
data$mycat <- factor(data$mycat, levels = c("Non-SD", "SD"))


plot_data <- data %>%
  group_by(mys, mycat) %>%
  summarise(n = n(), .groups = 'drop') %>%
  group_by(mys) %>%
  mutate(freq = n / sum(n),
         label = paste0(round(freq * 100), "%"))


# Add a label above each bar group
label_df <- data.frame(
  mys = 0,
  y = 1.1,         
  label = "1e-05"
)


# Set the order of categories
plot_data$mycat <- factor(plot_data$mycat, levels = c("Non-SD", "SD"))
library(ggalluvial)

# install.packages("ggpattern")   # if not installed
library(ggpattern)


# Now plot
mycols <- c("SD" = "grey45", "Non-SD" = "white")
library(ggalluvial)



mygg2<-ggplot(plot_data,aes(y = freq, x = mys, fill = mycat)) +
  # ---- PATTERNED COLUMNS ----
geom_col_pattern(
  aes(
    pattern = mycat
  ),
  color = "grey45",        # bar border
  pattern_fill   = "black",   # color of stripes
  pattern_colour = "white",
  pattern_density = 0.7,
  pattern_spacing = 0.01,
  pattern_angle = 45
) +
  scale_pattern_manual(name = "Response", values=c('stripe', 'wave')) +
  # geom_col(width = .7, color = "white") +
  # geom_flow(aes(alluvium = mycat), alpha= .5, color = "white",curve_type = "linear",   width = .7) +
  labs(title='Ate+Bev',x = "TP53 score", y = "Proportion", fill = "Response") +
  scale_fill_manual(values = mycols) +
  guides(fill = "none") + 
  scale_y_continuous(NULL, expand = c(0,0)) +
  # cowplot::theme_minimal_hgrid() +
  theme_classic() +   # classic theme
  geom_text(
    aes(label = label),
    position = position_stack(vjust = 0.5),  # center inside each segment
    size = 4,
    fontface = "bold",
    color = "black"
  ) + 
  theme(
    axis.text.x = element_text(size = 12),       
    axis.text.y = element_text(size = 12),       
    axis.title.x = element_text(size = 14),      
    axis.title.y = element_text(size = 14),      
    plot.title = element_text(size = 16, face = "bold")
    # ,
    # legend.position = "none"
  )+ylim(NA, 1.05)


posx = 0.5
posy = 0.95
mytxt <- grobTree(textGrob(txt.myp, x=posx,  y=posy, hjust=0.5, vjust=0.5, gp=gpar(col="black", fontsize=16)))
mygg2 <- mygg2 + annotation_custom(mytxt)
mygg2

fig3c<-mygg2
graphics.off()
pdf(file = paste0(FigDir,"3C.pdf"), width = 3, height = 4) 
mygg2
dev.off()
##########

# ============================================================
# Prepare data
# ============================================================

data = raw.data
info = raw.info

dat <- data.frame(
  P53 = data$TP53__MUT,
  Resp = info$Confirmed.Response_IRF,
  OS.time = info$OS.in.days/30.0,
  OS.event = ifelse(info$OS.censoring..1.cens.0.evt. == 0, 1, 0)  # 1 = death, 0 = censored
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
    palette = c("coral1","steelblue"),
    risk.table = FALSE,     # no risk table in subplots
    xlab = "OS(months)",
    ylab = "Survival probability",
    censor.shape = "|",
    censor.size = 2,
    ggtheme = theme_classic(base_size = 12)+
      theme(
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        axis.text.x  = element_text(size = 12),
        axis.text.y  = element_text(size = 12),
        plot.title = element_text(size = 16, face = "bold")  # center title
      ))$plot  # extract only the KM curve
}

# ============================================================
# Generate each subplot
# ============================================================
# 
# p_CR = km_plot("CR", dat)
# p_PR = km_plot("PR", dat)
p_SD = km_plot("SD", dat,"SD category")
# p_PD = km_plot("PD", dat)

p_ALL = km_plot(c("CR","PR","SD","PD"), dat,"All categories")
fig3d<-p_ALL

graphics.off()
pdf(paste0(FigDir,"3D.pdf"),height = 4,width = 3)
print(p_ALL)
dev.off()

graphics.off()
pdf(paste0(FigDir,"3D2.pdf"),height = 4,width = 3)
print(p_SD)
dev.off()


dat_SD<-dat[which(dat$Resp=="SD"),]
dat_CR<-dat[which(dat$Resp=="CR"),]
dat_PR<-dat[which(dat$Resp=="PR"),]
dat_PD<-dat[which(dat$Resp=="PD"),]


mycox_SD = coxph(Surv(OS.time, OS.event)~P53, dat_SD)
mycox_SD  = summary(mycox_SD)

mycox_CR = coxph(Surv(OS.time, OS.event)~P53, dat_CR)
mycox_CR  = summary(mycox_CR)

mycox_PR = coxph(Surv(OS.time, OS.event)~P53, dat_PR)
mycox_PR  = summary(mycox_PR)

mycox_PD = coxph(Surv(OS.time, OS.event)~P53, dat_PD)
mycox_PD  = summary(mycox_PD)




PV_SD = mycox_SD$coefficient[,5]
HR1_SD = mycox_SD$conf.int[,1]
HR.lo_SD = mycox_SD$conf.int[,3]
HR.hi_SD = mycox_SD$conf.int[,4]

PV_CR = mycox_CR$coefficient[,5]
HR1_CR = mycox_CR$conf.int[,1]
HR.lo_CR = mycox_CR$conf.int[,3]
HR.hi_CR = mycox_CR$conf.int[,4]

PV_PR = mycox_PR$coefficient[,5]
HR1_PR = mycox_PR$conf.int[,1]
HR.lo_PR = mycox_PR$conf.int[,3]
HR.hi_PR = mycox_PR$conf.int[,4]

PV_PD = mycox_PD$coefficient[,5]
HR1_PD = mycox_PD$conf.int[,1]
HR.lo_PD = mycox_PD$conf.int[,3]
HR.hi_PD = mycox_PD$conf.int[,4]

raw.data = data.frame(PV1=c(PV_SD,PV_CR,PV_PR,PV_PD), HR1=c(HR1_SD,HR1_CR,HR1_PR,HR1_PD), 
                      HR.lo=c(HR.lo_SD,HR.lo_CR,HR.lo_PR,HR.lo_PD), HR.hi=c(HR.hi_SD,HR.hi_CR,HR.hi_PR,HR.hi_PD))
raw.data<-raw.data[-2,]
rownames(raw.data)<-c("SD","PR","PD")
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

res = raw.data
colnames(res) = c("PV1", "HR1", "HR.lo", "HR.hi")
res

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
library(tidyverse)
library(forestplot)

ids2 = row.names(res)

mhr = res$HR1
mup = res$HR.hi
mlw = res$HR.lo
myp = res$PV1
myp = ifelse(myp<0.001, formatC(myp, format = "e", digits = 0), signif(myp, 1))
hr = paste(formatC(mhr, 2, format = 'f'), ' (', formatC(mlw, 2, format = 'f') ,' to ',formatC(mup, 2, format = 'f'),')' , sep="")

data = tibble(Feature = ids2, mean  = mhr, lower = mlw, upper = mup,  hr=hr, coxphP = myp)
header <- tibble(Feature = c("Feature"), hr='HR (95%CI)',coxphP = c('P-value'))
header$coxphP <- as.character(header$coxphP)
data$coxphP <- as.character(data$coxphP)

data_combined <- bind_rows(header, data)

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
cairo_pdf(paste0(FigDir,"1Supp.pdf"), width = 5, height = 2) 
P2
dev.off()




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

p_PD = km_plot("PD", dat,"PD response category")
p_ALL_PFS = km_plot(c("CR","PR","SD","PD"), dat,"ALL response categories")
fig3g<-p_PD

graphics.off()
pdf(paste0(FigDir,"3D2_2.pdf"),height = 4,width = 3)
print(p_ALL_PFS)
dev.off()


##############
#########################
###########
library(ggplot2)
library(clusterProfiler)
library(org.Hs.eg.db)
library(EnhancedVolcano)
library(dplyr)
myinf1 = "./data/Zhu_GO30140_IMbrave150_atezo_bev.rda"
load(file= myinf1)
##### enriched pathways

res<- readRDS("./data/High_P53_CRPR_vs_PD_DE.rds") 
sig_up <- res[res$adj.P.Val < 0.1 & res$logFC > 0.1, ]
rownames(sig_up)

sig_up <- res[res$adj.P.Val < 0.1 & res$logFC < 0.1, ]
rownames(sig_up)

Sig_immune_genes <- c( "ARHGDIG", "MIR181C", "TRGJ2", "TRAJ23", "TRAJ39", "TRAJ45", 
                       "TRGV3", "IGKV2.14", "IGKV1OR.2", "IGKV2.10", "IGLJ6", 
                       "IGHV3OR16.17", "IGHJ3P", "TBX21", "IKZF3", "CD70", "CCL4", "CXCL9", 
                       "IFNG", "NCR1",  "PDCD1", "TIGIT", "FCRL3", 
                       "ICOS", "IL26",  "MS4A15") ## according to chatgpt

genes <- c("PDCD1","TIGIT","ICOS","TBX21","IKZF3","IFNG","IL26","CCL4","CXCL9",
           "NCR1","CD70","FCRL3",
           "IGKV1OR.2","IGLJ6","TRAJ23","TRAJ39","TRAJ45",
           "TRGV3","MUC16","MUC21","ELF5","OVOL1","VGLL1",
           "PROK1","BRS3","NCAN","STOML3","DEGS1","CORO2A","MYH4","ENSG00000223911","CNTFR","TMEM125") ### label genes -acordng to chatgpt

#-----------------------------
# Volcano plot
#-----------------------------

# create custom key-value pairs for 'high', 'low', 'mid' expression by fold-change
# this can be achieved with nested ifelse statements
keyvals <- ifelse(
  res$P.Value > 0.01, 'grey80',                 #
  ifelse(res$logFC < -0.1, 'indianred4',        # PD up
         ifelse(res$logFC > 0.1, 'steelblue3', # CR/PR up
                'grey80')))

# Fix any NAs
keyvals[is.na(keyvals)] <- 'grey80'

# Add names for the legend
names(keyvals)[keyvals == 'indianred4'] <- 'Upregulated\n in PD'
names(keyvals)[keyvals == 'steelblue3'] <- 'Upregulated\n in CR/PR'
names(keyvals)[keyvals == 'grey80'] <- 'NS'

res[1,]

fig4e2_p<-EnhancedVolcano(res,
                          lab = rownames(res),
                          x = 'logFC',
                          y = 'P.Value',
                          title = 'High TP53: CR/PR vs PD',
                          selectLab=union(genes,Sig_immune_genes),
                          pCutoff = 0.01,
                          colCustom = keyvals,
                          labSize = 3,
                          labCol = 'black',
                          labFace = 'bold',
                          # boxedLabels = TRUE,
                          FCcutoff = 0.1,
                          # legendPosition = 'right',
                          colAlpha = 4/5,
                          pointSize = 0.8,
                          drawConnectors = TRUE,
                          widthConnectors = 0.1,
                          colConnectors = 'black',
                          legendLabSize = 10,
                          xlim = c(-1.2, 1.2),
                          ylim = c(NA, 6))


graphics.off()
pdf(paste0(FigDir,"1EE.pdf"),height = 7,width = 5)
print(fig4e2_p)
dev.off()



# Fig0_GroupBoxplot_6Genes_RvsPD_comb
# rm(list=ls())
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


load(file= myinf1)
info = info	##  samples

se = which(info$Visit =="Pre-treatment")
info = info[se,]
table(info$Confirmed.Response_IRF, info$Response)

comxx = intersect(row.names(data), row.names(info))
data = data[comxx,]
info = info[comxx,]
info = info[, 1:9]

score = data$mul.adj__TP53__MUT
info$score = score

se = which(info$score>median(info$score) & info$Confirmed.Response_IRF %in% c("CR", "PR"))
sam1 = row.names(info)[se]
se = which(info$score>median(info$score) & info$Confirmed.Response_IRF %in% c("PD"))
sam2= row.names(info)[se]

#------------------------
data = log2(rna+1)
se = which(colnames(data) %in% sam1)
dat1 = data[,se]
se = which(colnames(data) %in% sam2)
dat2 = data[,se]

subset = c("PDCD1", "CD274", "PDCD1LG2", "ICOS", "IDO1", "IFNG")
dat1 = dat1[subset,]
dat2 = dat2[subset,]

Group = c(rep("CR/PR", ncol(dat1)), rep("PD", ncol(dat2)) )
mydat = rbind(t(dat1), t(dat2) )
mydat = as.data.frame(mydat)
mydat = cbind(Group, mydat)

raw.data = mydat

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
data = raw.data


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Reshape to long format
df_long <- data %>%
  pivot_longer(cols = subset, names_to = "Gene", values_to = "Expression")
df_long$Gene <- factor(df_long$Gene, levels = subset)

df_long<-as.data.frame(df_long)
mygg4 <- ggplot(df_long, aes(x = Group, y = Expression, color = Group)) +
  geom_boxplot(outlier.shape = NA, width = 0.6, linewidth = 0.6) +
  # geom_jitter(width = 0.15, size = 0.5, alpha = 0.6) +
  facet_wrap(~Gene, scales = "free_y", nrow = 2) +
  # theme_classic() +
  
  scale_fill_manual(values = "white") +
  scale_color_manual(values = c("CR/PR" = "steelblue3", "PD" = "indianred3")) +
  labs(title = "", x = "", y = "Expression") +
  theme_classic() +
  stat_compare_means(
    comparisons = list(c("CR/PR", "PD")),
    method = "wilcox.test",
    label = "p.signif",
    size=6
  )+
  ylim(NA, 7.5)+
  theme(
    axis.text.x = element_text(size = 12),       
    axis.text.y = element_text(size = 12),       
    axis.title.x = element_text(size = 14),      
    axis.title.y = element_text(size = 14),      
    plot.title = element_text(size = 16, face = "bold"),
    strip.text = element_text(face = "bold"),
    legend.position = "none"
  )

fig3f<-mygg4

graphics.off()
pdf(file = paste0(FigDir,"3F.pdf"), width = 4.2, height = 4.2) 
mygg4
dev.off()


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

gsea_df_sig$Description
# Sort by normalized enrichment score (NES)
immune_gsea_sig <- immune_gsea_sig[order((immune_gsea_sig$NES)), ]
gsea_df_sig<-gsea_df_sig[order((gsea_df_sig$NES)), ]
# View top results
head(immune_gsea_sig)
gsea_df_sig$Description[which(gsea_df_sig$NES<0)]

library(enrichplot)
# Enrichment plot for top pathway
gseaplot2(gsea_res, geneSetID = 1, title = gsea_res@result$Description[1])

# Dotplot of top enriched terms
dotplot(gsea_res, showCategory = 20)

hist(immune_gsea_sig$NES)
hist(gsea_df_sig$NES)

immune_gsea_sig$Description

selected_pathways <- c(
  "adaptive immune response",                             # overall adaptive immunity
  "T cell activation",                                    # T cell function
  "T cell proliferation",                                  # T cell expansion
  "B cell activation",                                    # B cell function
  "positive regulation of cytokine production",           # immune signaling
  "leukocyte chemotaxis",                                 # migration/trafficking
  "regulation of lymphocyte proliferation",               # lymphocyte control
  "regulation of immune effector process",               # immune effector regulation
  # "innate immune response-activating signaling pathway", # innate immunity activation
  "macrophage activation"    ,                             # myeloid cell function
  "negative regulation of humoral immune response",
  "negative regulation of protein ubiquitination",
  "response to unfolded protein",
  "glutathione metabolic process",
  "oxidative phosphorylation",
  "tricarboxylic acid cycle",
  "negative regulation of coagulation",
  "negative regulation of fibrinolysis",
  # "programmed cell death in response to reactive oxygen species",
  "ribosome assembly")



gsea_df_sig_selected<-gsea_df_sig[which(gsea_df_sig$Description%in%selected_pathways),]


# import package
library("ggplot2")
data <- data.frame("GOs" = gsea_df_sig_selected$Description, 
                   "NES" = gsea_df_sig_selected$NES,
                   "Size" = gsea_df_sig_selected$setSize,
                   "p.adjust" = gsea_df_sig_selected$p.adjust)
# Order GOs based on NES
data$GOs <- factor(data$GOs, levels = data$GOs[order(data$NES)])

# Plot

fig3g<-ggplot(data = data, aes(x = NES, y = GOs, 
                               color = p.adjust, size = Size)) + 
  geom_point() +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey30",linewidth=1.5) + 
  # scale_color_gradientn(colors = c( "salmon3", "lightsalmon2", "lightsalmon","grey")) +
  scale_color_gradientn(colors = c( "brown", "indianred3","grey")) +

  theme_bw() + ylab("") + xlim(-2.6,2.6)+
  xlab("NES") + 
  ggtitle("GSEA (GO)")+
  theme(
    axis.text.x = element_text(size = 12),       
    axis.text.y = element_text(size = 14),       
    axis.title.x = element_text(size = 14),      
    axis.title.y = element_text(size = 14),      
    plot.title = element_text(size = 16, face = "bold")
    # plot.margin = margin(10, 10, 10, 10)  # Add margins
  )

graphics.off()
pdf(paste0(FigDir,"3G.pdf"),height = 4.2,width = 6.5)
print(fig3g)
dev.off()

###############################################################
### 3. IMMUNE CELLS – CR/PR vs PD 
###############################################################
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

myinf1 = "./data/Zhu_GO30140_IMbrave150_atezo_bev_Timer.txt"

data = read.table(myinf1, sep="\t", header=T, row.names=1, check.names=F)

mysam1 = c("A_1", "A_3", "A_5", "A_7", "A_11", "A_12", "A_13", "A_17", "A_19", "A_21", "A_23", "A_24", "A_29", "A_49", "A_51", "A_53", "A_58", "A_60", "A_61", "A_65", "A_69", "A_75", "A_83", "A_100", "A_102", "A_119", "A_121", "A_123", "A_133", "A_141", "A_142", "A_143", "A_147", "A_148", "A_163", "A_175", "A_184", "A_187", "A_191", "A_214", "F1_2", "F1_24", "F1_32", "F1_37", "F1_38", "F1_40", "F1_45")
mysam2 = c("A_20", "A_26", "A_27", "A_28", "A_33", "A_37", "A_38", "A_45", "A_48", "A_62", "A_70", "A_72", "A_84", "A_85", "A_86", "A_116", "A_118", "A_125", "A_126", "A_140", "A_144", "A_153", "A_160", "A_173", "A_176", "A_189", "A_195", "A_200", "A_204", "F1_1", "F1_9", "F1_12", "F1_20", "F1_22", "F1_26", "F1_31")


se = which(row.names(data) %in% mysam1)
dat1 = data[se,]

se = which(row.names(data) %in% mysam2)
dat2 = data[se,]


mycat = c(rep("CR/PR", nrow(dat1)), rep("PD", nrow(dat2)) )
mys = c(as.numeric(dat1[, "T cell CD8+"]),  as.numeric(dat2[, "T cell CD8+"]) )
data = data.frame(mys, mycat)
se = which(data$mycat !="")
data = data[se,]
raw.data = data



table(data$mycat)

xx1 = data$mys[data$mycat=="CR/PR"]
xx2 = data$mys[data$mycat=="PD"]
wilcox.test(xx1, xx2)
boxplot(list(xx1, xx2))

myp = wilcox.test(xx1, xx2, alternative="g")$p.value
myp = ifelse(myp<0.001, formatC(myp, format = "e", digits = 0), signif(myp, 1))
txt.myp = paste("P=", myp, sep="")


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
library(ggplot2)
library(grid)


data$mycat <- factor(data$mycat, levels = c("CR/PR", "PD") )

mycol = c('steelblue3', 'indianred3')[1:2]
mygg5 <- ggplot(data, aes(x=mycat, y=mys,  color=mycat)) + 
  # geom_violin(fill = "lightgrey", color = NA, trim = FALSE, width = 1) +
  # geom_boxplot(outlier.size = 0.5, alpha = 0.7, width = 0.6) +
  geom_boxplot(size=1) +
  # geom_jitter(width = 0.15, size = 0.7, alpha = 0.6) +
  scale_fill_manual(values="white") +
  stat_compare_means(comparisons = list(c("CR/PR","PD")), method = "wilcox.test",label = "p.signif",
                     size = 6 ,label.y = 0.8  ) +
  scale_color_manual(values=mycol) +
  labs(title='',x= 'Response', y = 'CD8+ Tcell')+
  scale_y_continuous(limits = c(0, 0.9)) +
  theme_classic() + 
  theme(legend.position="none") +
  theme(
    axis.text.x = element_text(size = 12),       
    axis.text.y = element_text(size = 12),       
    axis.title.x = element_text(size = 14),      
    axis.title.y = element_text(size = 14),      
    plot.title = element_text(size = 16, face = "bold")  
  )




fig3h<-mygg5
graphics.off()
pdf(file = paste0(FigDir,"3H.pdf"), width = 2, height = 4) 
mygg5
dev.off()


#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
##Immune clusters
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
## determine the immue-cold and Immune-hot  -- Imm cohort
# rm(list=ls())

library(reshape)
library(ggplot2)
library(ggpubr)
library(grid)
library(gridExtra)
library(stringr)
library(ComplexHeatmap)
library(ggsci)
library(circlize)
myinf1 = "./data/Zhu_GO30140_IMbrave150_atezo_bev.rda"

load(file= myinf1)
data = rna	
xx = apply(data>0, 1, sum)
se = which(xx>=10)## >0 in at least 10 samples
data = data[se,]
dim(data)

data = log2(data+1)
data= as.data.frame((data))
myavg = apply(data, 1, mean)
mystd = apply(data, 1, sd)
data = (data-myavg)/mystd


gene_set<-read.table("./data/immune_genesets.txt",sep="\t",header=T,row.names=1)

set_list<-c()
for (i in row.names(gene_set)){
  if (i %in% c("Translation","Proliferation") ){
    next
  }
  set_list<- c(set_list,str_trim(unlist(strsplit(gene_set[i,],","))))
}
library(ComplexHeatmap)
genes<-intersect(set_list,row.names(data))
immune_expression<-data[genes,]
mat<-as.matrix(immune_expression)
scaled_mat = t(scale(t(mat)))
scaled_mat=mat
gene_ss<-data.frame(marker=rep(0,length(genes)))
rownames(gene_ss)<-genes
for ( i in genes){
  for (j in row.names(gene_set)){
    gene_ss[i,"marker"]=ifelse(i %in% str_trim(unlist(strsplit(gene_set[j,],","))),j,0)
    if (gene_ss[i,"marker"] !=0){
      break
    }
  }
}

set.seed(123)
gene_ss$marker<-as.factor(gene_ss$marker)
col_fun = colorRamp2(c(-2, 0, 2), c("#3C5488", "white", "#E64B35"))
lgd=Legend(col_fun = colorRamp2(c(-2, 0, 2), c("#3C5488", "white", "#E64B35")))
pal<-pal_frontiers("default")(8)
names(pal)<-levels(gene_ss$marker)
row_ha = rowAnnotation( df=gene_ss,col=list(marker=row_cols),show_annotation_name=FALSE)
ht<-Heatmap(scaled_mat, right_annotation=row_ha,
            show_column_names = FALSE, cluster_rows = FALSE,name="expression",
            show_row_dend = FALSE,column_km=2, column_gap = unit(5,"mm"),
            show_row_names=FALSE,column_title = c("immune hot", "immune cold"),
            column_title_gp = gpar(fill = c( "#3C5488", "#E64B35")[2:1],col="white"))



graphics.off()
pdf(paste0(FigDir,"ImmuneHeatmap.pdf"),height = 3.5,width=4.5)
draw(ht)
dev.off()

a<-column_dend(ht)
b<-column_order(ht)
clust <- lapply(names(b), function(i){
  out <- data.frame(sampleID = colnames(scaled_mat[,b[[i]]]),
                    Cluster = paste0("cluster", i),
                    stringsAsFactors = FALSE)
  return(out)
}) %>%
  do.call(rbind, .)
xx= table(clust)
apply(xx, 2, sum)


#---


res = data.frame(clust$sampleID, clust$Cluster)
res$cluster.hot = ifelse(res$clust.Cluster=="cluster2", 1, 0)
res$cluster.col = ifelse(res$clust.Cluster=="cluster1", 1, 0)
row.names(res) = res[,1]
res = res[, 3:4]

# write.table(res, myoutf1, sep="\t", quote=F)

apply(res, 2, sum)





#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
## test
# rm(list=ls())
myinf1 = "./data/Zhu_GO30140_IMbrave150_atezo_bev.rda"
myinf2 = "./data/Zhu_GO30140_IMbrave150_atezo_bev__GenomicEvent_iRAS.txt"
myinf3 = "./data/Zhu_GO30140_IMbrave150_AteBev_ImmType.txt"


imm = read.table(myinf3, header=T, sep="\t", row.names=1, quote="")
imm.hot = row.names(imm)[imm$cluster.hot==1]
imm.col = row.names(imm)[imm$cluster.col==1]


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



load(file= myinf1)
info = info	##  samples



expr = log2(rna+1)
tis_genes = c('CCL5', 'CD27', 'CD274', 'CD276', 'CD8A', 'CMKLR1', 'CXCL9', 'CXCR6',
              'HLA-DQA1', 'HLA-DRB1', 'IDO1', 'LAG3', 'NKG7', 'PDCD1LG2', 'PSMB10',
              'STAT1', 'TIGIT', 'TNFRSF9')  ## PMID: 28650338

se = which(row.names(expr) %in% tis_genes)
expr = expr[se,]
tis.score = apply(expr, 2, mean)

se = which(info$Visit =="Pre-treatment")
info = info[se,]
table(info$Confirmed.Response_IRF, info$Response)

comxx = intersect(row.names(data), row.names(info))
data = data[comxx,]
info = info[comxx,]
tis.score = tis.score[comxx]

info = info[, 1:9]

raw.data = data
raw.info = info

#--------------------
#--------------------
score = data$mul.adj__TP53__MUT
info$score = score

info$imm = ifelse(row.names(info)%in%imm.hot, "H", "C")
table(info$imm, info$Response)
table(info$imm, info$Confirmed.Response_IRF)


info$Immune <- ifelse(row.names(info) %in% imm.hot, "Hot", "Cold")
info$Resp2 <- ifelse(info$Response==1, "CR/PR", "Non-responder")
info$Resp4 <- info$Confirmed.Response_IRF



resp2_col = c("CR/PR" = "steelblue3",
  "Non-responder" = "white")


mycol <-c(
  "CR" = "steelblue4",
  "PR" = "steelblue2", 
  "SD" = "grey60", 
  "PD" = "indianred4")

ha = HeatmapAnnotation(
  Response = info[colnames(scaled_mat),"Resp2"],

  col = list(
    Response = resp2_col
  ),
  annotation_label = c("CR/PR")   ,
  annotation_height = unit(3, "mm"),
  show_legend = FALSE
)


top_ha <- HeatmapAnnotation(
  group = anno_block(
    labels = c("Immune hot", "Immune cold"),
    gp = gpar(fill = c( "#E64B35","#3C5488"), col = NA),  # block colors, no border
    labels_gp = gpar(col = "white", fontsize = 12, fontface = "bold")
  )
)


library(circlize)

col_fun <- colorRamp2(
  c(-2, 0, 2),
  c("#3B4CC0", "white", "#B40426")
)

library(RColorBrewer)

# vector of 8 colors
row_cols <- brewer.pal(8, "Set2")
names(row_cols)<-levels(gene_ss$marker)


ht <- Heatmap(
  scaled_mat,
   col = col_fun,
   right_annotation = row_ha,
  top_annotation = top_ha,
  show_column_names = FALSE,
  cluster_rows = FALSE,
  name = "expression",
  bottom_annotation = ha,
  show_row_dend = FALSE,
  column_km = 2,
  column_gap = unit(5,"mm"),
  show_row_names = FALSE
)


graphics.off()
pdf(paste0(FigDir,"ImmuneHeatmap_2.pdf"),height = 4.5,width=5)
draw(ht)
dev.off()

graphics.off()
pdf(paste0(FigDir,"ImmuneHeatmap_22.pdf"),height = 4.5,width=5)
draw(ht)
dev.off()
library(dplyr)
library(ggplot2)
library(cowplot)
library(ggalluvial)


table(info$Confirmed.Response_IRF)
table(info$Immune)
info$score

med <- median(info$score, na.rm = TRUE)

df_bar <- data.frame(
  Response = info$Confirmed.Response_IRF,
  Immune = info$Immune,
  Score_group=ifelse(info$score >= med, "High", "Low"))


df_bar$Immune <- factor(df_bar$Immune, levels = c("Hot", "Cold"))
df_bar$Response <- factor(df_bar$Response, levels = c("CR","PR","SD","PD"))
df_bar$Score_group <- factor(df_bar$Score_group, levels = c("High", "Low"))

# Remove NAs
df_bar <- df_bar[!is.na(df_bar$Response) & !is.na(df_bar$Immune), ]

# Table of Response × Immune
tab <- table(df_bar$Response, df_bar$Immune)
tab
fisher.test(tab)

mycol <-c(
  "CR" = "steelblue4",
  "PR" = "steelblue2", 
  "SD" = "grey60", 
  "PD" = "indianred4")

# Summarize counts and frequencies
plot_data <- df_bar %>%
  group_by(Immune, Response) %>%
  summarise(n = n(), .groups = 'drop') %>%
  group_by(Immune) %>%
  mutate(freq = n / sum(n),
         label = paste0(round(freq * 100), "%"))

label_df <- data.frame(
  mys = 1.5,    # center above Hot/Cold bars
  y = 1.05,     # above the stacked bar
  label = "P = 0.03"
)
mygg3 <- ggplot(plot_data, aes(y = freq, x = Immune, fill = Response)) +
  
  geom_col(width = 0.7, color = "white") +
  # optional flow lines (alluvium)
  geom_flow(aes(alluvium = Response), alpha = 0.5, color = "white", curve_type = "linear", width = 0.7) +
  
  scale_fill_manual(values = mycol) + 
  scale_y_continuous(NULL, expand = c(0,0)) +
  
  labs(title='Ate+Bev',x = "Immune cluster", y = "Propotion", fill = "Response") +
  # segment labels inside bars
  geom_text(aes(label = label),
            position = position_stack(vjust = 0.5),
            size = 4,
            fontface = "bold",
            color = "black") +
  
  # overall p-value label
  geom_text(data = label_df, aes(x = mys, y = y, label = label),
            inherit.aes = FALSE,
            vjust = 0, hjust = 0.5,
            size = 5)  +
  
  # cowplot::theme_minimal_hgrid() +
  theme_classic()+
  theme(
    axis.text.x = element_text(size = 12),       
    axis.text.y = element_text(size = 12),       
    axis.title.x = element_text(size = 14),      
    axis.title.y = element_text(size = 14),      
    plot.title = element_text(size = 16, face = "bold")
    # ,
    # legend.position = "none"
  )+ylim(NA, 1.05)


fig3b<-mygg3
graphics.off()
pdf(file = paste0(FigDir,"response_imm_bar.pdf"), width = 3.2, height = 4) 
mygg3
dev.off()



df_bar$Group <- paste(df_bar$Immune, df_bar$Score_group, sep = "/")

df_bar$Group <- factor(
  df_bar$Group,
  levels = c("Hot/High","Cold/High","Hot/Low","Cold/Low")
)

plot_data <- df_bar %>%
  group_by(Group, Response) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(Group) %>%
  mutate(
    freq = n / sum(n),
    label = paste0(round(freq * 100), "%")
  )


mygg3 <- ggplot(plot_data, aes(x = Group, y = freq, fill = Response)) +
  
  geom_col(width = 0.7, color = "white") +
  
  geom_flow(aes(alluvium = Response),
            alpha = .5,
            color = "white",
            curve_type = "linear",
            width = .7) +
  
  scale_fill_manual(values = mycol) +
  
  scale_y_continuous(NULL, expand = c(0,0)) +
  
  labs(title="Ate+Bev",x = "Immune cluster / TP53 score", y = "", fill = "Response") +
  
  geom_text(aes(label = label),
            position = position_stack(vjust = 0.5),
            size = 4,
            fontface = "bold") +
  
  # cowplot::theme_minimal_hgrid() +
  theme_classic()+
  
  theme(
    axis.text.x = element_text(size = 12),       
    axis.text.y = element_text(size = 0),       
    axis.title.x = element_text(size = 14),      
    axis.title.y = element_text(size = 0),      
    plot.title = element_text(size = 16, face = "bold")
    # ,
    # legend.position = "none"
  )+ylim(NA, 1.05)

tab <- table(df_bar$Response, df_bar$Group)
tab
fisher.test(tab)

# mygg3 <- mygg3 +
#   scale_y_continuous(limits = c(0,1.15), expand = c(0,0))

mygg3



fig3b<-mygg3
graphics.off()
pdf(file = paste0(FigDir,"response_imm_bar_4cat.pdf"), width = 4.5, height = 4) 
mygg3
dev.off()



df_bar$Group <- paste(df_bar$Immune, df_bar$Score_group, sep="/")
df_bar$Group <- factor(df_bar$Group,
                       levels=c("Hot/High","Cold/High","Hot/Low","Cold/Low"))
df_bar$Response <- factor(df_bar$Response, levels = c("CR","PR","SD","PD"))

tab <- table(df_bar$Response, df_bar$Group)
tab


row_tot <- rowSums(tab)
col_tot <- colSums(tab)
N <- sum(tab)
exp <- outer(row_tot, col_tot) / N
oe <- tab / exp

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

names(pvals) <- rownames(tab)
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
plot_df$logOE <- log2(plot_df$OE)
plot_df$signedP <- plot_df$logP * sign(log2(plot_df$OE))

plot_df$Group <- factor(plot_df$Group,
                       levels=c("Hot/High","Cold/High","Hot/Low","Cold/Low"))
plot_df$Response <- factor(plot_df$Response, levels = c("CR","PR","SD","PD"))




P2 <- ggplot(plot_df, aes(x = Response, y = Group)) +
  
  geom_point(aes(color = logOE, size = logP)) +
  
  scale_size_continuous(
    name = "-log10(P)",
    range = c(2, 8),      # min and max point size
    breaks = scales::pretty_breaks(n = 3)
  ) +
  
  scale_color_gradientn(
    colors = c("steelblue4","#8DB0FE","white","indianred3","#B40426"),
    values = scales::rescale(c(min(plot_df$logOE), 0, max(plot_df$logOE))),
    limits = c(min(plot_df$logOE), max(plot_df$logOE)),
    name = "log2(OE)"
  ) +
  
  theme_bw() +
  theme(
    axis.text.x = element_text(size = 12, angle = 90),
    axis.text.y = element_text(size = 12),
    axis.title.x = element_text(size = 14),
    axis.title.y = element_text(size = 0),
    plot.title = element_text(size = 16, face = "bold")
  )

P<-ggplot(plot_df, aes(x = Response, y = Group)) +
  
  geom_point(aes(color = signedP), size = 6, show.legend = TRUE) +  # fixed size here
  
  scale_color_gradientn(
    colors = c("steelblue4","#8DB0FE","white","indianred3","#B40426"),  
    values = scales::rescale(c(min(plot_df$signedP), 0, max(plot_df$signedP))),
    limits = c(min(plot_df$signedP), max(plot_df$signedP)),
    name = "Signed\n-log10(p)"
  ) +
  
  theme_bw()+
  theme(
    axis.text.x = element_text(size = 12,angle = 90),       
    axis.text.y = element_text(size = 12),       
    axis.title.x = element_text(size = 14),      
    axis.title.y = element_text(size = 0),      
    plot.title = element_text(size = 16, face = "bold")  
  )


graphics.off()
pdf(file = paste0(FigDir,"dot_imm_bar_4cat_2.pdf"), width = 3.3, height = 3.5) 
P2
dev.off()

graphics.off()
pdf(file = paste0(FigDir,"dot_imm_bar_4cat.pdf"), width = 3.3, height = 3.5) 
P
dev.off()



info$tis=tis.score
info$imm = ifelse(row.names(info)%in%imm.hot, "H", "C")
table(info$imm, info$Response)
table(info$imm, info$Confirmed.Response_IRF)

fisher.test(table(info$imm, info$Response))

#--------------
xx1 = info$tis[row.names(info)%in%imm.hot]
xx2 = info$tis[row.names(info)%in%imm.col]
boxplot(list(xx1, xx2))
wilcox.test(xx1, xx2)

df <- data.frame(
  score = c(xx1, xx2),
  group = factor(c(rep("Hot", length(xx1)),
                   rep("Cold", length(xx2)))))

my_comparisons <- list(c("Hot", "Cold"))

# Plot with paired significance
P <- ggplot(df, aes(x = group, y = score, color = group)) +
  geom_boxplot(size=1, outlier.shape = NA) +
  # geom_boxplot(outlier.size = 0.5, alpha = 0.7, width = 0.6) +
  # geom_jitter(width = 0.15, size = 0.5, alpha = 0.7) +
  scale_color_manual(values = c("#3C5488", "#E64B35"))  +
  stat_compare_means(
    comparisons = my_comparisons,
    method = "wilcox.test",
    size=4.5,
    # paired = T,
    label = "p.signif"
  )+
  theme_classic() +
  labs(
    title = "Ate+Bev",
    # subtitle = "Each point = mean expression of one gene",
    x = "Immune\ncluster",
    y = "TIS score"
  ) +
  theme(
    axis.text.x = element_text(size = 12),       
    axis.text.y = element_text(size = 12),       
    axis.title.x = element_text(size = 14),      
    axis.title.y = element_text(size = 14),      
    plot.title = element_text(size = 16,face = "bold"),
    legend.position = "none"
  ) + coord_cartesian(ylim = c(NA,6))
graphics.off()
pdf(file = paste0(FigDir,"TP53_tis_Cluster.pdf"), width = 1.8, height = 3) 
P
dev.off()

graphics.off()
pdf(file = paste0("./Figure5/","TP53_tis_Cluster.pdf"), width = 1.8, height = 3.5) 
P
dev.off()




#--------------
xx1 = info$score[row.names(info)%in%imm.hot]
xx2 = info$score[row.names(info)%in%imm.col]
boxplot(list(xx1, xx2))
wilcox.test(xx1, xx2)

df <- data.frame(
  score = c(xx1, xx2)/sd(c(c(xx1, xx2))),
  group = factor(c(rep("Hot", length(xx1)),
                   rep("Cold", length(xx2)))))

my_comparisons <- list(c("Hot", "Cold"))

median_score <- median(df$score, na.rm = TRUE)

# 
# geom_boxplot(outlier.size = 0.5, alpha = 0.7, width = 0.6) +

P <- ggplot(df, aes(x = group, y = score, color = group)) +
  geom_violin(fill = "grey90", color = NA, trim = FALSE, width = 1) +
  geom_boxplot(size = 0.8, outlier.shape = NA,fill=NA) +
  geom_jitter(width = 0.15, size = 0.7, alpha = 0.7) +
  geom_hline(yintercept = median_score, linetype = "dashed", color = "black", size = 1) +  # median line
  scale_color_manual(values = c("#3C5488","#E64B35")) +
  stat_compare_means(
    comparisons = my_comparisons,
    method = "wilcox.test",
    size = 4.5,
    label = "p.signif",
    label.y = 2.5
  ) +
  theme_classic() +
  labs(
    title = "",
    x = "Immune cluster",
    y = "TP53 score"
  ) +
  theme(
    axis.text.x = element_text(size = 12),       
    axis.text.y = element_text(size = 12),       
    axis.title.x = element_text(size = 14),      
    axis.title.y = element_text(size = 14),      
    plot.title = element_text(size = 16,face = "bold"),
    legend.position = "none"
  ) +
  coord_cartesian(ylim = c(NA, 3))

graphics.off()
pdf(file = paste0(FigDir,"TP53_Score_Cluster.pdf"), width = 2.8, height = 4) 
P
dev.off()


# par(mfrow=c(2,2))

#-----------------
se = which(row.names(info)%in%imm.hot)
data = info[se,]

myList = list(NULL)
myList[[1]] = data$score[data$Response==1]
myList[[2]] = data$score[data$Response==0]
boxplot(myList)
wilcox.test(myList[[1]], myList[[2]])


df_hot <- data.frame(
  score = c(myList[[1]],myList[[2]]),
  group = factor(c(rep("Responding", length(myList[[1]])),
                   rep("Non-Responding", length(myList[[2]])))))

my_comparisons <- list(c("Responding", "Non-Responding"))

# Plot with paired significance
P_Hot <- ggplot(df_hot, aes(x = group, y = score, color = group)) +
  geom_boxplot(size=0.6, outlier.shape = NA) +
  scale_color_manual(values = c("#E64B35","#3C5488")) +
  stat_compare_means(
    comparisons = my_comparisons,
    method = "wilcox.test",
    size = 4.5,
    label = "p.signif"
  ) +
  theme_classic() +
  labs(
    title = "IMbrave150\n(Hot Immune)",
    x = NULL,
    y = "TP53 score"
  ) +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text.y = element_text(size = 12),
    axis.title.y = element_text(size = 14),
    plot.title = element_text(size = 14),
    legend.position = "none"
  ) +
  coord_cartesian(ylim = c(NA, 65))

graphics.off()

pdf(file = paste0(FigDir,"TP53_Score_Response.pdf"), width = 2, height = 3.5) 
P_Hot
dev.off()




myList = list(NULL)
myList[[1]] = data$score[data$Confirmed.Response_IRF=="CR"]
myList[[2]] = data$score[data$Confirmed.Response_IRF=="PR"]
myList[[3]] = data$score[data$Confirmed.Response_IRF=="SD"]
myList[[4]] = data$score[data$Confirmed.Response_IRF=="PD"]
boxplot(myList)

wilcox.test(myList[[3]], myList[[4]])		# 
wilcox.test(myList[[3]], myList[[2]])		#  
wilcox.test(myList[[4]], myList[[2]])		#  


df<-data.frame(score=c(myList[[1]],myList[[2]],myList[[3]],myList[[4]]),
               response=c(rep("CR",length(myList[[1]])),rep("PR",length(myList[[2]])),rep("SD",length(myList[[3]])),rep("PD",length(myList[[4]]))))
df$response <- factor(df$response, levels = c("CR", "PR", "SD", "PD"))


mycol <-c(
  "CR" = "steelblue4",
  "PR" = "steelblue2", 
  "SD" = "grey60", 
  "PD" = "indianred4")


my_comparisons <- list(
  c("CR", "SD"),
  c("PR", "SD"),
  c("PD", "SD"),
  c("PD", "PR"))

mygg1 <- ggplot(df, aes(x=response, y=score, color= response)) + 
  geom_boxplot(size=1,outlier.size = 0.3,linewidth=1) +
  # geom_jitter(width = 0.15, size = 0.3, alpha = 0.7) +
  # scale_fill_manual(values="white") +
  scale_color_manual(values=mycol) +
  
  # scale_fill_manual(values=mycol) +
  labs(title='Ate+Bev\n(Hot Immune)', x= 'Response', y = 'TP53 score') +
  theme_classic() + 
  theme(legend.position="none") +
  theme(
    axis.text.x = element_text(size = 12),       
    axis.text.y = element_text(size = 12),       
    axis.title.x = element_text(size = 14),      
    axis.title.y = element_text(size = 14),      
    plot.title = element_text(size = 16, face = "bold")
    # ,
    # plot.margin = margin(10, 10, 10, 10)  # Add margins
  ) +
  stat_compare_means(
    comparisons = my_comparisons,
    method = "wilcox.test",      
    label = "p.signif",
    size = 4.5   # Increase size of p-values
  )+coord_cartesian(ylim = c(NA, 95))


fig3a<-mygg1

graphics.off()
pdf(paste0(FigDir,"3A_hot.pdf"), width = 2.6, height = 3.5) 
mygg1
dev.off()


#-----------------
se = which(row.names(info)%in%imm.col)
data = info[se,]

myList = list(NULL)
myList[[1]] = data$score[data$Response==1]
myList[[2]] = data$score[data$Response==0]

# myList = list(NULL)
# myList[[1]] = data$score[data$Confirmed.Response_IRF%in%c("CR","PR")]
# myList[[2]] = data$score[data$Confirmed.Response_IRF%in%c("PD")]


boxplot(myList)
wilcox.test(myList[[1]], myList[[2]])

df_cold <- data.frame(
  score = c(myList[[1]],myList[[2]]),
  group = factor(c(rep("Responding", length(myList[[1]])),
                   rep("Non-Responding", length(myList[[2]])))))

my_comparisons <- list(c("Responding", "Non-Responding"))


df_cold <- data.frame(
  score = c(myList[[1]], myList[[2]]),
  group = c(rep("Responding", length(myList[[1]])),
            rep("Non-Responding", length(myList[[2]]))),
  Immune = "Cold"
)
se = which(row.names(info) %in% imm.hot)
data_hot = info[se,]

myList_hot = list(NULL)
myList_hot[[1]] = data_hot$score[data_hot$Response==1]
myList_hot[[2]] = data_hot$score[data_hot$Response==0]

# myList_hot[[1]] = data_hot$score[data_hot$Confirmed.Response_IRF%in%c("CR","PR")]
# myList_hot[[1]] = data_hot$score[data_hot$Confirmed.Response_IRF%in%c("PD")]

df_hot <- data.frame(
  score = c(myList_hot[[1]], myList_hot[[2]]),
  group = c(rep("Responding", length(myList_hot[[1]])),
            rep("Non-Responding", length(myList_hot[[2]]))),
  Immune = "Hot"
)

df_all <- rbind(df_hot, df_cold)

df_all$group <- factor(df_all$group,
                       levels = c("Responding","Non-Responding"))

df_all$Immune <- factor(df_all$Immune,
                        levels = c("Hot","Cold"))

my_comparisons <- list(c("Responding","Non-Responding"))

P_HotCold <- ggplot(df_all,
                    aes(x = group, y = score, color = group)) +
  geom_boxplot(size = 0.6, outlier.shape = NA) +
  
  facet_wrap(~Immune,
             nrow = 2,
             scales = "fixed") +   # shared y-scale
  
  scale_color_manual(values = c("royalblue3","lightsalmon3")) +
  
  stat_compare_means(
    comparisons = my_comparisons,
    method = "wilcox.test",
    label = "p.signif",
    size = 5
  ) +
  
  labs(
    x = "",
    y = "TP53 score"
  ) +
  
  theme_classic() +
  
  theme(
    axis.text.x = element_text(size = 12, angle = 45, hjust = 1),
    axis.text.y = element_text(size = 12),
    axis.title.y = element_text(size = 14),
    strip.text = element_text(size = 14, face = "bold"),
    legend.position = "none"
  ) +
  
  coord_cartesian(ylim = c(NA,70))

P_HotCold

graphics.off()
pdf(file = paste0(FigDir,"TP53_Score_Response_cold_hot1.pdf"), width = 2, height = 4.5) 
P_HotCold
dev.off()

P_Cold <- ggplot(df_cold, aes(x = group, y = score, color = group)) +
  geom_boxplot(size=0.6, outlier.shape = NA) +
  scale_color_manual(values = c("#E64B35","#3C5488")) +
  stat_compare_means(
    comparisons = my_comparisons,
    method = "wilcox.test",
    size = 4.5,
    label = "p.signif"
  ) +
  theme_classic() +
  labs(
    title = "(Cold Immune)",
    x = "",
    y = "TP53 score"
  ) +
  theme(
    axis.text.x = element_text(size = 12, angle = 45, hjust = 1),
    axis.text.y = element_text(size = 12),
    axis.title.x = element_text(size = 14),
    axis.title.y = element_text(size = 14),
    plot.title = element_text(size = 14),
    legend.position = "none"
  ) +
  coord_cartesian(ylim = c(NA, 70))

graphics.off()
pdf(file = paste0(FigDir,"TP53_Score_Response_cold.pdf"), width = 2, height = 3.5) 
P_Cold
dev.off()




myList = list(NULL)
myList[[1]] = data$score[data$Confirmed.Response_IRF=="CR"]
myList[[2]] = data$score[data$Confirmed.Response_IRF=="PR"]
myList[[3]] = data$score[data$Confirmed.Response_IRF=="SD"]
myList[[4]] = data$score[data$Confirmed.Response_IRF=="PD"]
boxplot(myList)

wilcox.test(myList[[3]], myList[[4]])		# 
wilcox.test(myList[[3]], myList[[2]])		#  
wilcox.test(myList[[4]], myList[[2]])		#  


df<-data.frame(score=c(myList[[1]],myList[[2]],myList[[3]],myList[[4]]),
               response=c(rep("CR",length(myList[[1]])),rep("PR",length(myList[[2]])),rep("SD",length(myList[[3]])),rep("PD",length(myList[[4]]))))
df$response <- factor(df$response, levels = c("CR", "PR", "SD", "PD"))


mycol <-c(
  "CR" = "steelblue4",
  "PR" = "steelblue2", 
  "SD" = "grey60", 
  "PD" = "indianred4")

my_comparisons <- list(
  c("CR", "SD"),
  c("PR", "SD"),
  c("PD", "SD"),
  c("PD", "PR"))

mygg1 <- ggplot(df, aes(x=response, y=score, color= response)) + 
  geom_boxplot(size=1,outlier.size = 0.3,linewidth=1) +
  # geom_jitter(width = 0.15, size = 0.3, alpha = 0.7) +
  # scale_fill_manual(values="white") +
  scale_color_manual(values=mycol) +
  
  # scale_fill_manual(values=mycol) +
  labs(title='Ate+Bev\n(Cold Immune)', x= 'Response', y = 'TP53 score') +
  theme_classic() + 
  theme(legend.position="none") +
  theme(
    axis.text.x = element_text(size = 12),       
    axis.text.y = element_text(size = 12),       
    axis.title.x = element_text(size = 14),      
    axis.title.y = element_text(size = 14),      
    plot.title = element_text(size = 16, face = "bold")
    # ,
    # plot.margin = margin(10, 10, 10, 10)  # Add margins
  ) +
  stat_compare_means(
    comparisons = my_comparisons,
    method = "wilcox.test",      
    label = "p.signif",
    size = 4.5   # Increase size of p-values
  )+coord_cartesian(ylim = c(NA, 95))


fig3a<-mygg1

graphics.off()
pdf(paste0(FigDir,"3A_cold.pdf"), width = 2.6, height = 3.5) 
mygg1
dev.off()




info$Immune
info$score

mean_score <- mean(info$score, na.rm = TRUE)
info$Score_group <- ifelse(info$score >= mean_score, "High", "Low")
info$Group4 <- paste(info$Immune, info$Score_group, sep = "_")

info$Group4 <- factor(
  info$Group4,
  levels = c(
    "Hot_High",
    "Hot_Low",
    "Cold_High",
    "Cold_Low"
  )
)


dat <- data.frame(
  combined  = info$Group4,
  P53.group = info$Score_group,
  Resp = info$Confirmed.Response_IRF,
  OS.time = info$OS.in.days/30.0,
  OS.event = ifelse(info$OS.censoring..1.cens.0.evt. == 0, 1, 0)  # 1 = death, 0 = censored
)


# filter response categories
dat = dat[dat$Resp %in% c("CR", "PR", "SD", "PD"), ]




# ============================================================
# KM plotting function for one category
# ============================================================

km_plot <- function(cat, dat,Title) {
  
  sub = dat[dat$Resp %in% cat, ]
  if(nrow(sub) < 5) {
    return(ggplot() + ggtitle(paste(cat, "(N too small)")))
  }
  
  fit = survfit(Surv(OS.time, OS.event) ~ combined, data=sub)
  
  ggsurvplot(
    fit,
    data = sub,
    pval = TRUE,
    title = Title,
    legend.title = "TP53 score",
    # legend.labs = c("High", "Low"),
    # palette = c("coral1","steelblue"),
    risk.table = FALSE,     # no risk table in subplots
    xlab = "OS(months)",
    ylab = "Survival probability",
    censor.shape = "|",
    censor.size = 2,
    ggtheme = theme_classic(base_size = 12)+
      theme(
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        axis.text.x  = element_text(size = 12),
        axis.text.y  = element_text(size = 12),
        plot.title   = element_text(size = 14, hjust = 0.5)  # center title
      ))$plot  # extract only the KM curve
}

# ============================================================
# Generate each subplot
# ============================================================

p_SD = km_plot(c("CR","PR","SD","PD"), dat,"SD response category")
km_plot("PD", dat,"SD response category")


graphics.off()
pdf(paste0(FigDir,"3D2.pdf"),height = 4,width = 3)
print(p_SD)
dev.off()





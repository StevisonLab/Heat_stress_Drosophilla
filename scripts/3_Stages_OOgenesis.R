## Load data
Ova_A <- read_csv("rawdata/Ovarioles.csv", col_types = cols(Age = col_number(), Number_Oocyte = col_number(), TUNEL_Cell = col_number(), Ovary = col_number(), Ovariole = col_number(), Replicate = col_number()))
str(Ova_A)

Ova_A$Replicate_Cor=as.factor(Ova_A$Replicate_Cor)
Ova_A$Species=as.factor(Ova_A$Species)
Ova_A$Treatment=as.factor(Ova_A$Treatment)
Ova_A$Stage=as.factor(Ova_A$Stage)
Ova_A$Age=as.factor(Ova_A$Age)
Ova_A$Ovary=as.factor(Ova_A$Ovary)
Ova_A$Ovariole=as.factor(Ova_A$Ovariole)

## Create unique replica ID
Ova_A$ID = paste(Ova_A$Replicate_Cor, Ova_A$Ovary, Ova_A$Ovariole, sep="_")
Ova_A$Ovary_cor = paste(Ova_A$Replicate_Cor, Ova_A$Ovary, sep="_")
Ova_A$Ovariole_cor = paste(Ova_A$Replicate_Cor, Ova_A$Ovary, Ova_A$Ovariole, sep="_")
Ova_A$Number_Oocyte=as.numeric(Ova_A$Number_Oocyte)

## Make TUNEL binary
Ova_A$TUNEL_Ovariole=ifelse(Ova_A$TUNEL_Cell==1, "Positive", "Negative") 
Ova_A$TUNEL_Ovariole=as.factor(Ova_A$TUNEL_Ovariole)


##Statistical analysis
#Oocytes
fit_Oocyte <- glmer.nb(Number_Oocyte ~ (1|Replicate_Cor)  + (1|Ovariole_cor) + (1|Replicate_Cor:Ovariole_cor) + Stage * Timepoint * Treatment * Species, data = Ova_A)
lsm <- pairs(emmeans(fit_Oocyte, ~ Stage * Timepoint * Treatment * Species))
summary(fit_Oocyte)$coefficients
summary(lsm, type = "response")
emmeans(res.aov7, ~ Species * Treatment)
anova_fitOoocyte <- Anova(fit_Oocyte)
anova_fitOoocyte
summary(fit_Oocyte)
cld9 <- cld(lsm, Letters = letters)
fit_Oocyte

##Save model tables
write_csv(as.data.frame(summary(fit_Oocyte)$coef), "Stats/Oocyte_glmer.csv", col_names = T)
write_csv(as.data.frame(summary(lsm, type = "response")), "Stats/Oocyte_emmeans.csv")

#TUNEL positive
res.aov7 <- lmer(TUNEL_Cell ~ (1|Replicate_Cor:Ovariole_cor) + Timepoint + Treatment + Species + Species:Treatment + Stage + Treatment:Timepoint:Stage:Species + Species:Treatment:Timepoint, data = Ova_A)
summary(res.aov7)
lsmT <- emmeans(res.aov7, ~ Timepoint * Treatment * Species)
summary(res.aov7)$coefficients
summary(lsmT, type = "response")
anova_TUNEL <- Anova(res.aov7)
cldTS <- cld(lsmT, Letters = letters)
summary(res.aov7)

##Save model tables
write.csv(summary(res.aov7)$coefficients, "Stats/TUNEL_Lmer.csv")
write_csv(summary(lsm), "Stats/TUNEL_emmeans.csv")

##Sub-setting for plotting
## Group data for visualization
New_Ovaries= summaryBy(Number_Oocyte+TUNEL_Cell~Species+Timepoint+Treatment+Stage+Replicate+Ovariole,data=Ova_A,FUN=sum,na.rm=T)

## Subset species and Stages
DMEL<-subset(Ova_A, Ova_A$Species=="DMEL", na.rm=TRUE, select = c(Species,Age, Stage, Number_Oocyte, Treatment, Ovariole, Ovary, TUNEL_Cell, Replicate_Cor, Timepoint, TUNEL_Ovariole))
DPSE<-subset(Ova_A, Ova_A$Species=="DPSE", na.rm=TRUE, select = c(Species,Age, Stage, Number_Oocyte, Treatment, Ovariole, Ovary, TUNEL_Cell, Replicate_Cor, Timepoint, TUNEL_Ovariole))
Stage11 <- subset(Ova_A, Ova_A$Stage=="Stage 11", na.rm=TRUE, select = c(Species,Age, Stage, Number_Oocyte, Treatment, Ovariole, Ovary, TUNEL_Cell, Replicate_Cor, Timepoint, TUNEL_Ovariole))
Stage1_7 <- subset(Ova_A, Ova_A$Stage=="Stages 1-7", na.rm=TRUE, select = c(Species,Age, Stage, Number_Oocyte, Treatment, Ovariole, Ovary, TUNEL_Cell, Replicate_Cor, Timepoint, TUNEL_Ovariole))
DPSEStage1_7 <- subset(Stage1_7, Stage1_7$Species=="DPSE", na.rm=TRUE, select = c(Species,Age, Stage, Number_Oocyte, Treatment, Ovariole, Ovary, TUNEL_Cell, Replicate_Cor, Timepoint, TUNEL_Ovariole))
DMELStage1_7 <- subset(Stage1_7, Stage1_7$Species=="DMEL", na.rm=TRUE, select = c(Species,Age, Stage, Number_Oocyte, Treatment, Ovariole, Ovary, TUNEL_Cell, Replicate_Cor, Timepoint, TUNEL_Ovariole))
Ova_A$TUNEL_Cell=as.numeric((Ova_A$TUNEL_Cell))

##Subset data from models to include in plots
cld9$group=ifelse(cld9$.group==" a       ","a", ifelse(cld9$.group=="       gh","gh", ifelse(cld9$.group=="  bcdef  ","bcdef", ifelse(cld9$.group=="        h","h", ifelse(cld9$.group=="  bcde   ","bcde", ifelse(cld9$.group=="    de   ","de", ifelse(cld9$.group=="  b      ","b", ifelse(cld9$.group=="  bc     ","bc", ifelse(cld9$.group=="   cd    ","cd", ifelse(cld9$.group=="     e   ","e", ifelse(cld9$.group=="      fg ","fg", ifelse(cld9$.group=="   cde   ","cde", ifelse(cld9$.group=="  bcd    ","bcd", cld9$.group)))))))))))))
cld1_7 <- subset(cld9, cld9$Stage=="Stages 1-7", na.rm=TRUE, select = c(Species, Treatment, Stage, Timepoint, .group, emmean, SE))
cldDMEL9 <- subset(cld9, cld9$Species=="DMEL", na.rm=TRUE, select = c(Species, Treatment, Stage, Timepoint, .group, emmean, SE))
cldDPSE9 <- subset(cld9, cld9$Species=="DPSE", na.rm=TRUE, select = c(Species, Treatment, Stage, Timepoint, .group, emmean, SE))

##Summary for data merge with significant letters
OOS=summaryBy(Number_Oocyte~Treatment+Species+Timepoint+Stage,data=Ova_A,na.rm=T, fun.max=max(Ova_A$Number_Oocyte), stringsAsFactors = FALSE)
TS=summaryBy(TUNEL_Cell~Treatment+Species+Timepoint+Stage+Ovariole_cor,data=Ova_A,na.rm=T,fun=max(Ova_A$TUNEL_Cell), stringsAsFactors = FALSE)
TS=summaryBy(TUNEL_Cell.1~Treatment+Species+Timepoint+Stage,data=TS,na.rm=T,FUN=sum, stringsAsFactors = FALSE)
OOS.summarized=merge(OOS,cld9, stringsAsFactors = FALSE, na.rm=TRUE, check.names = T)
TS.summarized=merge(TS,cldTS, stringsAsFactors = FALSE, na.rm=TRUE, check.names = F)

OOS.summarized$group=ifelse(OOS.summarized$.group==" a       ","a", ifelse(OOS.summarized$.group=="       gh","gh", ifelse(OOS.summarized$.group=="  bcdef  ","bcdef", ifelse(OOS.summarized$.group=="        h","h", ifelse(OOS.summarized$.group=="  bcde   ","bcde", ifelse(OOS.summarized$.group=="    de   ","de", ifelse(OOS.summarized$.group=="  b      ","b", ifelse(OOS.summarized$.group=="  bc     ","bc", ifelse(OOS.summarized$.group=="   cd    ","cd", ifelse(OOS.summarized$.group=="     e   ","e", ifelse(OOS.summarized$.group=="      fg ","fg", ifelse(OOS.summarized$.group=="   cde   ","cde", ifelse(OOS.summarized$.group=="  bcd    ","bcd", OOS.summarized$.group)))))))))))))
OOSdmel = subset(OOS.summarized, OOS.summarized$Species=="DMEL" & OOS.summarized$Stage=="Stages 1-7", na.rm=TRUE, select = c(Species, Timepoint, Stage, Treatment, Number_Oocyte.mean, group, emmean, SE))
OOSdpse = subset(OOS.summarized, OOS.summarized$Species=="DPSE" & OOS.summarized$Stage=="Stages 1-7", na.rm=TRUE, select = c(Species, Timepoint, Stage, Treatment, Number_Oocyte.mean, group, emmean, SE))
OOS1_7 = subset(OOS.summarized, OOS.summarized$Stage=="Stages 1-7", na.rm=TRUE, select = c(Species, Timepoint, Stage, Treatment, Number_Oocyte.mean, group, emmean, SE))
OOS11 = subset(OOS.summarized, OOS.summarized$Stage=="Stage 11", na.rm=TRUE, select = c(Species, Timepoint, Stage, Treatment, Number_Oocyte.mean, group, emmean, SE))

TSEr = subset(TS.summarized, TS.summarized$Timepoint=="Early", na.rm=TRUE, select = c(Species, Timepoint, Stage, Treatment,TUNEL_Cell.1.sum, .group, emmean, SE))
TSEr$group=ifelse(TSEr$.group==" ab ","ab", ifelse(TSEr$.group==" a  ","a", ifelse(TSEr$.group=="  b ","b", ifelse(TSEr$.group=="   c","c", TSEr$.group))))
OOS1_7$group2=ifelse(OOS1_7$group=="gh","ab", ifelse(OOS1_7$group=="h","a", ifelse(OOS1_7$group=="fg","bc",  OOS1_7$group)))
OOS11$SE2=ifelse(OOS11$Number_Oocyte.mean==0, 0, OOS11$SE)
Oocytes_stats

## Figure 5 Barplots
StagesM<-ggplot(OOS1_7, aes(x=Timepoint, y=Number_Oocyte.mean, fill=Treatment)) +
  ylab("Number of Oocytes per ovariole") +
  theme(axis.text.y = element_text(size = 14)) +
  rremove("xlab") +
  rremove("ylab") +
  theme(axis.text.x=element_blank()) + 
  theme(legend.position = "none") +
  scale_fill_manual(breaks = c("Control", "HighTemp"), values=c("blue","red")) +
  stat_summary(fun= "identity", geom="bar", position = position_dodge(0.9)) + 
  theme(panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  facet_wrap(~Species) + theme(strip.text.x = element_blank()) + geom_text(aes(y = Number_Oocyte.mean+0.4, x= Timepoint, label = str_trim(group2)), size = 4, position = position_dodge(0.8),
                                                                         hjust = 0.7) +
  geom_errorbar(aes(
                  ymin = Number_Oocyte.mean - SE,
                  ymax = Number_Oocyte.mean + SE,
                  x = Timepoint,
                  
                ),
                width = 0.1, position = position_dodge(0.8)) 
  
StagesM

StagesP<-ggplot(OOS11, aes(x=Timepoint, y=Number_Oocyte.mean, fill=Treatment)) +
  ylab("Number of Oocytes per ovariole") +
  theme(axis.text.y = element_text(size = 14)) +
  theme(axis.text.x = element_text(size = 14)) +
  rremove("xlab") +
  rremove("ylab") +
  theme(legend.position = "none") +
  scale_fill_manual(breaks = c("Control", "HighTemp"), values=c("blue","red")) +
  stat_summary(fun= "identity", geom="bar", position = position_dodge(0.9)) + 
  theme(panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  facet_wrap(~Species) +theme(strip.text.x = element_blank())+ geom_text(aes(y = Number_Oocyte.mean+0.06, x= Timepoint, label = str_trim(group)), size = 4, position = position_dodge(0.85),
                                                                         hjust = 0.7) +
  geom_errorbar(aes(
    ymin = Number_Oocyte.mean-SE2/10,
    ymax = Number_Oocyte.mean+SE2/10,
    x = Timepoint,
    
  ),
  width = 0.1, position = position_dodge(0.8)) 
StagesP

pdf("Figures/Fig5.pdf")
Combined6=ggarrange(ncol=1, StagesM, StagesP,common.legend = TRUE, font.label = list(size = 1, color = "black", face = "bold"), align = "hv", legend = "none")
Combined6=annotate_figure(Combined6, left = textGrob("Number of Oocytes per ovariole", rot = 90, vjust = 1, gp = gpar(cex = 1.3)))
Combined6=annotate_figure(Combined6, top = text_grob("D. melanogaster                         D. pseudoobscura", color = "black", face="bold.italic", family="sans", size = 15))
Combined6
dev.off()

##Figure 6
pdf("Figures/Fig6.pdf")
TUNEL <- ggplot(data=TSEr, aes(x=Treatment, y=TUNEL_Cell.1.sum, fill=Treatment, alpha=Stage)) +
  #scale_fill_manual(breaks = c("Stages 1-7", "Stages 8-10", "Stage 11", "Stage 12-14"), values=c("#eff3ff","#bdd7e7","#6baed6", "#2171b5")) +
  scale_fill_manual(breaks = c("Control", "HighTemp"), values=c("#2171b5","red")) +
  ylab("Oocytes TUNEL-positive per ovariole") +
  theme(axis.text.y = element_text(size = 14)) +
  theme(axis.text.x = element_text(size = 14)) +
  scale_alpha_manual(values=c(0.1, 0.25, 0.50, 0.95)) +
  theme(panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  stat_summary(fun= "identity", geom="bar", position = "stack") + facet_wrap(~Species) + theme(strip.text.x = element_blank()) 
  #geom_text(aes(y = TUNEL_Cell.1.sum, x= Treatment, label = str_trim(group)), size = 4, position = position_dodge(0.1),  hjust = 1, check_overlap = T) 
TUNEL                                                                                                                                                                                                                                                                                   
dev.off()
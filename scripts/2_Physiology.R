## Load data
CT<-read_csv("rawdata/CTmax.csv", col_types =cols(CT_max= col_number(), Comments = col_skip(), ...12 = col_skip()))
resp <- read_csv("rawdata/RQ.csv", col_types = cols(`Pre-Weight` = col_number(), `Post-Weight` = col_number(), Area_O2 = col_number(), 
                                  Area_CO2 = col_number(), Comment = col_skip()))
resp<-na.omit(resp)

## Calculate mass lost from Respirometry
resp$Mass_Lost = resp$`Pre-Weight` - resp$`Post-Weight`

## Convert raw data to ml
resp$Area_O2cor = (resp$Area_O2)/100
resp$Area_CO2cor = (resp$Area_CO2)/100

## Calculate O2 and CO2 per hour per grams and RQ
resp$Time = resp$Hours + (resp$Minutes/60)
resp$Vol_Correction = resp$`Chm Vol`+ 0.048 + 0.03 - (10.1*resp$`Pre-Weight`)
resp$O2pHr = (resp$Area_O2cor*((resp$Vol_Correction/resp$`Inj Vol`)/resp$Time))
resp$CO2pHr = (resp$Area_CO2cor*((resp$Vol_Correction/resp$`Inj Vol`)/resp$Time))
resp$O2pHrpGr = (resp$O2pHr/resp$`Pre-Weight`)
resp$CO2pHrpGr = (resp$CO2pHr/resp$`Pre-Weight`)
resp$RQ = resp$CO2pHrpGr/resp$O2pHrpGr

## Convert grams to mg and ml to ul
resp$Weight = resp$`Pre-Weight` * 1000
resp$CO2pHr2 = resp$CO2pHr *1000
resp$O2pHr2 = resp$O2pHr *1000

## Functions to add N to the boxplot
#CTmax Dmel
stat_box_data <- function(y, upper_limit = 33) {
  return( 
    data.frame(
      y =  upper_limit,
      label = paste('n =', length(y))
    )
  )
}
#CTmax Dpse
stat_box_data1 <- function(y, upper_limit = 42) {
  return( 
    data.frame(
      y =  upper_limit,
      label = paste('n =', length(y))
    )
  )
}

#RQ
stat_box_data2 <- function(y, upper_limit = 1.62) {
  return( 
    data.frame(
      y = 0.8 * upper_limit,
      label = paste('n =', length(y))
    )
  )
}
#O2
stat_box_data3 <- function(y, upper_limit = 5) {
  return( 
    data.frame(
      y = 0.8 * upper_limit,
      label = paste('n =', length(y))
    )
  )
}
#CO2
stat_box_data4 <- function(y, upper_limit = 4) {
  return( 
    data.frame(
      y = 0.8 * upper_limit,
      label = paste('n =', length(y))
    )
  )
}
#Weight
stat_box_data5 <- function(y, upper_limit = max(resp$Weight) * 0.1) {
  return( 
    data.frame(
      y = 0.8 * upper_limit,
      label = paste('n =', length(y))
    )
  )
}
CT$Sex=CT$`Sex `
## Subset by species 
CTDMEL<-subset(CT, CT$Species=="DMEL", na.rm=TRUE, select = c(Species, Sex, Treatment, CT_max))
CTDPSE<-subset(CT, CT$Species=="DPSE", na.rm=TRUE, select = c(Species, Sex, Treatment, CT_max))
respDPSE<-subset(resp, resp$Species=="DPSE", na.rm=TRUE, select = c(Species, Sex, Vial, Treatment, `Pre-Weight`, O2pHr2, CO2pHr2, RQ, Weight, Area_O2cor, Area_CO2cor))     
respDMEL<-subset(resp, resp$Species=="DMEL", na.rm=TRUE, select = c(Species, Sex, Vial, Treatment, `Pre-Weight`, O2pHr2, CO2pHr2, RQ, Weight, Area_O2cor, Area_CO2cor))

#Summary of the data
fun <- function(x){
  c(m=mean(x), v=var(x), n=length(x), min=min(x), max=max(x))
}
summary_by(cbind(O2pHr2,CO2pHr2, Weight, RQ) ~ Treatment, data=respDMEL, FUN=fun)
summary_by(cbind(O2pHr2,CO2pHr2, Weight, RQ) ~ Treatment, data=respDPSE, FUN=fun)
summary_by(CT_max ~ Treatment, data=CTDPSE, FUN=fun)

##Statistical analysis per measurement
##RQ, CO2 produced and O2 consumed 
res.aov4 <- lmer(RQ  ~ (1|Vial) + Treatment * Species * Sex, data = resp)
as.data.frame(summary(res.aov4)$coefficients)
write_csv2(as.data.frame(summary(res.aov4)$coefficients), "Stats/RQ_lmer.csv", col_names = T)
anova_res.aov4 <- Anova(res.aov4)
write_csv2(as.data.frame(anova_res.aov4), "Stats/RQ_aov.csv", col_names = T)
RQS=summaryBy(RQ~Treatment+Species+Sex,data=resp,FUN = mean,na.rm=T, stringsAsFactors = FALSE)
lsmRQ<- emmeans(res.aov4, ~ Species + Treatment)
summary(lsmRQ)
cldRQ <-cld(lsmRQ, Letters = letters, alpha = 0.05)


res.aov5 <- lmer(CO2pHr2  ~ (1|Vial) +Treatment * Species * Sex, data = resp)
as.data.frame(summary(res.aov5)$coefficients)
write_csv2(as.data.frame(summary(res.aov5)$coefficients), "Stats/CO2_lmer.csv", col_names = T)
anova_res.aov5 <- Anova(res.aov5)
write_csv2(as.data.frame(anova_res.aov5), "Stats/CO2_aov.csv", col_names = T)
CO2S=summaryBy(CO2pHr2~Treatment+Species+Sex,data=resp,FUN = mean,na.rm=T, stringsAsFactors = FALSE)

res.aov6 <- lmer(O2pHr2  ~ (1|Vial) +Treatment * Species * Sex, data = resp)
as.data.frame(summary(res.aov6)$coefficients)
write_csv2(as.data.frame(summary(res.aov6)$coefficients), "Stats/O2_lmer.csv", col_names = T)
anova_res.aov6 <- Anova(res.aov6)
write_csv2(as.data.frame(anova_res.aov6), "Stats/O2_aov.csv", col_names = T)
O2S=summaryBy(O2pHr2~Treatment+Species+Sex,data=resp,FUN = mean,na.rm=T, stringsAsFactors = FALSE)

res.aov7 <- lmer(Weight  ~ (1|Vial) +Treatment * Species * Sex, data = resp)
as.data.frame(summary(res.aov7)$coefficients)
write_csv2(as.data.frame(summary(res.aov7)$coefficients), "Stats/Weight_lmer.csv", col_names = T)
anova_res.aov7 <- Anova(res.aov7)
write_csv2(as.data.frame(anova_res.aov7), "Stats/Weight_aov.csv", col_names = T)
lsm7 <- emmeans(res.aov7, ~ Species * Treatment)
write_csv2(summary(lsm7), "Stats/emmeans_Weight.csv")

##Extract significant letters
cld7 <-cld(lsm7, Letters = letters, alpha = 0.05)

##Subset significant letters for plots
cldDMEL7 <- subset(cld7, cld7$Species=="DMEL", na.rm=TRUE, select = c(Species, Sex, Treatment, .group, emmean, SE))
cldDPSE7 <- subset(cld7, cld7$Species=="DPSE", na.rm=TRUE, select = c(Species, Sex, Treatment, .group, emmean, SE))

WS=summaryBy(Weight~Treatment+resp$Species+Sex,data=resp,FUN = mean,na.rm=T, stringsAsFactors = FALSE)
WS.summarized=merge(WS,cld7, stringsAsFactors = FALSE, na.rm=TRUE, check.names = F)
WS.summarized$group=ifelse(WS.summarized$.group==" abc","abc", ifelse(WS.summarized$.group==" ab ","ab", ifelse(WS.summarized$.group=="  bc","bc", ifelse(WS.summarized$.group==" a  ","a", ifelse(WS.summarized$.group=="   c","c",  WS.summarized$.group)))))
WSdmel = subset(WS.summarized, WS.summarized$Species=="DMEL", na.rm=TRUE, select = c(Species, Sex, Treatment, Weight.mean, group))
WSdpse = subset(WS.summarized, WS.summarized$Species=="DPSE", na.rm=TRUE, select = c(Species, Sex, Treatment, Weight.mean, group))

##CT max
res.aov8 <- lmer(CT_max  ~ (1|Vial) +Treatment * Species * Sex, data = CT)
as.data.frame(summary(res.aov8)$coefficients)
write_csv2(as.data.frame(summary(res.aov8)$coefficients), "Stats/CTmax_lmer.csv", col_names = T)
anova_res.aov8 <- Anova(res.aov8)
write_csv2(as.data.frame(anova_res.aov8), "Stats/CTmax_aov.csv", col_names = T)
lsm8<- emmeans(res.aov8, ~ Species+Treatment+Sex)
write_csv2(summary(lsm8), "Stats/emmeans_CTmax.csv")
sumlsm8 <- summary(lsm8)
cld8 <-cld(lsm8, Letters = letters, alpha = 0.05)

##Subset significant letters for plots
CTS=summaryBy(CT_max~Treatment+Species+Sex,data=CT,FUN = mean,na.rm=T, stringsAsFactors = FALSE)
CT.summarized=merge(CTS,cld8, stringsAsFactors = FALSE, na.rm=TRUE, check.names = F)
CT.summarized$group=ifelse(CT.summarized$.group=="        h ","h", ifelse(CT.summarized$.group=="       gh ","gh", ifelse(CT.summarized$.group=="   c      ","c", ifelse(CT.summarized$.group=="      fg  ","fg", ifelse(CT.summarized$.group=="      fgh ","fgh", ifelse(CT.summarized$.group==" a        ","a", ifelse(CT.summarized$.group=="     ef   ","ef", ifelse(CT.summarized$.group=="    d     ","d", ifelse(CT.summarized$.group=="         i","i", ifelse(CT.summarized$.group=="     e    ","e", ifelse(CT.summarized$.group==" ab       ","ab", ifelse(CT.summarized$.group=="  bc      ","bc", CT.summarized$.group))))))))))))
CTdmel = subset(CT.summarized, CT.summarized$Species=="DMEL", na.rm=TRUE, select = c(Species, Sex, Treatment, CT_max.mean, group))
CTdpse = subset(CT.summarized, CT.summarized$Species=="DPSE", na.rm=TRUE, select = c(Species, Sex, Treatment, CT_max.mean, group))

##Graphical abstract plot
GA = ggplot(sumlsm8, aes(x= Treatment, y= emmean, group = Species, color= Species))+
  rremove("xlab") + 
  theme(axis.title.y = element_blank()) +
  theme(axis.text.y = element_blank()) +
  theme(axis.text.x=element_blank()) +
  scale_color_manual(values=c('darkgreen', 'orange')) +
  geom_line(size=3) + facet_wrap(~Sex) 
GA

### Figure 4
## Panel A (CTmax)
# D. melanogaster
CT_DMEL=ggplot(CTDMEL, aes(x=Treatment, y=CT_max, fill=Treatment)) + 
  ylab(bquote(Critical~thermal~maximum~(CT[max]~"°C"))) +
  ylim(33,43) +
  rremove("xlab") +
  theme(axis.title.y = element_text(size = 11)) +
  theme(axis.text.y = element_text(size = 10)) +
  theme(axis.text.x=element_blank()) +
  theme(legend.position = "none") +
  stat_summary(fun.data = stat_box_data, geom = "text", fun = median,
               position = position_dodge(width = 0.75), size =4) +
  geom_text(data = CTdmel, aes(y = CT_max.mean, x = Treatment, label = group), vjust = -3, position = position_dodge(width = 1), size = 4) +
  geom_boxplot(position = position_dodge(width = 1)) +
  scale_fill_manual(breaks = c("C-C", "C-H", "H-C", "H-H"), values=c("#ffffb2","#fecc5c","#fd8d3c", "#e31a1c")) +
  geom_boxplot() + facet_wrap(~Sex) + theme( strip.background = element_rect(color="black", fill="white", size=1.5, linetype="solid"), strip.text.x = element_text(size = 10, color = "black", face = "bold"))
CT_DMEL


## D. pseudoobscura
CT_DPSE=ggplot(CTDPSE, aes(x=Treatment, y= CT_max, fill=Treatment)) + 
  ylab(bquote(Critical~thermal~maximum~(CT[max]~"°C"))) +
  ylim(33,43) +
  theme(axis.title.y = element_text(size = 10)) +
  theme(axis.text.y = element_text(size = 10)) +
  theme(axis.text.x=element_blank()) +
  rremove("ylab") +
  rremove("xlab") +
  theme(legend.position = "none") +
  theme(axis.text.x=element_blank()) +
  theme(axis.text.y=element_blank()) +
  stat_summary(fun.data = stat_box_data1, geom = "text", fun = median,
               position = position_dodge(width = 0.75), size=4) +
  scale_fill_manual(breaks = c("C-C", "C-H", "H-C", "H-H"), values=c("#ffffb2","#fecc5c","#fd8d3c", "#e31a1c")) +
  geom_text(data = CTdpse, aes(y = CT_max.mean, x = Treatment, label = group), vjust = -3, position = position_dodge(width = 1), size = 4) +
  geom_boxplot() + facet_wrap(~Sex) + theme( strip.background = element_rect(color="black", fill="white", size=1.5, linetype="solid"), strip.text.x = element_text(size = 10, color = "black", face = "bold"))
CT_DPSE


## Panel B (RQ)
# D. melanogaster
RQ_DMEL=ggplot(respDMEL, aes(x=Treatment, y=RQ, fill=Treatment)) + 
  ylim(0.6,1.3) +
  ylab("Respiratory Quotient (RQ)") +
  theme(axis.text.x=element_blank()) +
  theme(axis.title.y = element_text(size = 11)) +
  theme(axis.text.y = element_text(size = 10)) +
  rremove("xlab") +
  theme(legend.position = "none") +
  scale_fill_manual(breaks = c("C-C", "C-H", "H-C", "H-H"), values=c("#ffffb2","#fecc5c","#fd8d3c", "#e31a1c")) +
  stat_summary(fun.data = stat_box_data2, geom = "text", fun = median,
               position = position_dodge(width = 0.75), size=4) +
  geom_boxplot() + facet_wrap(~Sex) + theme( strip.background = element_rect(color="white", fill="white", size=1.5, linetype="blank"), strip.text.x = element_text(size = 14, color = "white", face = "bold"))
RQ_DMEL


# D. pseudoobscura
RQ_DPSE=ggplot(respDPSE, aes(x=Treatment, y=RQ, fill=Treatment)) + 
  ylim(0.6,1.3) +
  stat_summary(fun.data = stat_box_data2, geom = "text", fun = median,
               position = position_dodge(width = 0.75), size=4) +
  theme(plot.title = element_text(hjust = 0.5)) +
  rremove("xlab") +
  rremove("ylab") +
  theme(axis.text.x=element_blank()) +
  theme(axis.text.y=element_blank()) +
  theme(legend.position = "none") +
  scale_fill_manual(breaks = c("C-C", "C-H", "H-C", "H-H"), values=c("#ffffb2","#fecc5c","#fd8d3c", "#e31a1c")) +
  
  geom_boxplot() + facet_wrap(~Sex) + theme( strip.background = element_rect(color="white", fill="white", size=1.5, linetype="blank"), strip.text.x = element_text(size = 14, color = "white", face = "bold"))
RQ_DPSE

## Panel C (O2)
# D. melanogaster
O2DMEL=ggplot(respDMEL, aes(x=Treatment, y=O2pHr2, fill=Treatment)) + 
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(axis.text.x=element_blank()) +
  theme(axis.title.y = element_text(size = 11)) +
  theme(axis.text.y = element_text(size = 10)) +
  rremove("xlab") +
  theme(legend.position = "none") +
  scale_fill_manual(breaks = c("C-C", "C-H", "H-C", "H-H"), values=c("#ffffb2","#fecc5c","#fd8d3c", "#e31a1c")) +
  geom_boxplot() + facet_wrap(~Sex) + theme( strip.background = element_rect(color="white", fill="white", size=1.5, linetype="blank"), strip.text.x = element_text(size = 14, color = "white", face = "bold"))
O2DMEL

# D. pseudoobscura
O2DPSE=ggplot(respDPSE, aes(x=Treatment, y=O2pHr2, fill=Treatment)) + 
  ylim(0.5,4.2) +
  rremove("xlab") +
  rremove("ylab") +
  theme(axis.text.x=element_blank()) +
  theme(axis.text.y=element_blank()) +
  theme(legend.position = "none") +
  scale_fill_manual(breaks = c("C-C", "C-H", "H-C", "H-H"), values=c("#ffffb2","#fecc5c","#fd8d3c", "#e31a1c")) +
  geom_boxplot() + facet_wrap(~Sex) + theme( strip.background = element_rect(color="white", fill="white", size=1.5, linetype="blank"), strip.text.x = element_text(size = 14, color = "white", face = "bold"))
O2DPSE


## Panel D (CO2)
# D. melanogaster
CO2DMEL=ggplot(respDMEL, aes(x=Treatment, y=CO2pHr2, fill=Treatment)) + 
  theme(axis.title.y = element_text(size = 11)) +
  theme(axis.text.y = element_text(size = 10)) +
  ylab(bquote(ul ~ CO[2] ~ h^-1)) +
  ylim(0.7,3.8) +
  rremove("xlab") +
  theme(axis.text.x=element_blank()) + 
  theme(axis.text.x=element_blank()) +
  theme(legend.position = "none") +
  scale_fill_manual(breaks = c("C-C", "C-H", "H-C", "H-H"), values=c("#ffffb2","#fecc5c","#fd8d3c", "#e31a1c")) +
  geom_boxplot() + facet_wrap(~Sex) + theme( strip.background = element_rect(color="white", fill="white", size=1.5, linetype="blank"), strip.text.x = element_text(size = 14, color = "white", face = "bold"))
CO2DMEL


# D. pseudoobscura
CO2DPSE=ggplot(respDPSE, aes(x=Treatment, y=CO2pHr2, fill=Treatment)) + 
  ylim(0.7,3.8) +
  rremove("xlab") +
  rremove("ylab") +
  theme(axis.text.x=element_blank()) + 
  theme(axis.text.y=element_blank()) +
  theme(legend.position = "none") +
  scale_fill_manual(breaks = c("C-C", "C-H", "H-C", "H-H"), values=c("#ffffb2","#fecc5c","#fd8d3c", "#e31a1c")) +
  geom_boxplot() + facet_wrap(~Sex) + theme( strip.background = element_rect(color="white", fill="white", size=1.5, linetype="blank"), strip.text.x = element_text(size = 14, color = "white", face = "bold"))
CO2DPSE


## Panel E (Body mass)
# D. melanogaster
WDMEL=ggplot(respDMEL, aes(x=Treatment, y= Weight, fill=Treatment)) +
  ylim(0,1.5) +
  theme(axis.title.y = element_text(size = 11)) +
  theme(axis.text.y = element_text(size = 10)) +
  ylab(bquote(Body ~ mass ~ (mg))) +
  rremove("xlab") +
  theme(axis.text.x=element_blank()) + 
  theme(legend.position = "none") +
  geom_text(data = WSdmel, aes(y = Weight.mean, x = Treatment, label = group), vjust = -5, size = 4) +
  scale_fill_manual(breaks = c("C-C", "C-H", "H-C", "H-H"), values=c("#ffffb2","#fecc5c","#fd8d3c", "#e31a1c")) +
  geom_boxplot() + facet_wrap(~Sex) + theme( strip.background = element_rect(color="white", fill="white", size=1.5, linetype="blank"), strip.text.x = element_text(size = 14, color = "white", face = "bold"))
WDMEL

# D. pseudoobscura

WDPSE=ggplot(respDPSE, aes(x=Treatment, y=Weight, fill=Treatment)) + 
  ylim(0,1.5) +
  rremove("xlab") +
  rremove("ylab") +
  theme(axis.text.x=element_blank()) + 
  theme(axis.text.y=element_blank()) +
  theme(legend.position = "none") +
  scale_fill_manual(breaks = c("C-C", "C-H", "H-C", "H-H"), values=c("#ffffb2","#fecc5c","#fd8d3c", "#e31a1c")) +
  geom_text(data = WSdpse, aes(y = Weight.mean, x = Treatment, label = group), vjust = -6, size = 4) +
  geom_boxplot() + facet_wrap(~Sex) + theme( strip.background = element_rect(color="white", fill="white", size=1.5, linetype="blank"), strip.text.x = element_text(size = 14, color = "white", face = "bold"))
WDPSE

Combined2=ggarrange(ncol=1, CT_DMEL, RQ_DMEL, O2DMEL, CO2DMEL, WDMEL, common.legend = TRUE, font.label = list(size = 1, color = "black", face = "bold"), align = "hv", legend = "none")
Combined2=annotate_figure(Combined2, top = text_grob("D. melanogaster", color = "black", face="bold.italic", family="sans", size = 15))
Combined2

Combined4=ggarrange(ncol=1, CT_DPSE, RQ_DPSE, O2DPSE, CO2DPSE, WDPSE, common.legend = TRUE, font.label = list(size = 1, color = "black", face = "bold"), align = "hv", legend = "none")
Combined4=annotate_figure(Combined4, top = text_grob("D. pseudoobscura", color = "black", face="bold.italic", family="sans", size = 15))
Combined4

pdf("Figures/Fig4.pdf", width = 8.5, height = 18)
Combined5=ggarrange(ncol=2, Combined2, Combined4, common.legend = TRUE, font.label = list(size = 1, color = "black", face = "bold"), align = "hv", heights = c(0.3, 0.3, 0.3, 0.3), legend = "none")
Combined5
dev.off()



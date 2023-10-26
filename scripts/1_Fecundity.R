##Fecundity

##Import data set 
Eggs <- read_csv("rawdata/Eggs.csv", col_types = cols(Age = col_number(),
  Replicate = col_number(), Number_of_moms = col_number(),
  Total_eggs = col_number(), Number_of_adults = col_number(), 
  Number_necrosis_larvae = col_number(),
  Larvae = col_number(), Pupae = col_number()))
Eggs$Age <- as.factor(Eggs$Age)
View(Eggs)
str(Eggs)

##Adjust per fly
Eggs$Eggs_per_fly=Eggs$Total_eggs/Eggs$Number_of_moms
Eggs$Pupae_per_fly=Eggs$Pupae/Eggs$Number_of_moms
Eggs$Larvae_per_fly=Eggs$Larvae/Eggs$Number_of_moms
Eggs$Adults_per_fly=Eggs$Number_of_adults/Eggs$Number_of_moms


#And sum by day by replicate 
New_eggs=summaryBy(Total_eggs+Number_of_adults+Number_necrosis_larvae+Larvae+Pupae+Number_of_moms+Eggs_per_fly+Pupae_per_fly+Larvae_per_fly+Adults_per_fly~Replicate+Complete_Treatment+Species,data=Eggs,FUN=sum,na.rm=T)
New_eggs$Complete_Treatment <- as.factor(New_eggs$Complete_Treatment)
New_eggs$Replicate <- as.factor(New_eggs$Replicate)
New_eggs$Species <- as.factor(New_eggs$Species)
str(New_eggs)

##Subset by species for plots
EggsDPSE2<-subset(New_eggs, New_eggs$Species=="DPSE", na.rm=TRUE, select = c(Species, Complete_Treatment, Larvae_per_fly.sum, Pupae_per_fly.sum, Adults_per_fly.sum, Eggs_per_fly.sum))
EggsDMEL2<-subset(New_eggs, New_eggs$Species=="DMEL", na.rm=TRUE, select = c(Species, Complete_Treatment, Larvae_per_fly.sum, Pupae_per_fly.sum, Adults_per_fly.sum, Eggs_per_fly.sum))

##Statistical analysis per developmental stage
##Eggs
res.aov <- lmer(Eggs_per_fly.sum ~ (1|Replicate) + Species + Complete_Treatment + Species:Complete_Treatment, data = New_eggs)
as.data.frame(summary(res.aov)$coefficients)
write_csv2(as.data.frame(summary(res.aov)$coefficients), "Stats/Eggs_lmer.csv", col_names = T)
anova_res.aov <- Anova(res.aov)
write_csv2(as.data.frame(anova_res.aov), "Stats/Eggs_aov.csv", col_names = T)
lsm1 <- emmeans(res.aov, ~ Complete_Treatment:Species)
lsm1
summary(lsm1)
write_csv2(as.data.frame(lsm1), "Stats/emmeans_Eggs.csv")

##Extract significant letters
cld1 <- cld(lsm1, Letters = letters, alpha = 0.05)

##Subset significant letters for plots
cldDMEL1 <- subset(cld1, cld1$Species=="DMEL", na.rm=TRUE, select = c(Complete_Treatment, .group, emmean, SE))
cldDPSE1 <- subset(cld1, cld1$Species=="DPSE", na.rm=TRUE, select = c(Complete_Treatment, .group, emmean, SE))

##Pupae
res.aov2 <- lmer(Pupae_per_fly.sum ~ (1|Replicate) + Species + Complete_Treatment+ Species:Complete_Treatment, data = New_eggs)
as.data.frame(summary(res.aov2)$coefficients)
write_csv2(as.data.frame(summary(res.aov2)$coefficients), "Stats/Pupae_lmer.csv", col_names = T)
anova_res.aov2 <- Anova(res.aov2)
write_csv2(as.data.frame(anova_res.aov2), "Stats/Pupae_aov.csv", col_names = T)
lsm2 <- emmeans(res.aov2, ~ Complete_Treatment:Species)
write_csv2(as.data.frame(lsm2), "Stats/emmeans_Pupae.csv")

##Extract significant letters
cld2 <-cld(lsm2, Letters = letters, alpha = 0.05)

##Subset significant letters for plots
cldDMEL2 <- subset(cld2, cld2$Species=="DMEL", na.rm=TRUE, select = c(Complete_Treatment, .group, emmean, SE))
cldDPSE2 <- subset(cld2, cld2$Species=="DPSE", na.rm=TRUE, select = c(Complete_Treatment, .group, emmean, SE))

##Adults
res.aov3 <- lmer(Adults_per_fly.sum ~ (1|Replicate) + + Species + Complete_Treatment+ Species:Complete_Treatment, data = New_eggs)
as.data.frame(summary(res.aov2)$coefficients)
write_csv2(as.data.frame(summary(res.aov2)$coefficients), "Stats/Adult_lmer.csv", col_names = T)
anova_res.aov3 <- Anova(res.aov3)
write_csv2(as.data.frame(anova_res.aov3), "Stats/Adults_aov.csv", col_names = T)
lsm3 <- emmeans(res.aov3, ~ Species * Complete_Treatment)
write_csv2(summary(lsm3), "Stats/emmeans_Adults.csv")

##Extract significant letters
cld3 <-cld(lsm3, Letters = letters, alpha = 0.05)

##Subset significant letters for plots
cldDMEL3 <- subset(cld3, cld3$Species=="DMEL", na.rm=TRUE, select = c(Complete_Treatment, .group, emmean, SE))
cldDPSE3 <- subset(cld3, cld3$Species=="DPSE", na.rm=TRUE, select = c(Complete_Treatment, .group, emmean, SE))

###Plots
#Function to add N to boxplot Eggs
stat_box_data1 <- function(y, upper_limit = 310) {
  return( 
    data.frame(
      y = 0.8 * upper_limit,
      label = paste('n=', length(y))
    )
  )
}

##Subset for significant letters summarized
ES=summaryBy(Eggs_per_fly.sum~Complete_Treatment+Species,data=New_eggs,FUN = mean,na.rm=T, stringsAsFactors = FALSE)
E.summarized=merge(ES,cld1, stringsAsFactors = FALSE, na.rm=TRUE, check.names = F)
E.summarized$group=ifelse(E.summarized$.group==" a   ","a", ifelse(E.summarized$.group=="  bc ","bc", ifelse(E.summarized$.group=="  b  ","b", ifelse(E.summarized$.group=="   cd","cd", ifelse(E.summarized$.group=="    d","d", E.summarized$.group)))))
ESdmel = subset(E.summarized, E.summarized$Species=="DMEL", na.rm=TRUE, select = c(Complete_Treatment, Eggs_per_fly.sum.mean, group))
ESdpse = subset(E.summarized, E.summarized$Species=="DPSE", na.rm=TRUE, select = c(Complete_Treatment, Eggs_per_fly.sum.mean, group))

PS=summaryBy(Pupae_per_fly.sum~Complete_Treatment+Species,data=New_eggs,FUN = mean,na.rm=T)
P.summarized=merge(PS,cld2, all = T)
P.summarized$group=ifelse(P.summarized$.group=="    d ","d", ifelse(P.summarized$.group=="  bc  ","bc", ifelse(P.summarized$.group=="   c  ","c", ifelse(P.summarized$.group==" a    ","a", ifelse(P.summarized$.group==" ab   ","ab", ifelse(P.summarized$.group=="     e","e", P.summarized$.group))))))
PSdmel = subset(P.summarized, P.summarized$Species=="DMEL", na.rm=TRUE, select = c(Complete_Treatment, Pupae_per_fly.sum.mean, group))
PSdpse = subset(P.summarized, P.summarized$Species=="DMEL", na.rm=TRUE, select = c(Complete_Treatment, Pupae_per_fly.sum.mean, group))

AS=summaryBy(Adults_per_fly.sum~Complete_Treatment+Species,data=New_eggs,FUN = mean,na.rm=T)
A.summarized=merge(AS,cld3, all = T)
A.summarized$group=ifelse(A.summarized$.group=="   c ","c", ifelse(A.summarized$.group=="    d","d", ifelse(A.summarized$.group==" ab  ","ab", ifelse(A.summarized$.group=="  b  ","b", ifelse(A.summarized$.group==" a   ","a", A.summarized$.group)))))
ASdmel = subset(A.summarized, A.summarized$Species=="DMEL", na.rm=TRUE, select = c(Complete_Treatment, Adults_per_fly.sum.mean, group))
ASdpse = subset(A.summarized, A.summarized$Species=="DMEL", na.rm=TRUE, select = c(Complete_Treatment, Adults_per_fly.sum.mean, group))

##Figure 3 Panel A Left
DMEL3=ggplot(EggsDMEL2, aes(x=Complete_Treatment, y=Eggs_per_fly.sum, fill=Complete_Treatment))+ 
  scale_fill_manual(breaks = c("CC", "CHT", "HTC", "HTHT"), values=c("#ffffb2","#fecc5c","#fd8d3c", "#e31a1c")) +
  rremove("ylab") +
  rremove("xlab") +
  theme(axis.text.x=element_blank()) +
  ylim(0,280) +
  theme(legend.position = "none") +
  stat_summary(fun.data = stat_box_data1, geom = "text", fun = median, position = position_dodge(width = 1), size = 5) +
  geom_text(data = ESdmel, aes(y = Eggs_per_fly.sum.mean, x = Complete_Treatment, label = group), vjust = -2, position = position_dodge(width = 1), size = 5) +
  geom_boxplot(position = position_dodge(width = 1)) 
DMEL3

##Figure 3 Panel A Right
DPSE3=ggplot(EggsDPSE2, aes(x=Complete_Treatment, y=Eggs_per_fly.sum, fill=Complete_Treatment)) + 
  scale_fill_manual(breaks = c("CC", "CHT", "HTC", "HTHT"), values=c("#ffffb2","#fecc5c","#fd8d3c", "#e31a1c")) +
  rremove("ylab") +
  rremove("xlab") +
  theme(legend.position = "none") +
  theme(axis.text.x=element_blank()) +
  theme(axis.text.y=element_blank()) +
  ylim(0,280) +
  stat_summary(fun.data = stat_box_data1, geom = "text", fun = median, position = position_dodge(width = 0.8), size = 5) +
  geom_text(data = ESdpse, aes(y = Eggs_per_fly.sum.mean, x = Complete_Treatment, label = group), vjust = -3.5, size = 5) +
  geom_boxplot() 
DPSE3

##Figure 3 Panel B Left
Pupae1DMEL2=ggplot(EggsDMEL2, aes(x=Complete_Treatment, y=Pupae_per_fly.sum, fill=Complete_Treatment)) + 
  scale_fill_manual(breaks = c("CC", "CHT", "HTC", "HTHT"), values=c("#ffffb2","#fecc5c","#fd8d3c", "#e31a1c")) +
  rremove("ylab") +
  rremove("xlab") +
  theme(axis.text.x=element_blank()) +
  ylim(0,110) +
  theme(legend.position = "none") +
  geom_text(data = PSdmel, aes(y = Pupae_per_fly.sum.mean, x = Complete_Treatment, label = group), vjust = -2, size = 5) +
  geom_boxplot()
Pupae1DMEL2

##Figure 3 Panel B Right
Pupae1DPSE2=ggplot(EggsDPSE2, aes(x=Complete_Treatment, y=Pupae_per_fly.sum, fill=Complete_Treatment)) + 
  scale_fill_manual(breaks = c("CC", "CHT", "HTC", "HTHT"), values=c("#ffffb2","#fecc5c","#fd8d3c", "#e31a1c")) +
  rremove("ylab") +
  rremove("xlab") +
  theme(axis.text.x=element_blank()) +
  theme(axis.text.y=element_blank()) +
  ylim(0,110) +
  theme(legend.position = "none") +
  geom_text(data = PSdpse, aes(y = Pupae_per_fly.sum.mean, x = Complete_Treatment, label = group), vjust = -8, size = 5) +
  geom_boxplot() 
Pupae1DPSE2

##Figure 3 Panel C Left
Adult1DMEL2=ggplot(EggsDMEL2, aes(x=Complete_Treatment, y=Adults_per_fly.sum, fill=Complete_Treatment)) + 
  scale_fill_manual(breaks = c("CC", "CHT", "HTC", "HTHT"), values=c("#ffffb2","#fecc5c","#fd8d3c", "#e31a1c")) +
  rremove("ylab") +
  rremove("xlab") +
  theme(axis.text.x=element_blank()) +
  ylim(0,90) +
  theme(legend.position = "none") +
  geom_text(data = ASdmel, aes(y = Adults_per_fly.sum.mean+10, x = Complete_Treatment, label = group), vjust = -0.5, size = 5) +
  geom_boxplot() 
Adult1DMEL2

##Figure 3 Panel C Right
Adult1DPSE2=ggplot(EggsDPSE2, aes(x=Complete_Treatment, y=Adults_per_fly.sum, fill=Complete_Treatment)) + 
  scale_fill_manual(breaks = c("CC", "CHT", "HTC", "HTHT"), values=c("#ffffb2","#fecc5c","#fd8d3c", "#e31a1c")) +
  rremove("ylab") +
  rremove("xlab") +
  theme(axis.text.x=element_blank()) +
  theme(axis.text.y=element_blank()) +
  ylim(0,90) +
  theme(legend.position = "none") +
  geom_text(data = ASdpse, aes(y = Adults_per_fly.sum.mean, x = Complete_Treatment, label =group), vjust = -9, size = 5) +
  geom_boxplot() 
Adult1DPSE2

##Combined figure 

Combined3=ggarrange(ncol=1, DMEL3, Pupae1DMEL2, Adult1DMEL2, common.legend = TRUE, font.label = list(size = 1, color = "black", face = "bold"), align = "hv", heights = c(1.1, 0.8, 0.8), legend = "none")
Combined3=annotate_figure(Combined3, top = text_grob("D. melanogaster", color = "black", face="bold.italic", family="sans", size = 15))
Combined3

Combined1=ggarrange(ncol=1, DPSE3, Pupae1DPSE2, Adult1DPSE2, common.legend = TRUE, font.label = list(size = 1, color = "black", face = "bold"), align = "hv", heights = c(1.1, 0.8, 0.8), legend = "none")
Combined1=annotate_figure(Combined1, top = text_grob("D. pseudoobscura", color = "black", face="bold.italic", family="sans", size = 15))
Combined1

pdf("Figures/Fig3.pdf", width = 8.5, height = 8.5)
Combined2=ggarrange(ncol=2, Combined3, Combined1, common.legend = TRUE, font.label = list(size = 1, color = "black", face = "bold"), align = "hv", heights = c(0.3, 0.3, 0.3, 0.3), legend = "none")
Combined2
dev.off()


#####Supplemental Data with fecundity by day
###Data by day of both species
EggsDM<-subset(Eggs, Eggs$Species=="DMEL", na.rm=TRUE, select = c(Species, Complete_Treatment, Larvae_per_fly, Pupae_per_fly, Adults_per_fly, Eggs_per_fly, Age))
EggsDP<-subset(Eggs, Eggs$Species=="DPSE", na.rm=TRUE, select = c(Species, Complete_Treatment, Larvae_per_fly, Pupae_per_fly, Adults_per_fly, Eggs_per_fly, Age))

##DMEL
pdf("Figures/Supplementals1.pdf")
Day=ggplot(EggsDM, aes(x=Age, y=Eggs_per_fly, fill=Complete_Treatment)) + 
  scale_fill_manual(breaks = c("CC", "CHT", "HTC", "HTHT"), values=c("#ffffb2","#fecc5c","#fd8d3c", "#e31a1c")) +
  ggtitle("D. melanogaster")+
  geom_boxplot() 
Day
Day2=ggplot(EggsDM, aes(x=Age, y=Pupae_per_fly, fill=Complete_Treatment)) + 
  scale_fill_manual(breaks = c("CC", "CHT", "HTC", "HTHT"), values=c("#ffffb2","#fecc5c","#fd8d3c", "#e31a1c")) +
  ggtitle("D. melanogaster")+
  geom_boxplot() 
Day2
Day3=ggplot(EggsDM, aes(x=Age, y=Adults_per_fly, fill=Complete_Treatment)) + 
  scale_fill_manual(breaks = c("CC", "CHT", "HTC", "HTHT"), values=c("#ffffb2","#fecc5c","#fd8d3c", "#e31a1c")) +
  ylim(0,12) + ggtitle("D. melanogaster")+
  geom_boxplot()
Day3

##DPSE
Day=ggplot(EggsDP, aes(x=Age, y=Eggs_per_fly, fill=Complete_Treatment)) + 
  scale_fill_manual(breaks = c("CC", "CHT", "HTC", "HTHT"), values=c("#ffffb2","#fecc5c","#fd8d3c", "#e31a1c")) +
  ggtitle("D. pseudoobscura")+geom_boxplot() 
Day
Day2=ggplot(EggsDP, aes(x=Age, y=Pupae_per_fly, fill=Complete_Treatment)) + 
  scale_fill_manual(breaks = c("CC", "CHT", "HTC", "HTHT"), values=c("#ffffb2","#fecc5c","#fd8d3c", "#e31a1c")) +
  ggtitle("D. pseudoobscura") +
  geom_boxplot() 
Day2
Day3=ggplot(EggsDP, aes(x=Age, y=Adults_per_fly, fill=Complete_Treatment)) + 
  scale_fill_manual(breaks = c("CC", "CHT", "HTC", "HTHT"), values=c("#ffffb2","#fecc5c","#fd8d3c", "#e31a1c")) +
  ggtitle("D. pseudoobscura") +
  geom_boxplot()
Day3
dev.off()

######Load required packages
library(ggplot2)
library(ggpubr)
theme_set(theme_pubr())
library(doBy)
library(dplyr)
library(tidyr)
library(readr)
library(ggthemes)
library(ggprism)
library(multcompView)
library(multcomp)
library(ragg)
library(agricolae)
library(car)
library(plyr)
library(memisc)
library(tidyverse) # ggplot & helper functions
library(scales)
library(data.table)
library(ggnewscale)
library(lme4)
library(emmeans)

##### RUN ANALYSIS
###Codes are independent from each other, and do not need to be run in this specific order
##Fecundity
source("scripts/1_Fecundity.R")

##Physiology
source("scripts/2_Physiology.R")

##Oogenesis
source("scripts/3_Stages_OOgenesis.R")

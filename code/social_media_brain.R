library(dplyr)
library(lme4)
library(interactions)
library(jtools)
library(ggsci)
library(weights)
library(EMAtools)
library(ggplot2)
library(Rcpp)
library(lm.beta)
library(psych)
library(corrplot)

rm(list = ls())

#load data
df_wide <- read.csv("../data/df_wide.csv")
df_long <- read.csv("../data/df_long.csv")

length(unique(df_wide$pID)) #54
length(unique(df_long$pID)) #54

df_wide$condition <- relevel(df_wide$condition, ref='control')
df_long$condition <- relevel(df_long$condition, ref='control')

###################
##### METHODS #####
###################

#demography
summary(df_wide$age)
sd(df_wide$age)
table(df_wide$gender)
table(df_wide$race)
table(df_wide$race_numeric)
summary(df_wide$race_numeric) #no NAs


#standardize social media variable (time-varying predictor variable)
df_long$SocialMediaz = with(df_long, ave(SocialMedia, pID, FUN=stdz)) #raw survey scores
df_long$social_media_minz = with(df_long, ave(social_media_min, pID, FUN=stdz)) #converted to indicate minutes


###################
##### RESULTS #####
###################


###########################################
##### 4.1. Time spent on social media #####
###########################################

## days of at least one use out of 28 days
summary(df_wide$soc_med_freq_daily_total)  
sd(df_wide$soc_med_freq_daily_total, na.rm=T)

## average time spent per use
summary(df_wide$soc_med_min_mean)  
sd(df_wide$soc_med_min_mean, na.rm=T)
###range
summary(df_long$social_media_min)


#############################################################################################################
##### 4.2. Associations between functional connectivity, internalizing symptoms, and emotion regulation #####
#############################################################################################################

#depression and anxiety prevelance
table(df_wide$prescan_CESD_sum)
summary(df_wide$prescan_CESD_sum)
sd(df_wide$prescan_CESD_sum)

table(df_wide$prescan_STAI_sum)
summary(df_wide$prescan_STAI_sum)
sd(df_wide$prescan_STAI_sum, na.rm=T)

# frontoparietal connectivity, depression, and anxiety
test = lm(prescan_CESD_sum ~ Fronto.parietal_Task_Control_withinFC
           + age + as.factor(gender) + as.factor(race_numeric) + rung_community + condition, df_wide)
summ(test, digits = 3, confint = TRUE)

test = lm(prescan_STAI_sum ~ Fronto.parietal_Task_Control_withinFC
          + age + as.factor(gender) + as.factor(race_numeric) + rung_community + condition, df_wide)
summ(test, digits = 3, confint = TRUE)

test = lm(prescan_DERS_mean ~ Fronto.parietal_Task_Control_withinFC
          + age + as.factor(gender) + as.factor(race_numeric) + rung_community + condition, df_wide)
summ(test, digits = 3, confint = TRUE)

####################
##### Figure 2 #####
####################

# CESD ~ Frontoparietal
ggplot (df_wide,aes(Fronto.parietal_Task_Control_withinFC, prescan_CESD_sum))+
  geom_point(shape=21,size=3,colour="white",fill="#2080F6")+
  geom_smooth(method=lm,color='red', size=.5)+
  theme(text=element_text(size=15))+
  xlab("Resting state frontoparietal FC (z')") + ylab ("Depressive symptoms (CES-D)" ) +
  theme_bw()+
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(), axis.line=element_line(colour="black"))

# STAI~ Frontoparietal
ggplot (df_wide,aes(Fronto.parietal_Task_Control_withinFC, prescan_STAI_sum))+
  geom_point(shape=21,size=3,colour="white",fill="#2080F6")+
  geom_smooth(method=lm,color='red', size=.5)+
  theme(text=element_text(size=15))+
  xlab("Resting state frontoparietal FC (z')") + ylab ("Anxiety symptoms (STAI)" ) +
  theme_bw()+
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(), axis.line=element_line(colour="black"))


# DERS ~ Frontoparietal
ggplot (df_wide,aes(Fronto.parietal_Task_Control_withinFC, prescan_DERS_mean))+
  geom_point(shape=21,size=3,colour="white",fill="#2080F6")+
  geom_smooth(method=lm,color='red', size=.5)+
  theme(text=element_text(size=15))+
  xlab("Resting state frontoparietal FC (z')") + ylab ("Difficulty regulating emotion (DERS)" ) +
  theme_bw()+
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(), axis.line=element_line(colour="black"))


##############################################################################################################
##### 4.3. Relationships among social media use, functional connectivity, and subsequent negative affect #####
##############################################################################################################

test = lmer(NegativeEMA ~ social_media_minz  * Fronto.parietal_Task_Control_withinFC + (social_media_minz  | pID) 
             + age + as.factor(gender) + as.factor(race_numeric) + rung_community + condition, df_long)
summ(test, digits = 3, confint = TRUE)
summary(test)


#simple slopes analysis
interactions::sim_slopes(model=test, pred=social_media_minz, modx=Fronto.parietal_Task_Control_withinFC,
                         digits = getOption("jtools-digits", default = 3),
                         confint = T)


#Figure 3

test = lmer(NegativeEMA ~ social_media_minz  * Fronto.parietal_Task_Control_withinFC + (social_media_minz  | pID), df_long)


interact_plot(test, pred = social_media_minz,  modx = Fronto.parietal_Task_Control_withinFC, interval = TRUE,
              x.label = "Social media use", y.label = "Negative affect",
              int.type = "confidence", int.width = 0.95)




#####################################
##### Supplementary Information #####
#####################################

#SI2. Subregions within the frontoparietal system



#SI3. Results using raw scores for the time spent on social media measure

test = lmer(NegativeEMA ~ SocialMediaz  * Fronto.parietal_Task_Control_withinFC + (SocialMediaz  | pID) 
            + age + as.factor(gender) + as.factor(race_numeric) + rung_community + condition, df_long)
summ(test, digits = 3, confint = TRUE)
summary(test)

#SI4 Current positive affect

test = lmer(PositiveEMA ~ social_media_minz  * Fronto.parietal_Task_Control_withinFC + (social_media_minz  | pID) 
            + age + as.factor(gender) + as.factor(race_numeric) + rung_community + condition, df_long)
summ(test, digits = 3, confint = TRUE)
summary(test)


#SI5. Difficulties in Emotion Regulation (DERS) subscales

df_wide$FC = df_wide$Fronto.parietal_Task_Control_withinFC 
ders = df_wide %>%
  dplyr::select("FC", grep("prescan_DERS_", colnames(df_wide)))

names(ders) = sub("prescan_DERS_*", "", names(ders))

M <- Hmisc::rcorr(as.matrix(ders))
corrplot(M$r, p.mat = M$P, insig = "label_sig",
         sig.level = c(.001, .05), pch.cex=0.9, pch.col = "white",
         method="color", type="lower")

cor.test(ders$FC, ders$goals)
cor.test(ders$FC, ders$strategies)

#SI6. Temporal relationships between time spent on social media and affect

df_long$NegativeEMA_previousz = with(df_long, ave(NegativeEMA_previous, pID, FUN=stdz)) 
df_long$PositiveEMA_previousz = with(df_long, ave(PositiveEMA_previous, pID, FUN=stdz)) 

test = lmer(SocialMedia  ~ NegativeEMA_previousz * Fronto.parietal_Task_Control_withinFC   + (1 | pID) 
            + age + as.factor(gender) + as.factor(race_numeric) + rung_community + condition, df_long)
summ(test, digits = 3, confint = TRUE)


test = lmer(SocialMedia  ~ PositiveEMA_previousz * Fronto.parietal_Task_Control_withinFC  + (PositiveEMA_previousz | pID) 
            + age + as.factor(gender) + as.factor(race_numeric) + rung_community + condition, df_long)
summ(test, digits = 3, confint = TRUE)


#SI7.Results not controlling for potential covariates
test = lm(prescan_CESD_sum ~ Fronto.parietal_Task_Control_withinFC, df_wide)
summ(test, digits = 3, confint = TRUE)

test = lm(prescan_STAI_sum ~ Fronto.parietal_Task_Control_withinFC, df_wide)
summ(test, digits = 3, confint = TRUE)

test = lm(prescan_DERS_mean ~ Fronto.parietal_Task_Control_withinFC, df_wide)
summ(test, digits = 3, confint = TRUE)

test = lmer(NegativeEMA ~ social_media_minz  * Fronto.parietal_Task_Control_withinFC + (social_media_minz  | pID), df_long)
summ(test, digits = 3, confint = TRUE)



#extra

#define a function to get standardized beta from lmer
lmer.beta <- function(object) {
  sdy <- sd(attr(object, "resp")$y) # the y values are now in the 'y' slot 
  ###                                 of the resp attribute
  sdx <- apply(attr(object, "pp")$X, 2, sd) # And the X matriz is in the 'X' slot of the pp attr
  sc <- fixef(object)*sdx/sdy
  #mimic se.ranef from pacakge "arm"
  se.fixef <- function(obj) as.data.frame(summary(obj)[10])[,2] # last change - extracting 
  ##             the standard errors from the summary
  se <- se.fixef(object)*sdx/sdy
  return(data.frame(stdcoef=sc, stdse=se))
}
###### Comp Lakes Biotic Models:

# load needed packages 
library(vegan)
library(sp) #space
library(ggplot2)
library(lme4)
library(lmerTest)# for p-value
library(MuMIn) # forr squared
library(PerformanceAnalytics) # correlation analysis 

setwd("~/Documents/Zoop data for papers")

# read in data file 
d <- read.csv("comp_lakes_master_data_2.csv", header=T)
k <- read.csv("COMPLAKES_2016data_a.csv", header=T) # no inlet outlet data
l <- read.csv("2016_Comp_Lakes_nutrients.csv", header=T)
cl <- read.csv("cl_averages.csv", header=T)

d$Site
# exclude mud lake because it is shallow seepage lake 
d2 <- d[c(1:37, 41:55),]
summary(d2$Site)

#correlate lake size variables
d$Circum_m
hist(d2$Circum_m) 
hist(log10(d2$Circum_m)) #normal

hist(d2$Sur_Area_m2) 
hist(log10(d2$Sur_Area_m2))#somewhat normal

hist(d2$Max_Depth) 
hist(log10(d2$Max_Depth))#somewhat normal

hist(d2$NDVI) 
hist(log10(d2$NDVI +1)) # worse

cor.test(log10(d2$Circum_m), (log10(d2$Sur_Area_m2)))
# cor=0.8955285
cor.test(log10(d$Max_Depth), (log10(d$Sur_Area_m2)))
# cor=0.5509976  
cor.test(log10(d$Max_Depth), log10(d$Circum_m))
# cor= 0.5389244 

cor.test(log10(d$NDVI), d$Elevation)
# -0.6601449

d$rho_change.hypo_sur.


######################################
###### Zooplankton Richness #########
#####################################
# univariate 
hist(d2$zoop_rich)
uni_zrich.1 <- glm(zoop_rich ~ log10(Circum_m), data = d2)
summary(uni_zrich.1) # better

uni_zrich.2 <- glm(zoop_rich ~ log10(Sur_Area_m2), data = d2)
summary(uni_zrich.2) # sig

uni_zrich.3 <- glm(zoop_rich ~ log10(Max_Depth), data = d2)
summary(uni_zrich.3) #  

uni_zrich.3 <- glm(zoop_rich ~ NDVI, data = d2)
summary(uni_zrich.3) # sig

uni_zrich.3 <- glm(zoop_rich ~ Elevation, data = d2)
summary(uni_zrich.3) # sig

uni_zrich.3 <- glm(zoop_rich ~ Fish, data = d2)
summary(uni_zrich.3)

uni_zrich.3 <- glm(zoop_rich ~ lake_order, data = d2)
summary(uni_zrich.3) # sig 

uni_zrich.3 <- glm(zoop_rich ~ sample_number, data = d2)
summary(uni_zrich.3) # sig 


#### model with RE: site and watershed
# build null model (w/o elevation):

#null
zrich.null <- lmer(zoop_rich ~  Fish  + scale(lake_order) + scale(log10(Circum_m))
                   + scale(log10(Max_Depth)) + (1| Site) +  (1| Watershed) + (1| visit_re), data = d2, REML=FALSE)
summary(zrich.null)

#global model
zrich.global <- lmer(zoop_rich ~  scale(Elevation) + Fish + scale(lake_order) + scale(log10(Circum_m))
                     + scale(log10(Max_Depth)) + (1| Site) +  (1| Watershed) + (1| visit_re), data = d2, REML=FALSE)

summary(zrich.global)

#elevation
zrich.mod <- lmer(zoop_rich ~  scale(Elevation) + scale(lake_order)
                  + (1| Site) +  (1| Watershed) + (1| visit_re), data = d2, REML=FALSE)

summary(zrich.mod)

AIC(zrich.null, zrich.global, zrich.mod)
hist(residuals(zrich.mod))

 # better with elevation 
r.squaredGLMM(zrich.null) 
r.squaredGLMM(zrich.global)
r.squaredGLMM(zrich.mod)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Describe:
# Zooplankton taxa richness was negatively correlated with:
# elevation (-0.65 ± 0.62) 
# and posivielty correlated with
# lake order (2.08 ± 0.81, p<0.03) 
# R2m:0.312, R2c:0.884 




######## total density ############
##################################

hist(d2$total_rich)
uni_zrich.1 <- glm(total_rich ~ log10(Circum_m), data = d2)
summary(uni_zrich.1) # better

uni_zrich.2 <- glm(total_rich ~ log10(Sur_Area_m2), data = d2)
summary(uni_zrich.2) 

uni_zrich.3 <- glm(total_rich ~ log10(Max_Depth), data = d2)
summary(uni_zrich.3) #  

uni_zrich.3 <- glm(total_rich ~ NDVI, data = d2)
summary(uni_zrich.3) # ***

uni_zrich.3 <- glm(total_rich ~ Elevation, data = d2)
summary(uni_zrich.3) # sig

uni_zrich.3 <- glm(total_rich ~ Fish, data = d2)
summary(uni_zrich.3)

uni_zrich.3 <- glm(total_rich ~ lake_order, data = d2)
summary(uni_zrich.3) # sig 

uni_zrich.3 <- glm(total_rich ~ sample_number, data = d2)
summary(uni_zrich.3) # sig 

#### model with RE: site and watershed
# build null model (w/o elevation):

#null
trich.null <- lmer(total_rich ~  Fish  + scale(lake_order) + scale(log10(Circum_m))
                   + scale(log10(Max_Depth)) + (1| Site) +  (1| Watershed) + (1| visit_re), data = d2, REML=FALSE)
summary(trich.null)

#global model
trich.global <- lmer(total_rich ~  scale(Elevation) + Fish + scale(lake_order) + scale(log10(Circum_m))
                     + scale(log10(Max_Depth)) + (1| Site) +  (1| Watershed) + (1| visit_re), data = d2, REML=FALSE)

summary(trich.global)

#elevation
trich.mod <- lmer(total_rich ~  scale(Elevation) + Fish  + scale(log10(Circum_m))
                  + scale(log10(Max_Depth)) + (1| Site) +  (1| Watershed) + (1| visit_re), data = d2, REML=FALSE)

summary(trich.mod)

AIC(trich.null, trich.global, trich.mod)
hist(residuals(trich.mod))
hist(residuals(trich.global))

# better with elevation 
r.squaredGLMM(trich.null) 
r.squaredGLMM(trich.global) # best model 
r.squaredGLMM(trich.mod)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Describe:
# Total taxa richness was negatively correlated with:
# elevation (-4.27 ± 0.79, p<4.36e-06)
# Fish (-8.43 ± 1.76, p< 2.97e-05)
# log 10 (Max depth) (-3.55 ± 1.01, p< 0.002)
# and positively correlated with:
# lake order (1.79 ± 1.27) 
# log 10 (Circumference) (3.95 ± 0.76, p< 8.10e-06)

# R2m:0.713, R2c:0.713 







####################################
###### Zooplankton density #########
#####################################
hist(d$zoop_den)
hist(log10(d$zoop_den +1))

uni_den.1 <- glm(log10(zoop_den +1) ~ scale(log10(Circum_m)), data = d)
summary(uni_den.1) 

uni_den.1 <- glm(log10(zoop_den +1) ~ scale(log10(Sur_Area_m2)), data = d)
summary(uni_den.1) 

uni_den.1 <- glm(log10(zoop_den +1) ~ scale(log10(Max_Depth)), data = d)
summary(uni_den.1) #  

uni_den.1 <- glm(log10(zoop_den +1) ~ scale(NDVI), data = d)
summary(uni_den.1) # ***

uni_den.1 <- glm(log10(zoop_den +1) ~ scale(Elevation), data = d)
summary(uni_den.1) # sig

uni_den.1 <- glm(log10(zoop_den +1) ~ Fish, data = d)
summary(uni_den.1) # **

uni_den.1 <- glm(log10(zoop_den +1) ~ scale(lake_order), data = d)
summary(uni_den.1) # sig 

uni_den.1 <- glm(log10(zoop_den +1) ~ scale(sample_number), data = d)
summary(uni_den.1) # sig 


#### model with RE: site and watershed
#null
zden.null <- lmer(log10(zoop_den +1) ~  Fish  + scale(lake_order) + scale(log10(Circum_m))
                   + scale(log10(Max_Depth)) + (1| Site) +  (1| Watershed) + (1| visit_re), data = d2, REML=FALSE)
summary(zden.null)

#global model
zden.global <- lmer(log10(zoop_den +1) ~  scale(Elevation) + Fish + scale(lake_order) + scale(log10(Circum_m))
                     + scale(log10(Max_Depth)) + (1| Site) +  (1| Watershed) + (1| visit_re), data = d2, REML=FALSE)

summary(zden.global)

#elevation
zden.mod <- lmer(log10(zoop_den +1) ~  scale(Elevation) + Fish  + scale(log10(Circum_m))
                 + scale(log10(Max_Depth)) + (1| Site) +  (1| Watershed) + (1| visit_re), data = d2, REML=FALSE)

summary(zden.mod)

AIC(zden.null, zden.global, zden.mod)
hist(residuals(zden.null))
hist(residuals(zden.global))
hist(residuals(zden.mod))

# better with elevation 
r.squaredGLMM(zden.null) 
r.squaredGLMM(zden.global) # best model 
r.squaredGLMM(zden.mod)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Describe:
# Zooplankton density was negatively correlated with:
# elevation (-4.27 ± 0.79, p<4.36e-06)
# Fish (-8.43 ± 1.76, p< 2.97e-05)
# lake order (1.79 ± 1.27) 
# log 10 (Circumference) (3.95 ± 0.76, p< 8.10e-06)
# log 10 (Max depth) (-3.55 ± 1.01, p< 0.002)
# R2m:0.713, R2c:0.713 
# N= 52




###### Zooplankton biovolume #########
######################################
hist(d$BV_zoop)
hist(log10(d$BV_zoop +1))

uni_den.1 <- glm(log10(BV_zoop +1) ~ scale(log10(Circum_m)), data = d2)
summary(uni_den.1) 

uni_den.1 <- glm(log10(BV_zoop +1) ~ scale(log10(Sur_Area_m2)), data = d2)
summary(uni_den.1) 

uni_den.1 <- glm(log10(BV_zoop +1) ~ scale(log10(Max_Depth)), data = d2)
summary(uni_den.1) #  

uni_den.1 <- glm(log10(BV_zoop +1) ~ scale(NDVI), data = d2)
summary(uni_den.1) # ***

uni_den.1 <- glm(log10(BV_zoop +1) ~ scale(Elevation), data = d2)
summary(uni_den.1) # sig

uni_den.1 <- glm(log10(BV_zoop +1) ~ Fish, data = d2)
summary(uni_den.1) # **

uni_den.1 <- glm(log10(BV_zoop +1) ~ scale(lake_order), data = d2)
summary(uni_den.1) # sig 

uni_den.1 <- glm(log10(BV_zoop +1) ~ scale(sample_number), data = d2)
summary(uni_den.1) # sig 


# null model
zBV.null <- lmer(log10(BV_zoop +1) ~  Fish  + scale(lake_order) + scale(log10(Circum_m))
                     + scale(log10(Max_Depth)) + (1| Site) +  (1| Watershed) + (1| visit_re), data = d2, REML=FALSE)
summary(zBV.null)

# global model
zBV.global <- lmer(log10(BV_zoop +1) ~ scale(Elevation)+ Fish  + scale(lake_order) + scale(log10(Circum_m))
                 + scale(log10(Max_Depth)) + (1| Site) +  (1| Watershed) + (1| visit_re), data = d2, REML=FALSE)
summary(zBV.global)

#elevation model
zBV.mod <- lmer(log10(BV_zoop +1) ~ scale(Elevation) + scale(lake_order)
                   + (1| Site) +  (1| Watershed) + (1| visit_re), data = d2, REML=FALSE)
summary(zBV.mod)

AIC(zBV.null, zBV.global, zBV.mod)
hist(residuals(zBV.null))
hist(residuals(zBV.global))
hist(residuals(zBV.mod))

# better with elevation 
r.squaredGLMM(zBV.null) 
r.squaredGLMM(zBV.global) # best model 
r.squaredGLMM(zBV.mod)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Describe:
# Zooplankton density was negatively correlated with:
# elevation (-0.02 ± 0.09)
# and positively with
# lake order (0.09 ± 0.07) 
# R2m:0.053, R2c:0.731 




###### Ratio Zooplankton/Phytoplankton biovolume #########
##########################################################

hist(d$ratio_ZBV_PBV)
hist(asin(sqrt(d$ratio_ZBV_PBV)))

uni_Z_PBV.1 <- glm(asin(sqrt(ratio_ZBV_PBV)) ~ scale(log10(Circum_m)), data = d)
summary(uni_Z_PBV.1) 

uni_Z_PBV.1 <- glm(asin(sqrt(ratio_ZBV_PBV)) ~ scale(log10(Sur_Area_m2)), data = d)
summary(uni_Z_PBV.1) 

uni_Z_PBV.1 <- glm(asin(sqrt(ratio_ZBV_PBV)) ~ scale(log10(Max_Depth)), data = d)
summary(uni_Z_PBV.1) #  

uni_Z_PBV.1 <- glm(asin(sqrt(ratio_ZBV_PBV)) ~ scale(NDVI), data = d)
summary(uni_Z_PBV.1) # ***

uni_Z_PBV.1 <- glm(asin(sqrt(ratio_ZBV_PBV)) ~ scale(Elevation), data = d)
summary(uni_Z_PBV.1) # sig

uni_Z_PBV.1 <- glm(asin(sqrt(ratio_ZBV_PBV)) ~ Fish, data = d)
summary(uni_Z_PBV.1) # **

uni_Z_PBV.1 <- glm(asin(sqrt(ratio_ZBV_PBV)) ~ scale(lake_order), data = d)
summary(uni_Z_PBV.1) # sig 

uni_Z_PBV.1 <- glm(asin(sqrt(ratio_ZBV_PBV)) ~ scale(sample_number), data = d)
summary(uni_Z_PBV.1) # sig 


## Build null model 
d$Sampling_interval
zpBV.null <- lmer(asin(sqrt(ratio_ZBV_PBV)) ~ Fish  + scale(lake_order) + scale(log10(Circum_m))
                      + scale(log10(Max_Depth)) + (1| Site) +  (1| Watershed) + (1| visit_re), data = d, REML=FALSE)
summary(zpBV.null)

zpBV.global <- lmer(asin(sqrt(ratio_ZBV_PBV)) ~ scale(Elevation)+ Fish  + scale(lake_order) + scale(log10(Circum_m))
                      + scale(log10(Max_Depth)) + (1| Site) +  (1| Watershed) + (1| visit_re), data = d, REML=FALSE)
summary(zpBV.global)

zpBV.mod <- lmer(asin(sqrt(ratio_ZBV_PBV)) ~ scale(Elevation)+ Fish 
                    + (1| Site) +  (1| Watershed) + (1| visit_re), data = d, REML=FALSE)
summary(zpBV.mod)


AIC(zpBV.null, zpBV.global, zpBV.mod)
hist(residuals(zBV.null))
hist(residuals(zBV.global))
hist(residuals(zpBV.mod)) # no difference between null and elevation models


r.squaredGLMM(zpBV.null) 
r.squaredGLMM(zpBV.global) # best model 
r.squaredGLMM(zpBV.mod)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#



###### Proportion of gravid zooplankton #########
#################################################

hist(d$prop_gravid)
hist(asin(sqrt(d$prop_gravid)))

uni_propg <- glm(asin(sqrt(prop_gravid)) ~ scale(log10(Circum_m)), data = d)
summary(uni_propg) 

uni_propg <- glm(asin(sqrt(prop_gravid)) ~ scale(log10(Sur_Area_m2)), data = d)
summary(uni_propg) 

uni_propg <- glm(asin(sqrt(prop_gravid)) ~ scale(log10(Max_Depth)), data = d)
summary(uni_propg) #  

uni_propg <- glm(asin(sqrt(prop_gravid)) ~ scale(NDVI), data = d)
summary(uni_propg) # ***

uni_propg <- glm(asin(sqrt(prop_gravid)) ~ scale(Elevation), data = d)
summary(uni_propg) # sig

uni_propg <- glm(asin(sqrt(prop_gravid)) ~ Fish, data = d)
summary(uni_propg) # **

uni_propg <- glm(asin(sqrt(prop_gravid)) ~ scale(lake_order), data = d)
summary(uni_propg) # sig 

uni_propg <- glm(asin(sqrt(prop_gravid)) ~ scale(sample_number), data = d)
summary(uni_propg) # sig 


# null model
zgrad.null <- lmer(asin(sqrt(prop_gravid)) ~ Fish  + scale(lake_order) + scale(log10(Circum_m))
                       + scale(log10(Max_Depth)) + (1| Site) +  (1| Watershed) + (1| visit_re), data = d, REML=FALSE)
summary(zgrad.null)


zgrad.global <- lmer(asin(sqrt(prop_gravid)) ~  scale(Elevation) + Fish  + scale(lake_order) + 
                         scale(log10(Circum_m)) + scale(log10(Max_Depth)) + (1| Site) +  
                         (1| Watershed) + (1| visit_re), data = d, REML=FALSE)
summary(zgrad.global)


zgrad.mod <- lmer(asin(sqrt(prop_gravid)) ~  scale(Elevation) + Fish + (1| Site) +  
                       (1| Watershed) + (1| visit_re), data = d, REML=FALSE)
summary(zgrad.mod)

AIC(zgrad.null, zgrad.global, zgrad.mod)
hist(residuals(zgrad.null))
hist(residuals(zgrad.global))
hist(residuals(zgrad.mod)) # no difference between null and elevation models


r.squaredGLMM(zgrad.null) 
r.squaredGLMM(zgrad.global) # best model 
r.squaredGLMM(zgrad.mod)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Describe:
# Zooplankton fecundity was positively correlated with:
# elevation (0.08 ± 0.03, p<0.03)
# and negatively with
# fish presense (-0.12 ± 0.06, p<0.09) 
# R2m:0.199, R2c:0.432 




###### Average zoop size  #########
##############################################
hist(d$Ave_size)
hist(log10(d$Ave_size + 1))

uni_Ave_sz <- glm(log10(Ave_size +1) ~ scale(log10(Circum_m)), data = d)
summary(uni_Ave_sz) 

uni_Ave_sz <- glm(log10(Ave_size +1) ~ scale(log10(Sur_Area_m2)), data = d)
summary(uni_Ave_sz) 

uni_Ave_sz <- glm(log10(Ave_size +1) ~ scale(log10(Max_Depth)), data = d)
summary(uni_Ave_sz) #  

uni_Ave_sz <- glm(log10(Ave_size +1) ~ scale(NDVI), data = d)
summary(uni_Ave_sz) # ***

uni_Ave_sz <- glm(log10(Ave_size +1) ~ scale(Elevation), data = d)
summary(uni_Ave_sz) # sig

uni_Ave_sz <- glm(log10(Ave_size +1) ~ Fish, data = d)
summary(uni_Ave_sz) # **

uni_Ave_sz <- glm(log10(Ave_size +1) ~ scale(lake_order), data = d)
summary(uni_Ave_sz) # sig 

uni_Ave_sz <- glm(log10(Ave_size +1) ~ scale(sample_number), data = d)
summary(uni_Ave_sz) # sig 


# null model
zsize.null <- lmer(log10(Ave_size +1) ~ Fish  + scale(lake_order) + scale(log10(Circum_m))
                   + scale(log10(Max_Depth)) + (1| Site) +  (1| Watershed) + (1| visit_re), data = d, REML=FALSE)
summary(zsize.null)


zsize.global <- lmer(log10(Ave_size +1) ~  scale(Elevation) + Fish  + scale(lake_order) + 
                       scale(log10(Circum_m)) + scale(log10(Max_Depth)) + (1| Site) +  
                       (1| Watershed) + (1| visit_re), data = d, REML=FALSE)
summary(zsize.global)

zsize.mod <- lmer(log10(Ave_size +1) ~  scale(Elevation) + Fish +  scale(lake_order) +
                    scale(log10(Max_Depth)) + (1| Site) + (1| Watershed) + (1| visit_re), data = d, REML=FALSE)
summary(zsize.mod)

AIC(zsize.null, zsize.global, zsize.mod)
hist(residuals(zsize.null))
hist(residuals(zsize.global))
hist(residuals(zsize.mod)) # no difference between null and elevation models


r.squaredGLMM(zsize.null) 
r.squaredGLMM(zsize.global) # best model 
r.squaredGLMM(zsize.mod)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Describe:
# Zooplankton average size was positively correlated with:
# elevation (0.06 ± 0.01, p<4.66e-05)
# lake order (0.04 ± 0.01, p<0.003)
# max depth (0.03 ± 0.01, p<0.0003)
# and negatively with
# fish presense (-0.16 ± 0.02, p<2.32e-07) 
# R2m:0.853, R2c:0.946 



###### Phytoplankton genera richness #########
##############################################

hist(d$phy_gen_rich_ave)

uni_prich <- glm(phy_gen_rich_ave ~ scale(log10(Circum_m)), data = d)
summary(uni_prich) 

uni_prich <- glm(phy_gen_rich_ave ~ scale(log10(Sur_Area_m2)), data = d)
summary(uni_prich) 

uni_prich <- glm(phy_gen_rich_ave ~ scale(log10(Max_Depth)), data = d)
summary(uni_prich) #  

uni_prich <- glm(phy_gen_rich_ave ~ scale(NDVI), data = d)
summary(uni_prich) # ***

uni_prich <- glm(phy_gen_rich_ave ~ scale(Elevation), data = d)
summary(uni_prich) # sig

uni_prich <- glm(phy_gen_rich_ave ~ Fish, data = d)
summary(uni_prich) # **

uni_prich <- glm(phy_gen_rich_ave ~ scale(lake_order), data = d)
summary(uni_prich) # sig 

uni_prich <- glm(phy_gen_rich_ave ~ scale(sample_number), data = d)
summary(uni_prich) # sig 


# null model
prich.null <- lmer(phy_gen_rich_ave ~ Fish  + scale(lake_order) + scale(log10(Circum_m))
                   + scale(log10(Max_Depth)) + (1| Site) +  (1| Watershed) + (1| visit_re), data = d, REML=FALSE)
summary(prich.null)


prich.global <- lmer(phy_gen_rich_ave ~  scale(Elevation) + Fish  + scale(lake_order) + 
                       scale(log10(Circum_m)) + scale(log10(Max_Depth)) + (1| Site) +  
                       (1| Watershed) + (1| visit_re), data = d, REML=FALSE)
summary(prich.global)

prich.mod <- lmer(phy_gen_rich_ave ~  scale(Elevation) + scale(lake_order) + 
                    scale(log10(Circum_m)) + (1| Site) +  
                    (1| Watershed) + (1| visit_re), data = d, REML=FALSE)
summary(prich.mod)

AIC(prich.null, prich.global, prich.mod)
hist(residuals(prich.null))
hist(residuals(prich.global))
hist(residuals(prich.mod)) # no difference between null and elevation models


r.squaredGLMM(prich.null) 
r.squaredGLMM(prich.global) # best model 
r.squaredGLMM(prich.mod)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Describe:
# Phytoplankton genera richness was negatively correlated with:
# elevation (-2.04 ± 0.83, p<0.019)
# lake order (-2.33 ± 1.07, p<0.04)
# max depth ( -1.47 ± 0.64, p<0.03)
# and positively correlated with
# circumference (2.06 ± 0.54, p<0.0005) 
# R2m:0.853, R2c:0.946 






###### Phytoplankton biovolume ##############
##############################################
hist(d$BV_phyto_ave)
hist(log10(d$BV_phyto_ave))

uni_pbv <- glm(log10(d$BV_phyto_ave) ~ scale(log10(Circum_m)), data = d)
summary(uni_pbv) 

uni_pbv <- glm(log10(d$BV_phyto_ave) ~ scale(log10(Sur_Area_m2)), data = d)
summary(uni_pbv) 

uni_pbv <- glm(log10(d$BV_phyto_ave) ~ scale(log10(Max_Depth)), data = d)
summary(uni_pbv) #  

uni_pbv <- glm(log10(d$BV_phyto_ave) ~ scale(NDVI), data = d)
summary(uni_pbv) # ***

uni_pbv <- glm(log10(d$BV_phyto_ave) ~ scale(Elevation), data = d)
summary(uni_pbv) # sig

uni_pbv <- glm(log10(d$BV_phyto_ave) ~ Fish, data = d)
summary(uni_pbv) # **

uni_pbv <- glm(log10(d$BV_phyto_ave) ~ scale(lake_order), data = d)
summary(uni_pbv) # sig 

uni_pbv <- glm(log10(d$BV_phyto_ave) ~ scale(sample_number), data = d)
summary(uni_pbv) # sig 

# null model
pBV.null <- lmer(log10(d$BV_phyto_ave) ~ Fish  + scale(lake_order) + scale(log10(Circum_m))
                   + scale(log10(Max_Depth)) + (1| Site) +  (1| Watershed) + (1| visit_re), data = d, REML=FALSE)
summary(pBV.null)


pBV.global <- lmer(log10(d$BV_phyto_ave) ~  scale(Elevation) + Fish  + scale(lake_order) + 
                       scale(log10(Circum_m)) + scale(log10(Max_Depth)) + (1| Site) +  
                       (1| Watershed) + (1| visit_re), data = d, REML=FALSE)
summary(pBV.global)

pBV.mod <- lmer(log10(d$BV_phyto_ave) ~  scale(Elevation) + (1| Site) +  
                    (1| Watershed) + (1| visit_re), data = d, REML=FALSE)
summary(pBV.mod)

AIC(pBV.null, pBV.global, pBV.mod)
hist(residuals(pBV.null))
hist(residuals(pBV.global))
hist(residuals(pBV.mod)) # no difference between null and elevation models


r.squaredGLMM(pBV.null) 
r.squaredGLMM(pBV.global) # best model 
r.squaredGLMM(pBV.mod)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Describe:
# Phytoplankton genera richness was negatively correlated with:
# elevation (-0.25 ± 0.06, p<0.002)
#
# R2m:0.262, R2c:0.380 







#model anova's
anova(znrich.null1.a,zrich.global.a,zrich.test.1.a,zrich.null.a)
summary(zrich.global.a)

anova(zden.globala, zden.mod3, zden.null3, zden.nulla)
summary(zden.mod3)

anova(zBV.nullt.a, zBV.elev.a, zBV.null.a, zBV.global.a)
summary(zBV.nullt.a)


summary(zgrav.mod)

summary(prich.ele)

summary(zBV.elev)


data(BCI)
head(BCI)

r <- read.csv("COMP_LAKES_zoop_rich_tobin.csv", header=T)
names(r)
r$Lake
rblue <- r[c(1:3), c(3:39)]
rdia <- r[c(4:6), c(3:39)]
risa <- r[c(8:10), c(3:39)]
rjas <- r[c(11:12), c(3:39)]
rlion <- r[c(13:15), c(3:39)]
rlong <- r[c(16:18), c(3:39)]
rlost <- r[c(19:21), c(3:39)]
rmud <- r[c(22:24), c(3:39)]
rpear <- r[c(25:27), c(3:39)]
rred <- r[c(28:29), c(3:39)]
rsnow <- r[c(30:32), c(3:39)]
rupd <- r[c(33:35), c(3:39)]
ryand <- r[c(36:37), c(3:39)]





r2 <- r[, c(3:39)]
names(r2)

# all rare
S <- specnumber(r2) # observed number of species
(raremax <- min(rowSums(r2)))
Srare <- rarefy(r2, raremax)
plot(S, Srare, xlab = "Observed No. of Species", ylab = "Rarefied No. of Species")
abline(0, 1)
rarecurve(r2, step = 20, sample = raremax, col = "blue", cex = 0.6)

# blue rare
S <- specnumber(rblue) # observed number of species
(raremax <- min(rowSums(rblue)))
Srare <- rarefy(rblue, raremax)
plot(S, Srare, xlab = "Observed No. of Species", ylab = "Rarefied No. of Species")
abline(0, 1)
rarecurve(rblue, step = 20, sample = raremax, col = "blue", cex = 0.6)

# diamond rare
S <- specnumber(rdia) # observed number of species
(raremax <- min(rowSums(rdia)))
Srare <- rarefy(rdia, raremax)
plot(S, Srare, xlab = "Observed No. of Species", ylab = "Rarefied No. of Species")
abline(0, 1)
rarecurve(rdia, step = 20, sample = raremax, col = "blue", cex = 0.6)

?layout


sp1b <- specaccum(rblue)
sp1c <- specaccum(rdia)
sp1d <- specaccum(risa)
sp1e <- specaccum(rjas)
sp1f <- specaccum(rlion)
sp1g <- specaccum(rlong)
sp1h <- specaccum(rlost)
sp1j <- specaccum(rmud)
sp1k <- specaccum(rpear)
sp1l <- specaccum(rred)
sp1m <- specaccum(rsnow)
sp1n <- specaccum(rupd)
sp1o <- specaccum(ryand)



png("Species accumulation curve.png",  
    width = 3.75,
    height = 5.25,
    units = "in",
    res = 1200,
    pointsize = 4 )

nf <- layout(matrix(c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14), nrow = 7, ncol = 2))
layout.show(nf)


#plot1 blue
par(mar = c(0.3, 2.2, 0.3, 1))
par(oma=c(4,4,3,3))

mtext(text="A common x-axis label",side=1,line=0,outer=TRUE)
mtext(text="A common y-axis label",side=2,line=0,outer=TRUE)


plot(sp1b, xaxt="n")
legend("bottomright", legend=c("Blue"), cex=1.3, bty = "n")

plot(sp1c, xaxt="n")
legend("bottomright", legend=c("Diamond"), cex=1.3, bty = "n")

plot(sp1d, xaxt="n")
legend("bottomright", legend=c("Isabelle"), cex=1.3, bty = "n")

plot(sp1e,xaxt="n")
legend("bottomright", legend=c("Jasper"), cex=1.3, bty = "n")

plot(sp1f, xaxt="n")
legend("bottomright", legend=c("Lion"), cex=1.3, bty = "n")

plot(sp1g, xaxt="n")
legend("bottomright", legend=c("Long"), cex=1.3, bty = "n")

plot(sp1h, xaxt="n")
legend("bottomright", legend=c("Lost"), cex=1.3, bty = "n")
axis(1, at=c(1:3))

plot(sp1j, xaxt="n")
legend("bottomright", legend=c("Mud"), cex=1.3, bty = "n")

plot(sp1k, xaxt="n")
legend("bottomright", legend=c("Pear"), cex=1.3, bty = "n")

plot(sp1l, xaxt="n")
legend("bottomright", legend=c("Red Deer"), cex=1.3, bty = "n")

plot(sp1m, xaxt="n")
legend("bottomright", legend=c("Snowbank"), cex=1.3, bty = "n")

plot(sp1n, xaxt="n")
legend("bottomright", legend=c("Upper Diamond"), cex=1.3, bty = "n")

plot(sp1o, xlim=c(1,3), xaxt="n")
legend("bottomright", legend=c("Yankee Doodle"), cex=1.3, bty = "n")
axis(1, 1:3)


dev.off()









## Add Lomolino model using argument 'add'
plot(mod1, add = TRUE, col=2, lwd=2)
## Fit Arrhenius models to all random accumulations
mods <- fitspecaccum(sp2, "arrh")
plot(mods, col="hotpink")
boxplot(sp2, col = "yellow", border = "blue", lty=1, cex=0.3, add= TRUE)
## Use nls() methods to the list of models
sapply(mods$models, AIC)

library(lattice)
library(ggplot2)
library(dplyr)

############ plot bio models
d$Ave_size
d$prop_gravid

se <- function(x) {sd(x,na.rm=TRUE)/sqrt(length(x))}

t_rich <- aggregate(total_rich ~ Elevation, data=d2, FUN=mean) 
t_rich$SE <- aggregate(total_rich ~ Elevation, data=d2, FUN=se)[,2]

range(t_rich$total_rich)

z_rich <- aggregate(zoop_rich ~ Elevation, data=d2, FUN=mean) 
z_rich$SE <- aggregate(zoop_rich ~ Elevation, data=d2, FUN=se)[,2]
range(z_rich$Elevation)


p_rich <- aggregate(phy_gen_rich_ave ~ Elevation, data=d2, FUN=mean) 
p_rich$SE <- aggregate(phy_gen_rich_ave ~ Elevation, data=d2, FUN=se)[,2]
range(p_rich$phy_gen_rich_ave)


z_den <- aggregate(zoop_den ~ Elevation, data=d, FUN=mean) 
z_den$SE <- aggregate(zoop_den ~ Elevation, data=d, FUN=se)[,2]

z_BV <- aggregate(BV_zoop ~ Elevation, data=d, FUN=mean) 
z_BV$SE <- aggregate(BV_zoop ~ Elevation, data=d, FUN=se)[,2]

p_BV <- aggregate(log10(BV_phyto_ave +1) ~ Elevation, data=d2, FUN=mean) 
p_BV$SE <- aggregate(log10(BV_phyto_ave +1) ~ Elevation, data=d2, FUN=se)[,2]

rzp_BV <- aggregate(asin(sqrt(ratio_ZBV_PBV)) ~ Elevation, data=d, FUN=mean) 
rzp_BV$SE <- aggregate(asin(sqrt(ratio_ZBV_PBV)) ~ Elevation, data=d, FUN=se)[,2]

zgravid <- aggregate(prop_gravid ~ Elevation, data=d, FUN=mean) 
zgravid$SE <- aggregate(prop_gravid ~ Elevation, data=d, FUN=se)[,2]

ave_sz <- aggregate(Ave_size ~ Elevation, data=d2, FUN=mean) 
ave_sz$SE <- aggregate(Ave_size ~ Elevation, data=d2, FUN=se)[,2]


#### Plot ####
range(p_BV$BV_phyto_ave)
png("Biotic model results.png",  
        width = 6.4,
        height = 4.75,
        units = "in",
        res = 1200,
        pointsize = 4 )

nf <- layout(matrix(c(1, 2, 3, 4, 5, 6, 7, 8, 9), nrow = 3, ncol = 3))
layout.show(nf)

par(mar = c(2, 4.5, 2.5, 3.2) +0.1)
par(oma=c(4,4,3,3))

png("Comp_lakestrich.png",  
    width = 2,
    height = 1,
    units = "in",
    res = 1200,
    pointsize = 2 )
par(mar=c(5, 5, 4, 2.7),mgp=c(3.2,1.4,0))

plot(x= t_rich$Elevation, y=t_rich$total_rich, xlab="Elevation (m)", ylab="Total taxa richness", xlim=c(2950,3550),
     ylim=c(20, 60), pch =16, cex=2, cex.axis=2, cex.lab=2)
lablist.x<-as.vector(c(2983:3550))
axis(1, at=seq(2950, 3550, by=100), labels = F)
arrows(t_rich$Elevation, (t_rich$total_rich + t_rich$SE), 
       t_rich$Elevation, (t_rich$total_rich - t_rich$SE),
       length=0, code=3)
abline(glm(t_rich$total_rich ~ t_rich$Elevation), col="gray 40", lty=2)

dev.off()

png("Comp_lakeszrich.png",  
    width = 2,
    height = 1,
    units = "in",
    res = 1200,
    pointsize = 2 )
par(mar=c(5, 5, 4, 2.7),mgp=c(3.2,1.4,0))

plot(x= z_rich$Elevation, y=z_rich$zoop_rich, xlab="Elevation (m)", ylab="Zoo taxa richness", xlim=c(2950,3550),
     ylim=c(1, 15), pch = 16, cex=2, cex.axis=2, cex.lab=2)
lablist.x<-as.vector(c(2983:3550))
axis(1, at=seq(2950, 3550, by=100), labels = F)
arrows(z_rich$Elevation, (z_rich$zoop_rich + z_rich$SE), 
       z_rich$Elevation, (z_rich$zoop_rich - z_rich$SE),
       length=0, code=3)
abline(glm(z_rich$zoop_rich ~ z_rich$Elevation), col="gray 40", lty=2)

dev.off()


png("Comp_lakesprich.png",  
    width = 2,
    height = 1,
    units = "in",
    res = 1200,
    pointsize = 2 )
par(mar=c(5, 5, 4, 2.7),mgp=c(3.2,1.4,0))

plot(x= p_rich$Elevation, y=p_rich$phy_gen_rich_ave, xlab="Elevation (m)", ylab="Phyto taxa richness", xlim=c(2950,3550),
     ylim=c(15, 35), pch = 16,  xaxt="n", cex=2, cex.axis=2, cex.lab=2)
lablist.x<-as.vector(c(2983:3550))
axis(1, at=seq(2950, 3550, by=100), cex.axis=2)
arrows(p_rich$Elevation, (p_rich$phy_gen_rich_ave + p_rich$SE), 
       p_rich$Elevation, (p_rich$phy_gen_rich_ave - p_rich$SE),
       length=0, code=3)
abline(glm(p_rich$phy_gen_rich_ave ~ p_rich$Elevation), col="gray 40", lty=2)


dev.off()


png("Comp_lakeden.png",  
    width = 2,
    height = 1,
    units = "in",
    res = 1200,
    pointsize = 2 )
par(mar=c(5, 5, 4, 2.7),mgp=c(3.2,1.4,0))

plot(x= z_den$Elevation, y=z_den$zoop_den, xlab="Elevation (m)", ylab="Zoo density (perL)", xlim=c(2950,3550),
     ylim=c(0, 90), pch = 16,  xaxt="n", cex=2, cex.axis=2, cex.lab=2)
lablist.x<-as.vector(c(2983:3550))
axis(1, at=seq(2950, 3550, by=100), labels = F)
arrows(z_den$Elevation, (z_den$zoop_den + z_den$SE), 
       z_den$Elevation, (z_den$zoop_den - z_den$SE),
       length=0, code=3)
abline(glm(z_den$zoop_den ~ z_den$Elevation), col="gray 40", lty=2)

dev.off()

png("Comp_lakeszbv.png",  
    width = 2,
    height = 1,
    units = "in",
    res = 1200,
    pointsize = 2 )
par(mar=c(5, 5, 4, 2.7),mgp=c(3.2,1.4,0))

plot(x= z_BV$Elevation, y=z_BV$BV_zoop, xlab="Elevation (m)", ylab="Zoo bv", xlim=c(2950,3550),
     ylim=c(0, 65), pch = 16,  xaxt="n", cex=2, cex.axis=2, cex.lab=2)
lablist.x<-as.vector(c(2983:3550))
axis(1, at=seq(2950, 3550, by=100), labels = F)
arrows(z_BV$Elevation, (z_BV$BV_zoop + z_BV$SE), 
       z_BV$Elevation, (z_BV$BV_zoop - z_BV$SE),
       length=0, code=3)
abline(glm(z_BV$BV_zoop ~ z_BV$Elevation), col="gray 40", lty=2)

dev.off()

png("Comp_lakesBV.png",  
    width = 2,
    height = 1,
    units = "in",
    res = 1200,
    pointsize = 2 )
par(mar=c(5, 5, 4, 2.7),mgp=c(3.2,1.4,0))

plot(x= p_BV$Elevation, y=p_BV$`log10(BV_phyto_ave + 1)`, xlab="Elevation (m)", ylab="log10(Phyto bv + 1)", xlim=c(2950,3550),
     ylim=c(2, 5), pch = 16,  xaxt="n", cex=2, cex.axis=2, cex.lab=2)
lablist.x<-as.vector(c(2983:3550))
axis(1, at=seq(2950, 3550, by=100), cex.axis=2)
arrows(p_BV$Elevation, (p_BV$`log10(BV_phyto_ave + 1)` + p_BV$SE), 
       p_BV$Elevation, (p_BV$`log10(BV_phyto_ave + 1)` - p_BV$SE),
       length=0, code=3)
abline(glm(p_BV$`log10(BV_phyto_ave + 1)` ~ p_BV$Elevation), col="gray 40", lty=2)


dev.off()


png("Comp_lakesZ_PBV.png",  
    width = 2,
    height = 1,
    units = "in",
    res = 1200,
    pointsize = 2 )
par(mar=c(5, 5, 4, 2.7),mgp=c(3.2,1.4,0))

plot(x= rzp_BV$Elevation, y=rzp_BV$`asin(sqrt(ratio_ZBV_PBV))`, xlab="Elevation (m)", ylab="Asin(sqrt(Zoo:Phyto bv))", xlim=c(2950,3550),
     ylim=c(0, 0.25), pch = 16, cex=2, cex.axis=2, cex.lab=2)
lablist.x<-as.vector(c(2983:3550))
axis(1, at=seq(2950, 3550, by=100), cex.axis=2)
arrows(rzp_BV$Elevation, (rzp_BV$`asin(sqrt(ratio_ZBV_PBV))` + rzp_BV$SE), 
       rzp_BV$Elevation, (rzp_BV$`asin(sqrt(ratio_ZBV_PBV))` - rzp_BV$SE),
       length=0, code=3)
abline(glm(rzp_BV$`asin(sqrt(ratio_ZBV_PBV))` ~ rzp_BV$Elevation), col="gray 40", lty=2)


dev.off()

png("Comp_lakefecund.png",  
    width = 2,
    height = 1,
    units = "in",
    res = 1200,
    pointsize = 2 )
par(mar=c(5, 5, 4, 2.7),mgp=c(3.2,1.4,0))

plot(x= zgravid$Elevation, y=zgravid$prop_gravid, xlab="Elevation (m)", ylab="Asin(sqrt(proportion gravid zoo))", xlim=c(2950,3550),
     ylim=c(0, 0.52), pch = 16,  xaxt="n", cex=2, cex.axis=2, cex.lab=2)
lablist.x<-as.vector(c(2983:3550))
axis(1, at=seq(2950, 3550, by=100), cex.axis=2)
arrows(zgravid$Elevation, (zgravid$prop_gravid + zgravid$SE), 
       zgravid$Elevation, (zgravid$prop_gravid - zgravid$SE),
       length=0, code=3)
abline(glm(zgravid$prop_gravid ~ zgravid$Elevation), col="gray 40", lty=2)

dev.off()


png("Comp_lakezoopsize.png",  
    width = 2,
    height = 1,
    units = "in",
    res = 1200,
    pointsize = 2 )
par(mar=c(5, 5, 4, 2.7),mgp=c(3.2,1.4,0))
plot(x= ave_sz$Elevation, y=ave_sz$Ave_size, xlab="Elevation (m)", ylab="Zoo size (mm)", xlim=c(2950,3550),
     ylim=c(0.25, 2), pch = 16,  xaxt="n", cex=2, cex.axis=2, cex.lab=2)
lablist.x<-as.vector(c(2983:3550))
axis(1, at=seq(2950, 3550, by=100), cex.axis=2)
arrows(ave_sz$Elevation, (ave_sz$Ave_size + ave_sz$SE), 
       ave_sz$Elevation, (ave_sz$Ave_size - ave_sz$SE),
       length=0, code=3)
abline(glm(ave_sz$Ave_size ~ ave_sz$Elevation), col="gray 40", lty=2)

dev.off()

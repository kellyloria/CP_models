###################
# 2016 Comp Lakes 
# Univariate analysis of biotic charactersitics of alpine lakes along an elevation gradient
###################
# load needed packages 

library(lme4)
library(lmerTest)# for p-value
library(MuMIn) # forr squared
library(PerformanceAnalytics) 
library(nlme)


d <- read.csv("2016_Comp_Lakes_water_chem.csv", header=T)
names(d)
w <- read.csv("COMPLAKES_2016data_a.csv", header=T)
names(w)
cp <- read.csv("comp_lakes_master_data_2.csv", header=T)
names(cp)

# Subset df d for sur and hypo samples only (GLV lakes have meta too)
#surface measurements:
head(d)
d2 <- subset(d, sample_depth=="sur" | sample_depth=="hypo",
             select= Site:phyto_Bvcheck)
summary(d2$sample_depth)


#######################################
# Relevant factors for abiotic qualities of lakes and their dataframes:
#     1. Richness: Zooplankton taxa richness, phytoplankton richness, total richness
#     2. Abundance: Zoop density, phyto density, zoop biovolume, phyto biovolume, chlor-a
#     3. Population traits: Zoop size, zoop fecundity


##################
# Phyto richness #
##################
hist(cp$phy_gen_rich_ave) # fairly normal 

PHYTR.i.mod1 <- glmer(phy_gen_rich_ave ~ scale(SA_m2) + scale(max_depth) + scale(Elevation) + Fish + 
                   scale(Elevation)*Fish + scale(Elevation)*scale(max_depth) + (1|sample_num) + 
                     (1|Site), family = poisson, data = cp)

summary(PHYTR.i.mod1)
hist(residuals(PHYTR.i.mod1))
# no significant interactions 

PHYTR.mod1 <- glmer(phy_gen_rich_ave ~ scale(SA_m2) + scale(max_depth) + scale(Elevation) + Fish +
                     (1|sample_num) + (1|Site), family = poisson, data = cp)

summary(PHYTR.mod1)
hist(residuals(PHYTR.mod1))
vif(PHYTR.mod1)
# no significant drivers


##################
# Zoop richness #
##################
hist(cp$zoop_rich) # not normal 
ZOOR.i.mod1 <- glmer(zoop_rich ~ scale(SA_m2) + scale(max_depth) + scale(Elevation) + Fish + 
                        scale(Elevation)*Fish + scale(Elevation)*scale(max_depth) + (1|sample_num) + 
                        (1|Site), family = poisson, data = cp)

summary(ZOOR.i.mod1)
hist(residuals(ZOOR.i.mod1))

ZOOR.mod1 <- glmer(zoop_rich ~ scale(SA_m2) + scale(max_depth) + scale(Elevation) + Fish
                     + (1|sample_num) + (1|Site), family = poisson, data = cp)

summary(ZOOR.mod1)
hist(residuals(ZOOR.mod1))
vif(ZOOR.mod1)

##################
# Total richness #
##################
hist(cp$total_rich) # fairly normal 

TOTALR.i.mod1 <- glmer(total_rich ~ scale(SA_m2) + scale(max_depth) + scale(Elevation) + Fish + 
                       scale(Elevation)*Fish + scale(Elevation)*scale(max_depth) + (1|sample_num) + 
                       (1|Site), family = poisson, data = cp)

summary(TOTALR.i.mod1)
hist(residuals(TOTALR.i.mod1))
# no sig interactions

TOTALR.mod1 <- glmer(total_rich ~ scale(SA_m2) + scale(max_depth) + scale(Elevation) + Fish + 
                       (1|sample_num) + (1|Site), family = poisson, data = cp)

summary(TOTALR.mod1)
hist(residuals(TOTALR.mod1))
# cut total and phyto ricness-- phytoplankton 


#################
# PHYTO DENSITY #
#################
hist((d2$phyto_den))
hist(log10(d2$phyto_den)) # better

PHYTD.i.mod1 <- lmer(log10(phyto_den) ~ scale(SA_m2) + scale(max_depth) + scale(Elevation) + 
                    epi + Fish + scale(Elevation)*epi + scale(Elevation)*Fish + 
                    scale(Elevation)*scale(max_depth) + (1|sample_num) + (1|Site), data = d2)

summary(PHYTD.i.mod1)
hist(residuals(PHYTD.i.mod1))


PHYTD.mod1 <- lmer(log10(phyto_den) ~ scale(SA_m2) + scale(max_depth) + scale(Elevation) + 
                       epi + Fish + (1|sample_num) + (1|Site), data = d2)

summary(PHYTD.mod1)
hist(residuals(PHYTD.mod1))



PHYTD.mod1 <- lmer(log10(phyto_den) ~ scale(SA_m2) + scale(max_depth) + scale(Elevation) + 
                     epi + (1|sample_num) + (1|Site), data = d2)

summary(PHYTD.mod1)
hist(residuals(PHYTD.mod1))

################
# Zoop Density #
################
hist(cp$zoop_den) # not normal
hist(log10(cp$zoop_den +1)) # more normal 
hist(cp$zoop_den2)
hist(log10(cp$zoop_den2 +1))

ZOOD.i.mod1 <- lmer(log10(zoop_den2 +1) ~ scale(SA_m2) + scale(max_depth) + scale(Elevation) + Fish + 
                         scale(Elevation)*Fish + scale(Elevation)*scale(max_depth) + (1|sample_num) + 
                         (1|Site), data = cp)

summary(ZOOD.i.mod1)
hist(residuals(ZOOD.i.mod1))
#no sig interactions

ZOOD.mod1 <- lmer(log10(zoop_den2 +1) ~ scale(SA_m2) + scale(max_depth) + scale(Elevation) + Fish +
                    (1|sample_num) + (1|Site), data = cp)

summary(ZOOD.mod1)
hist(residuals(ZOOD.mod1))
vif(ZOOD.mod1)


ZOOD.mod1 <- lmer(log10(zoop_den +1) ~ scale(SA_m2) + scale(max_depth) + scale(Elevation) +
                    (1|sample_num) + (1|Site), data = cp)

summary(ZOOD.mod1)
hist(residuals(ZOOD.mod1))


#####################################
# Phyto Biovolume (total BV per mL) #
#####################################
# Calculated by totall BV / percent analyzed corrected for setting 
hist(d2$phyto_BV)# not normal
hist(log10(d2$phyto_BV)) # better

PHYTBV.i.mod1 <- lmer(log10(phyto_BV) ~ scale(SA_m2) + scale(max_depth) + scale(Elevation) + 
                       epi + Fish + scale(Elevation)*epi + scale(Elevation)*Fish + 
                       scale(Elevation)*scale(max_depth) + (1|sample_num) + (1|Site), data = d2)

summary(PHYTBV.i.mod1)
hist(residuals(PHYTBV.i.mod1))


PHYTBV.i.mod1 <- lmer(log10(phyto_BV) ~ scale(SA_m2) + scale(max_depth) + scale(Elevation) + 
                        epi + Fish + scale(Elevation)*Fish + 
                        scale(Elevation)*scale(max_depth) + (1|sample_num) + (1|Site), data = d2)

summary(PHYTBV.i.mod1)
hist(residuals(PHYTBV.i.mod1))


PHYTBV.i.mod1 <- lmer(log10(phyto_BV) ~ scale(SA_m2) + scale(max_depth) + scale(Elevation) + 
                        epi + Fish + scale(Elevation)*Fish + (1|sample_num) + (1|Site), data = d2)

summary(PHYTBV.i.mod1)
hist(residuals(PHYTBV.i.mod1))


PHYTBV.mod1 <- lmer(log10(phyto_BV) ~ scale(SA_m2) + scale(max_depth) + scale(Elevation) + 
                       epi + Fish +  (1|sample_num) + (1|Site), data = d2)

summary(PHYTBV.mod1)
hist(residuals(PHYTBV.mod1))

##################
# Zoop Biovolume #
##################
hist(cp$BV_zoop) # not normal
hist(log10(cp$BV_zoop +1)) # more normal 

ZOOBV.i.mod1 <- lmer(log10(BV_zoop +1) ~ scale(SA_m2) + scale(max_depth) + scale(Elevation) + Fish + 
                      scale(Elevation)*Fish + scale(Elevation)*scale(max_depth) + (1|sample_num) + 
                      (1|Site), data = cp)

summary(ZOOBV.i.mod1)
hist(residuals(ZOOBV.i.mod1))
#no sig interactions

ZOOBV.mod1 <- lmer(log10(BV_zoop +1) ~ scale(SA_m2) + scale(max_depth) + scale(Elevation) + Fish +
                    (1|sample_num) + (1|Site), data = cp)

summary(ZOOBV.mod1)
hist(residuals(ZOOBV.mod1))


#################
# Zoop Ave Size #
#################
hist(cp$Ave_size) # not normal
hist(log10(cp$Ave_size +1)) # maybe better 
ZOOSZ.i.mod1 <- lmer(Ave_size ~ scale(SA_m2) + scale(max_depth) + scale(Elevation) + Fish + 
                       scale(Elevation)*Fish + scale(Elevation)*scale(max_depth) + (1|sample_num) + 
                       (1|Site), data = cp)

summary(ZOOSZ.i.mod1)
hist(residuals(ZOOSZ.i.mod1))
#no sig interactions

ZOOSZ.mod1 <- lmer(Ave_size ~ scale(SA_m2) + scale(max_depth) + scale(Elevation) + Fish +
                     (1|sample_num) + (1|Site), data = cp)

summary(ZOOSZ.mod1)
hist(residuals(ZOOSZ.mod1))

#################
# Zoop Fecundity #
#################


hist(cp$prop_gravid) # not normal
hist(log10(cp$prop_gravid +1)) # better but not great
hist(asin(sqrt(cp$prop_gravid +1))) 

ZOOG.i.mod1 <- lmer(log10(prop_gravid +1)  ~ scale(SA_m2) + scale(max_depth) + scale(Elevation) + Fish + 
                       scale(Elevation)*Fish + scale(Elevation)*scale(max_depth) + (1|sample_num) + 
                       (1|Site), data = cp)

summary(ZOOG.i.mod1)
hist(residuals(ZOOG.i.mod1))
#no sig interactions

ZOOG.mod1 <- lmer(log10(prop_gravid +1)  ~ scale(SA_m2) + scale(max_depth) + scale(Elevation) + Fish +
                     (1|sample_num) + (1|Site), data = cp)

summary(ZOOG.mod1)
hist(residuals(ZOOG.mod1))



#ZOOG.i.mod1 <- glmer(cbind(#fecund, #not fecund)  ~ scale(SA_m2) + scale(max_depth) + scale(Elevation) + Fish + 
#                      scale(Elevation)*Fish + scale(Elevation)*scale(max_depth) + (1|sample_num) + 
#                      (1|Site), family=binomail data = cp)

summary(ZOOG.i.mod1)
hist(residuals(ZOOG.i.mod1))
cp$num_notgravid


ZOOGNG.i.mod1 <- glmer(cbind(num_gravid, num_notgravid)  ~ scale(SA_m2) + scale(max_depth) + scale(Elevation) + Fish + 
                      scale(Elevation)*Fish + scale(Elevation)*scale(max_depth) + (1|sample_num) + 
                      (1|Site), family =binomial, data = cp)

summary(ZOOGNG.i.mod1)
hist(residuals(ZOOGNG.i.mod1))

ZOOGNG.mod1 <- glmer(cbind(num_gravid, num_notgravid)  ~ scale(SA_m2) + scale(max_depth) + scale(Elevation) + 
                         Fish + (1|sample_num) + (1|Site), family =binomial, data = cp)

summary(ZOOGNG.mod1)
hist(residuals(ZOOGNG.i.mod1))

##################################################
##################################################
# BIOTIC FACTORS 
#     1. Richness: Zooplankton taxa richness, phytoplankton richness, total richness
#     2. Abundance: Zoop density, phyto density, zoop biovolume, phyto biovolume, chlor-a
#     3. Population traits: Zoop size, zoop fecundity

##################################################
# What aspects of biotic communities were associated with elevation?: 
#    1. Phytoplankton density: PHYTD.mod1 (trend negative for elevation and positive for surface area; sig negative for epi)
#    2. Zooplankton density: ZOOD.mod1 (sig negative for elevation and trend for fish)
#    3. Phytoplankton Biovolume: PHYTBV.i.mod1 (trend of interaction for epi*elevation, Fish*elevation, max depth*elevation)
#            (sig for epi, and max depth, and trend for elevation)
#    4. Zooplankton average size: ZOOSZ.mod1 (sig negative for fish, and positive elevation)
#    5. Zooplankton fecundity: ZOOG.mod1 (sig negative for fish, sig positive for elevation)


# Other sig drivers not related to elevation?
#    1. Total richness: TOTALR.mod1 (sig for surface area)
#    2. Zooplankton taxa richness: ZOOR.mod1 (sig for Surface area)

# not significant to anything 
#    1. Phytoplankton genera richness



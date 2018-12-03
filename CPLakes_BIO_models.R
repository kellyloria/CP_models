###################
# 2016 Comp Lakes 
# Univariate analysis of biotic charactersitics of alpine lakes along an elevation gradient
###################
# load needed packages 
library(vegan)
library(sp) #space
library(ggplot2)
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


####################
# Models to compare:
#   A. lmer(y ~ time + elevation + SA + depth + sample depth + (1|site))
#   B. lmer(y ~ time + elevation + SA + depth + (1|site))
#   C. lmer(y ~ time + elevation + SA + depth + (1|site) + (1|watershed))
#   D. lmer(y ~ time + elevation + SA + depth + sample depth + (1|site) + (1|watershed))




##################
# Phyto richness #
##################
hist(cp$phy_gen_rich_ave) # fairly normal 
names(cp)
cp$watershed

# Model B:
#   B. lmer(y ~ time + elevation + SA + depth + (1|site))
PHYT_R.B.mod1 <- lmer(phy_gen_rich_ave ~ scale(Elevation) + Fish + scale(log10(SA_m2)) 
                   + scale(log10(max_depth)) + sample_num
                   + (1 | Site), data = cp)
summary(PHYT_R.B.mod1)
hist(residuals(PHYT_R.B.mod1))
r.squaredGLMM(PHYT_R.B.mod1)

# Model D:
#   D. lmer(y ~ time + elevation + SA + depth + (1|site) + (1|watershed))
PHYT_R.D.mod1 <- lmer(phy_gen_rich_ave ~ scale(Elevation) + Fish + scale(log10(SA_m2)) 
                   + scale(log10(max_depth)) + sample_num
                   + (1 | Site) + (1| watershed), data = cp)
summary(PHYT_R.D.mod1)
hist(residuals(PHYT_R.D.mod1))
r.squaredGLMM(PHYT_R.D.mod1)

AIC(PHYT_R.B.mod1, PHYT_R.D.mod1)

##################
# Zoop richness #
##################
hist(cp$zoop_rich) # not normal 
names(cp)
cp$watershed

# Model B:
#   B. lmer(y ~ time + elevation + SA + depth + (1|site))
ZOO_R.B.mod1 <- lmer(zoop_rich ~ scale(Elevation) + Fish + scale(log10(SA_m2)) 
                       + scale(log10(max_depth)) + sample_num
                       + (1 | Site), data = cp)
summary(ZOO_R.B.mod1)
hist(residuals(ZOO_R.B.mod1))
r.squaredGLMM(ZOO_R.B.mod1)

# Model D:
#   D. lmer(y ~ time + elevation + SA + depth + (1|site) + (1|watershed))
ZOOP_R.D.mod1 <- lmer(zoop_rich ~ scale(Elevation) + Fish + scale(log10(SA_m2)) 
                       + scale(log10(max_depth)) + sample_num
                       + (1 | Site) + (1| watershed), data = cp)
summary(ZOOP_R.D.mod1)
hist(residuals(ZOOP_R.D.mod1))
r.squaredGLMM(ZOOP_R.D.mod1)

AIC(ZOO_R.B.mod1, ZOOP_R.D.mod1)

##################
# Total richness #
##################
hist(cp$total_rich) # fairly normal 
names(cp)
cp$watershed

# Model B:
#   B. lmer(y ~ time + elevation + SA + depth + (1|site))
TOTAL_R.B.mod1 <- lmer(total_rich ~ scale(Elevation) + Fish + scale(log10(SA_m2)) 
                       + scale(log10(max_depth)) + sample_num
                       + (1 | Site), data = cp)
summary(TOTAL_R.B.mod1)
hist(residuals(TOTAL_R.B.mod1))
r.squaredGLMM(TOTAL_R.B.mod1)

# Model D:
#   D. lmer(y ~ time + elevation + SA + depth + (1|site) + (1|watershed))
TOTAL_R.D.mod1 <- lmer(total_rich ~ scale(Elevation) + Fish + scale(log10(SA_m2)) 
                       + scale(log10(max_depth)) + sample_num
                       + (1 | Site) + (1| watershed), data = cp)
summary(TOTAL_R.D.mod1)
hist(residuals(TOTAL_R.D.mod1))
r.squaredGLMM(TOTAL_R.D.mod1)

AIC(TOTAL_R.B.mod1, TOTAL_R.D.mod1)

###########
# CHLOR-A #
###########
hist((d2$chla))
hist(log10(d2$chla +1)) # better
# Model A:
#   A. lmer(y ~ time + elevation + SA + depth + sample depth + (1|site))
CHLA.A.mod1 <- lmer(log10(chla +1) ~ scale(Elevation) + Fish + scale(log10(SA_m2)) 
                    + scale(log10(max_depth)) + sample_depth +  sample_num
                    + (1 | Site), data = d2)
summary(CHLA.A.mod1)
hist(residuals(CHLA.A.mod1))
r.squaredGLMM(CHLA.A.mod1)

# Model B:
#   B. lmer(y ~ time + elevation + SA + depth + (1|site))
CHLA.B.mod1 <- lmer(log10(chla +1)  ~ scale(Elevation) + Fish + scale(log10(SA_m2)) 
                    + scale(log10(max_depth)) + sample_num
                    + (1 | Site), data = d2)
summary(CHLA.B.mod1)
hist(residuals(CHLA.B.mod1))
r.squaredGLMM(CHLA.B.mod1)

# Model C:
#   C. lmer(y ~ time + elevation + SA + depth + (1|site) + (1|watershed))
CHLA.C.mod1 <- lmer(log10(chla +1) ~ scale(Elevation) + Fish + scale(log10(SA_m2)) 
                    + scale(log10(max_depth)) + sample_depth +  sample_num
                    + (1 | Site) + (1| watershed), data = d2)

summary(CHLA.C.mod1)
hist(residuals(CHLA.C.mod1))
r.squaredGLMM(CHLA.C.mod1)

# Model D:
#   D. lmer(y ~ time + elevation + SA + depth + sample depth + (1|site) + (1|watershed))
CHLA.D.mod1 <- lmer(log10(chla +1) ~ scale(Elevation) + Fish + scale(log10(SA_m2)) 
                    + scale(log10(max_depth)) + sample_num
                    + (1 | Site) + (1| watershed), data = d2)

summary(CHLA.D.mod1)
hist(residuals(CHLA.D.mod1))
r.squaredGLMM(CHLA.D.mod1)

AIC(CHLA.A.mod1, CHLA.B.mod1, CHLA.C.mod1, CHLA.D.mod1)

#################
# PHYTO DENSITY #
#################
hist((d2$phyto_den))
hist(log10(d2$phyto_den)) # better
# Model A:
#   A. lmer(y ~ time + elevation + SA + depth + sample depth + (1|site))
PHYTD.A.mod1 <- lmer(log10(phyto_den) ~ scale(Elevation) + Fish + scale(log10(SA_m2)) 
                     + scale(log10(max_depth)) + sample_depth +  sample_num
                     + (1 | Site), data = d2)
summary(PHYTD.A.mod1)
hist(residuals(PHYTD.A.mod1))
r.squaredGLMM(PHYTD.A.mod1)

# Model B:
#   B. lmer(y ~ time + elevation + SA + depth + (1|site))
PHYTD.B.mod1 <- lmer(log10(phyto_den)  ~ scale(Elevation) + Fish + scale(log10(SA_m2)) 
                     + scale(log10(max_depth)) + sample_num
                     + (1 | Site), data = d2)
summary(PHYTD.B.mod1)
hist(residuals(PHYTD.B.mod1))
r.squaredGLMM(PHYTD.B.mod1)

# Model C:
#   C. lmer(y ~ time + elevation + SA + depth + (1|site) + (1|watershed))
PHYTD.C.mod1 <- lmer(log10(phyto_den) ~ scale(Elevation) + Fish + scale(log10(SA_m2)) 
                     + scale(log10(max_depth)) + sample_depth +  sample_num
                     + (1 | Site) + (1| watershed), data = d2)

summary(PHYTD.C.mod1)
hist(residuals(PHYTD.C.mod1))
r.squaredGLMM(PHYTD.C.mod1)

# Model D:
#   D. lmer(y ~ time + elevation + SA + depth + sample depth + (1|site) + (1|watershed))
PHYTD.D.mod1 <- lmer(log10(phyto_den) ~ scale(Elevation) + Fish + scale(log10(SA_m2)) 
                     + scale(log10(max_depth)) + sample_num
                     + (1 | Site) + (1| watershed), data = d2)

summary(PHYTD.D.mod1)
hist(residuals(PHYTD.D.mod1))
r.squaredGLMM(PHYTD.D.mod1)

AIC(PHYTD.A.mod1, PHYTD.B.mod1, PHYTD.C.mod1, PHYTD.D.mod1)

################
# Zoop Density #
################
hist(cp$zoop_den) # not normal
hist(log10(cp$zoop_den +1)) # more normal 
names(cp)

# Model B:
#   B. lmer(y ~ time + elevation + SA + depth + (1|site))
ZOO_D.B.mod1 <- lmer((log10(zoop_den +1)) ~ scale(Elevation) + Fish + scale(log10(SA_m2)) 
                       + scale(log10(max_depth)) + sample_num
                       + (1 | Site), data = cp)
summary(ZOO_D.B.mod1)
hist(residuals(ZOO_D.B.mod1))
r.squaredGLMM(ZOO_D.B.mod1)

# Model D:
#   D. lmer(y ~ time + elevation + SA + depth + (1|site) + (1|watershed))
ZOO_D.D.mod1 <- lmer((log10(zoop_den +1)) ~ scale(Elevation) + Fish + scale(log10(SA_m2)) 
                       + scale(log10(max_depth)) + sample_num
                       + (1 | Site) + (1| watershed), data = cp)
summary(ZOO_D.D.mod1)
hist(residuals(ZOO_D.D.mod1))
r.squaredGLMM(ZOO_D.D.mod1)

AIC(ZOO_D.B.mod1, ZOO_D.D.mod1)

#####################################
# Phyto Biovolume (total BV per mL) #
#####################################
# Calculated by totall BV / percent analyzed corrected for setting 
hist(cp$BV_phyto_ave) #not normal 
hist(d2$phyto_BV)# not normal
hist(log10(d2$phyto_BV)) # better
summary(d2$phyto_BV)

# Model A:
#   A. lmer(y ~ time + elevation + SA + depth + sample depth + (1|site))
PHYT_BV.A.mod1 <- lmer((log10(phyto_BV)) ~ scale(Elevation) + Fish + scale(log10(SA_m2)) 
                    + scale(log10(max_depth)) + sample_depth +  sample_num
                    + (1 | Site), data = d2)
summary(PHYT_BV.A.mod1)
hist(residuals(PHYT_BV.A.mod1))
r.squaredGLMM(PHYT_BV.A.mod1)

# Model B:
#   B. lmer(y ~ time + elevation + SA + depth + (1|site))
PHYT_BV.B.mod1 <- lmer((log10(phyto_BV))  ~ scale(Elevation) + Fish + scale(log10(SA_m2)) 
                    + scale(log10(max_depth)) + sample_num
                    + (1 | Site), data = d2)
summary(PHYT_BV.B.mod1)
hist(residuals(PHYT_BV.B.mod1))
r.squaredGLMM(PHYT_BV.B.mod1)

# Model C:
#   C. lmer(y ~ time + elevation + SA + depth + (1|site) + (1|watershed))
PHYT_BV.C.mod1 <- lmer((log10(phyto_BV)) ~ scale(Elevation) + Fish + scale(log10(SA_m2)) 
                    + scale(log10(max_depth)) + sample_depth +  sample_num
                    + (1 | Site) + (1| watershed), data = d2)

summary(PHYT_BV.C.mod1)
hist(residuals(PHYT_BV.C.mod1))
r.squaredGLMM(PHYT_BV.C.mod1)

# Model D:
#   D. lmer(y ~ time + elevation + SA + depth + sample depth + (1|site) + (1|watershed))
PHYT_BV.D.mod1 <- lmer((log10(phyto_BV)) ~ scale(Elevation) + Fish + scale(log10(SA_m2)) 
                    + scale(log10(max_depth)) + sample_num
                    + (1 | Site) + (1| watershed), data = d2)

summary(PHYT_BV.D.mod1)
hist(residuals(PHYT_BV.D.mod1))
r.squaredGLMM(PHYT_BV.D.mod1)

AIC(PHYT_BV.A.mod1, PHYT_BV.B.mod1, PHYT_BV.C.mod1, PHYT_BV.D.mod1)

##################
# Zoop Biovolume #
##################
hist(cp$BV_zoop) # not normal
hist(log10(cp$BV_zoop +1)) # more normal 
names(cp)

# Model B:
#   B. lmer(y ~ time + elevation + SA + depth + (1|site))
ZOO_BV.B.mod1 <- lmer((log10(BV_zoop +1)) ~ scale(Elevation) + Fish + scale(log10(SA_m2)) 
                     + scale(log10(max_depth)) + sample_num
                     + (1 | Site), data = cp)
summary(ZOO_BV.B.mod1)
hist(residuals(ZOO_BV.B.mod1))
r.squaredGLMM(ZOO_BV.B.mod1)

# Model D:
#   D. lmer(y ~ time + elevation + SA + depth + (1|site) + (1|watershed))
ZOO_BV.D.mod1 <- lmer((log10(BV_zoop +1)) ~ scale(Elevation) + Fish + scale(log10(SA_m2)) 
                     + scale(log10(max_depth)) + sample_num
                     + (1 | Site) + (1| watershed), data = cp)
summary(ZOO_BV.D.mod1)
hist(residuals(ZOO_BV.D.mod1))
r.squaredGLMM(ZOO_BV.D.mod1)

AIC(ZOO_BV.B.mod1, ZOO_BV.D.mod1)
# models not significant 

#################
# Zoop Ave Size #
#################
hist(cp$Ave_size) # not normal
hist(log10(cp$Ave_size +1)) # maybe better 
names(cp)

# Model B:
#   B. lmer(y ~ time + elevation + SA + depth + (1|site))
ZOO_SZ.B.mod1 <- lmer((log10(Ave_size +1)) ~ scale(Elevation) + Fish + scale(log10(SA_m2)) 
                     + scale(log10(max_depth)) + sample_num
                     + (1 | Site), data = cp)
summary(ZOO_SZ.B.mod1)
hist(residuals(ZOO_SZ.B.mod1))
r.squaredGLMM(ZOO_SZ.B.mod1)

# Model D:
#   D. lmer(y ~ time + elevation + SA + depth + (1|site) + (1|watershed))
ZOO_SZ.D.mod1 <- lmer((log10(Ave_size +1)) ~ scale(Elevation) + Fish + scale(log10(SA_m2)) 
                     + scale(log10(max_depth)) + sample_num
                     + (1 | Site) + (1| watershed), data = cp)
summary(ZOO_SZ.D.mod1)
hist(residuals(ZOO_SZ.D.mod1))
r.squaredGLMM(ZOO_SZ.D.mod1)

AIC(ZOO_SZ.B.mod1, ZOO_SZ.D.mod1)

#################
# Zoop Fecundity #
#################
hist(cp$prop_gravid) # not normal
hist(log10(cp$prop_gravid +1)) # better but not great
hist(asin(sqrt(cp$prop_gravid +1))) 
names(cp)

# Model B:
#   B. lmer(y ~ time + elevation + SA + depth + (1|site))
ZOO_G.B.mod1 <- lmer((log10(prop_gravid +1)) ~ scale(Elevation) + Fish + scale(log10(SA_m2)) 
                      + scale(log10(max_depth)) + sample_num
                      + (1 | Site), data = cp)
summary(ZOO_G.B.mod1)
hist(residuals(ZOO_G.B.mod1))
r.squaredGLMM(ZOO_G.B.mod1)

# Model D:
#   D. lmer(y ~ time + elevation + SA + depth + (1|site) + (1|watershed))
ZOO_G.D.mod1 <- lmer((log10(prop_gravid +1)) ~ scale(Elevation) + Fish + scale(log10(SA_m2)) 
                      + scale(log10(max_depth)) + sample_num
                      + (1 | Site) + (1| watershed), data = cp)
summary(ZOO_G.D.mod1)
hist(residuals(ZOO_G.D.mod1))
r.squaredGLMM(ZOO_G.D.mod1)

AIC(ZOO_G.B.mod1, ZOO_G.D.mod1)

##################################################
##################################################
# BIOTIC FACTORS 
#     1. Richness: Zooplankton taxa richness, phytoplankton richness, total richness
#     2. Abundance: Zoop density, phyto density, zoop biovolume, phyto biovolume, chlor-a
#     3. Population traits: Zoop size, zoop fecundity

##################################################
# What aspects of biotic communities were associated with elevation?: 
#    1. Zooplankton density: ZOO_D.B.mod1 (sig negative for fish)
#    2. Zooplankton average size: ZOO_SZ.B.mod1 (sig negative for fish, sig positive for sample number)
#    3. Zooplankton fecundity: ZOO_G.B.mod1 (sig negative for fish, sig positive for sample number)


# Other sig drivers not related to elevation?
#    1. Phytoplankton genera richness: PHYT_R.B.mod1 (sig for surface area)
#    2. Zooplankton taxa richness: ZOO_R.B.mod1 (sig for sample number so time)
#    3. Phytoplankton Biovolume: PHYT_BV.C.mod1 (sig negative for sample location at surface)



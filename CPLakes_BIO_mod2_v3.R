###################
# 2016 Comp Lakes 
# Univariate analysis of biotic charactersitics of alpine lakes along an elevation gradient
###################
# Post revisions # 

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

summary(cp$Site)

sz <- read.csv("CL_zoo_size_piet.csv", header=T)
names(sz)

# Subset df d for sur and hypo samples only (GLV lakes have meta too)
#surface measurements:
head(d)
d2 <- subset(d, sample_depth=="sur" | sample_depth=="hypo",
             select= Site:phyto_Bvcheck)
summary(d2$sample_depth)


summary(d2$Site)
d_v1 <- subset(d2, Site=="Albion" | Site=="GL4" | Site=="Lost Lake" 
               | Site=="Blue Lake" | Site=="Isabelle Lake" | Site=="Snowbank Lake" 
               | Site=="Diamond Lake" | Site=="Mud Lake" | Site=="Jasper Lake" 
               | Site=="Upper Diamond Lake "| Site=="Forest Lake"| Site=="Lion Lake 2"
               | Site=="Pear Reservoir" | Site=="Yankee Doodle " | Site=="GL1" 
               | Site=="Long Lake" | Site=="Red Deer Lake" | Site=="GL1"
               | Site=="Mitchell Lake" | Site=="Red Rock Lake",
               select= Site:phyto_Bvcheck)
summary(d_v1$Site)


w_v1 <- subset(w2, Site=="Albion" | Site=="GL4" | Site=="Lost Lake" 
               | Site=="Blue Lake" | Site=="Isabelle Lake" | Site=="Snowbank Lake" 
               | Site=="Diamond Lake" | Site=="Mud Lake" | Site=="Jasper Lake" 
               | Site=="Upper Diamond Lake "| Site=="Forest Lake"| Site=="Lion Lake 2"
               | Site=="Pear Reservoir" | Site=="Yankee Doodle " | Site=="GL1" 
               | Site=="Long Lake" | Site=="Red Deer Lake" | Site=="GL1"
               | Site=="Mitchell Lake" | Site=="Red Rock Lake", select= Date:Elevation)
summary(w_v1$Site)


#######################################
# Relevant factors for abiotic qualities of lakes and their dataframes:
#     1. Richness: Zooplankton taxa richness, phytoplankton richness, total richness
#     2. Abundance: Zoop density, phyto density, zoop biovolume, phyto biovolume, chlor-a
#     3. Population traits: Zoop size, zoop fecundity


##################
# Phyto richness #
##################
hist(cp$phy_gen_rich_ave)
PHYTR.i.mod1 <- glmer(phy_gen_rich_ave ~ scale(SA_m2) + scale(max_depth) + scale(Elevation) + Fish + 
                        scale(Elevation)*Fish + scale(Elevation)*scale(max_depth) + (1|sample_num) + 
                        (1|Site), family = poisson, data = cp)

summary(PHYTR.i.mod1)
hist(residuals(PHYTR.i.mod1))

PHYTR.i.mod2 <- glmer(phy_gen_rich_ave ~ scale(max_depth) + scale(Elevation) + Fish + 
                        scale(Elevation)*Fish + (1|sample_num) + 
                        (1|Site), family = poisson, data = cp)

summary(PHYTR.i.mod2)
hist(residuals(PHYTR.i.mod2))


PHYTR.mod1 <- glmer(phy_gen_rich_ave ~ scale(SA_m2) + scale(max_depth) + scale(Elevation) + Fish 
                      + (1|sample_num) + (1|Site), family = poisson, data = cp)

summary(PHYTR.mod1)
hist(residuals(PHYTR.mod1))
r.squaredGLMM(PHYTR.mod1)

#... not clear

##################
# Total richness #
##################
# zooplankton sp. richness
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
r.squaredGLMM(ZOOR.mod1)

#################
# PHYTO DENSITY #
#################
hist((d_v1$phyto_den))
hist(log10(d_v1$phyto_den)) # better

PHYTD.i.mod1 <- lmer(log10(phyto_den) ~ scale(SA_m2) + scale(max_depth) + scale(Elevation) + 
                       epi + Fish + scale(Elevation)*epi + scale(Elevation)*Fish + 
                       scale(Elevation)*scale(max_depth) + (1|sample_num) + (1|Site), data = d_v1)

summary(PHYTD.i.mod1)
hist(residuals(PHYTD.i.mod1))


PHYTD.mod1 <- lmer(log10(phyto_den) ~ scale(SA_m2) + scale(max_depth) + scale(Elevation) + 
                     epi + Fish + (1|sample_num) + (1|Site), data = d_v1)

summary(PHYTD.mod1)
hist(residuals(PHYTD.mod1)) # significant
r.squaredGLMM(PHYTD.mod1)


PHYTD.mod2 <- lmer(log10(phyto_den) ~ scale(SA_m2) + scale(max_depth) + scale(Elevation) + 
                     epi + Fish + scale(sample_num) + (1|Site), data = d_v1)

summary(PHYTD.mod2)
hist(residuals(PHYTD.mod2)) # significant
r.squaredGLMM(PHYTD.mod2)

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
r.squaredGLMM(ZOOD.mod1)


ZOOD.mod1 <- lmer(log10(zoop_den2 +1) ~ scale(SA_m2) + scale(max_depth) + scale(Elevation) + Fish +
                    scale(sample_num) + (1|Site), data = cp)

summary(ZOOD.mod1)
hist(residuals(ZOOD.mod1))

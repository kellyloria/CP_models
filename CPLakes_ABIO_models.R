###################
# 2016 Comp Lakes 
# Univariate analysis of abiotic charactersitics of alpine lakes along an elevation gradient
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
             select= Site:phyto_den)
summary(d2$sample_depth)

summary(d2$Site)
d3 <- subset(d2, Site=="Albion" | Site=="GL4" | Site=="Lost Lake" 
             | Site=="Blue Lake" | Site=="Isabelle Lake" | Site=="Snowbank Lake" 
             | Site=="Diamond Lake" | Site=="Mud Lake" | Site=="Jasper Lake" 
             | Site=="Upper Diamond Lake"| Site=="Forest Lake"| Site=="Lion Lake 2"
             | Site=="Pear Reservoir" | Site=="Yankee Doodle" | Site=="GL1" 
             | Site=="Long Lake" | Site=="Red Deer Lake" | Site=="GL1", 
             select= Site:phyto_den)
summary(d3$Site)

head(w)
w2 <- subset(w, sample_depth=="sur" | sample_depth=="hypo",
             select= Date:Elevation)
summary(w2$sample_depth)
w2$Site
w3 <- subset(w, sample_depth=="sur",
             select= Date:Elevation)
summary(w3$sample_depth)

w4 <- subset(w2, Site=="Albion" | Site=="GL4" | Site=="Lost Lake" 
             | Site=="Blue Lake" | Site=="Isabelle Lake" | Site=="Snowbank Lake" 
             | Site=="Diamond Lake" | Site=="Mud Lake" | Site=="Jasper Lake" 
             | Site=="Upper Diamond Lake"| Site=="Forest Lake"| Site=="Lion Lake 2"
             | Site=="Pear Reservoir" | Site=="Yankee Doodle" | Site=="GL1" 
             | Site=="Long Lake" | Site=="Red Deer Lake" | Site=="GL1", 
             select= Date:Elevation)
summary(w4$Site)


names(cp)
cp2 <- subset(cp, Site=="Albion" | Site=="GL4" | Site=="Lost Lake" 
             | Site=="Blue Lake" | Site=="Isabelle Lake" | Site=="Snowbank Lake" 
             | Site=="Diamond Lake" | Site=="Mud Lake" | Site=="Jasper Lake" 
             | Site=="Upper Diamond Lake"| Site=="Forest Lake"| Site=="Lion Lake 2"
             | Site=="Pear Reservoir" | Site=="Yankee Doodle" | Site=="GL1" 
             | Site=="Long Lake" | Site=="Red Deer Lake" | Site=="GL1", 
             select= Site:Sort)
summary(cp2$Site)

## Physical variables for lake "morphology" and elevation
#     - Size: surface area and max depth 

hist(cp$Sur_Area_m2) # not normal 
hist(cp$Max_Depth) # not normal 

# Are lake surface area and depth associated with elevation?
sur_area.mod <- lmer(log10(Sur_Area_m2) ~ scale(Elevation) + (1 | Site), data = cp)
summary(sur_area.mod) # no pattern
hist(residuals(sur_area.mod))
plot(cp$Elevation, log10(cp$Sur_Area_m2))
# no 

depth.mod <- lmer(log10(Max_Depth) ~ scale(Elevation) + (1 | Site), data = cp)
summary(depth.mod) # no pattern 
hist(residuals(depth.mod))
plot(cp$Elevation, log10(cp$Max_Depth))
# no

# Are max depth and surface area correlated to each other?
my_cor <- cp[, c(7, 9)]
chart.Correlation(my_cor, histogram=TRUE, pch=19) 
# corelation: r = 0.35


#######################################
# Relevant factors for abiotic qualities of lakes and their dataframes: 
# DF= w: Temp, chlor-a, conductivity, secchi, light at surface
# DF= d: doc, TDN, TDP, NO3, PO4, FI, SUVA
# DF= cp: rho change


####################
# Models to compare:
#   A. lmer(y ~ time + elevation + SA + depth + sample depth + (1|site))
#   B. lmer(y ~ time + elevation + SA + depth + (1|site))
#   C. lmer(y ~ time + elevation + SA + depth + (1|site) + (1|watershed))
#   D. lmer(y ~ time + elevation + SA + depth + sample depth + (1|site) + (1|watershed))

#########################################################################################

#######
# DOC #
#######
hist((d2$DOC_mg_L))
hist(log10(d2$DOC_mg_L +1)) #better

# Model A:
#   A. lmer(y ~ time + elevation + SA + depth + sample depth + (1|site))
doc.A.mod1 <- lmer((log10(DOC_mg_L +1)) ~ scale(Elevation) + Fish + scale(log10(SA_m2)) 
                   + scale(log10(max_depth)) + sample_depth +  sample_num
                   + (1 | Site), data = d2)
summary(doc.A.mod1)
hist(residuals(doc.A.mod1))
r.squaredGLMM(doc.A.mod1)

# Model B:
#   B. lmer(y ~ time + elevation + SA + depth + (1|site))
doc.B.mod1 <- lmer((log10(DOC_mg_L +1)) ~ scale(Elevation) + Fish + scale(log10(SA_m2)) 
                   + scale(log10(max_depth)) + sample_num
                   + (1 | Site), data = d2)
summary(doc.B.mod1)
hist(residuals(doc.B.mod1))
r.squaredGLMM(doc.B.mod1)

# Model C:
#   C. lmer(y ~ time + elevation + SA + depth + (1|site) + (1|watershed))
doc.C.mod1 <- lmer((log10(DOC_mg_L +1)) ~ scale(Elevation) + Fish + scale(log10(SA_m2)) 
                   + scale(log10(max_depth)) + sample_depth +  sample_num
                   + (1 | Site) + (1| watershed), data = d2)
summary(doc.C.mod1)
hist(residuals(doc.C.mod1))
r.squaredGLMM(doc.C.mod1)

# Model D:
#   D. lmer(y ~ time + elevation + SA + depth + sample depth + (1|site) + (1|watershed))
doc.D.mod1 <- lmer((log10(DOC_mg_L +1)) ~ scale(Elevation) + Fish + scale(log10(SA_m2)) 
                   + scale(log10(max_depth)) + sample_num 
                   + (1 | Site) + (1| watershed), data = d2)
summary(doc.D.mod1)
hist(residuals(doc.D.mod1))
r.squaredGLMM(doc.D.mod1)

AIC(doc.A.mod1, doc.B.mod1, doc.C.mod1, doc.D.mod1) 
#D = best model




#######
# TDN #
#######
hist((d2$TDN_mg_L))

# Model A:
#   A. lmer(y ~ time + elevation + SA + depth + sample depth + (1|site))
TDN.A.mod1 <- lmer(TDN_mg_L ~ scale(Elevation) + Fish + scale(log10(SA_m2)) 
                   + scale(log10(max_depth)) + sample_depth +  sample_num
                   + (1 | Site), data = d2)
summary(TDN.A.mod1)
hist(residuals(TDN.A.mod1))
r.squaredGLMM(TDN.A.mod1)

# Model B:
#   B. lmer(y ~ time + elevation + SA + depth + (1|site))
TDN.B.mod1 <- lmer(TDN_mg_L ~ scale(Elevation) + Fish + scale(log10(SA_m2)) 
                   + scale(log10(max_depth)) + sample_num
                   + (1 | Site), data = d2)
summary(TDN.B.mod1)
hist(residuals(TDN.B.mod1))
r.squaredGLMM(TDN.B.mod1)

# Model C:
#   C. lmer(y ~ time + elevation + SA + depth + (1|site) + (1|watershed))
TDN.C.mod1 <- lmer(TDN_mg_L ~ scale(Elevation) + Fish + scale(log10(SA_m2)) 
                   + scale(log10(max_depth)) + sample_depth +  sample_num
                   + (1 | Site) + (1| watershed), data = d2)

summary(TDN.C.mod1)
hist(residuals(TDN.C.mod1))
r.squaredGLMM(TDN.C.mod1)

# Model D:
#   D. lmer(y ~ time + elevation + SA + depth + sample depth + (1|site) + (1|watershed))
TDN.D.mod1 <- lmer(TDN_mg_L ~ scale(Elevation) + Fish + scale(log10(SA_m2)) 
                   + scale(log10(max_depth)) + sample_num
                   + (1 | Site) + (1| watershed), data = d2)

summary(TDN.D.mod1)
hist(residuals(TDN.D.mod1))
r.squaredGLMM(TDN.D.mod1)

AIC(TDN.A.mod1, TDN.B.mod1, TDN.C.mod1, TDN.D.mod1) 
# B = best model 




#######
# TDP #
#######
hist((d2$TDP_mg_L))
hist(log10(d2$TDP_mg_L +1)) # not better


# Model A:
#   A. lmer(y ~ time + elevation + SA + depth + sample depth + (1|site))
TDP.A.mod1 <- lmer(TDP_mg_L ~ scale(Elevation) + Fish + scale(log10(SA_m2)) 
                   + scale(log10(max_depth)) + sample_depth +  sample_num
                   + (1 | Site), data = d2)
summary(TDP.A.mod1)
hist(residuals(TDP.A.mod1))
r.squaredGLMM(TDP.A.mod1)

# Model B:
#   B. lmer(y ~ time + elevation + SA + depth + (1|site))
TDP.B.mod1 <- lmer(TDP_mg_L ~ scale(Elevation) + Fish + scale(log10(SA_m2)) 
                   + scale(log10(max_depth)) + sample_num
                   + (1 | Site), data = d2)
summary(TDP.B.mod1)
hist(residuals(TDP.B.mod1))
r.squaredGLMM(TDP.B.mod1)

# Model C:
#   C. lmer(y ~ time + elevation + SA + depth + (1|site) + (1|watershed))
TDP.C.mod1 <- lmer(TDP_mg_L ~ scale(Elevation) + Fish + scale(log10(SA_m2)) 
                   + scale(log10(max_depth)) + sample_depth +  sample_num
                   + (1 | Site) + (1| watershed), data = d2)

summary(TDP.C.mod1)
hist(residuals(TDP.C.mod1))
r.squaredGLMM(TDP.C.mod1)

# Model D:
#   D. lmer(y ~ time + elevation + SA + depth + sample depth + (1|site) + (1|watershed))
TDP.D.mod1 <- lmer(TDP_mg_L ~ scale(Elevation) + Fish + scale(log10(SA_m2)) 
                   + scale(log10(max_depth)) + sample_num
                   + (1 | Site) + (1| watershed), data = d2)

summary(TDP.D.mod1)
hist(residuals(TDP.D.mod1))
r.squaredGLMM(TDP.D.mod1)

AIC(TDP.A.mod1, TDP.B.mod1, TDP.C.mod1, TDP.D.mod1)
# b = best model but no obvious drivers 


#######
# NO3 #
#######
hist((d2$NO3_mg_L))
hist(log10(d2$NO3_mg_L +1)) # not better
hist(asin(sqrt(d2$NO3_mg_L)))


# Model A:
#   A. lmer(y ~ time + elevation + SA + depth + sample depth + (1|site))
NO3.A.mod1 <- lmer(NO3_mg_L ~ scale(Elevation) + Fish + scale(log10(SA_m2)) 
                   + scale(log10(max_depth)) + sample_depth +  sample_num
                   + (1 | Site), data = d2)
summary(NO3.A.mod1)
hist(residuals(NO3.A.mod1))
r.squaredGLMM(NO3.A.mod1)

# Model B:
#   B. lmer(y ~ time + elevation + SA + depth + (1|site))
NO3.B.mod1 <- lmer(NO3_mg_L ~ scale(Elevation) + Fish + scale(log10(SA_m2)) 
                   + scale(log10(max_depth)) + sample_num
                   + (1 | Site), data = d2)
summary(NO3.B.mod1)
hist(residuals(NO3.B.mod1))
r.squaredGLMM(NO3.B.mod1)

# Model C:
#   C. lmer(y ~ time + elevation + SA + depth + (1|site) + (1|watershed))
NO3.C.mod1 <- lmer(NO3_mg_L ~ scale(Elevation) + Fish + scale(log10(SA_m2)) 
                   + scale(log10(max_depth)) + sample_depth +  sample_num
                   + (1 | Site) + (1| watershed), data = d2)

summary(NO3.C.mod1)
hist(residuals(NO3.C.mod1))
r.squaredGLMM(NO3.C.mod1)

# Model D:
#   D. lmer(y ~ time + elevation + SA + depth + sample depth + (1|site) + (1|watershed))
NO3.D.mod1 <- lmer(NO3_mg_L ~ scale(Elevation) + Fish + scale(log10(SA_m2)) 
                   + scale(log10(max_depth)) + sample_num
                   + (1 | Site) + (1| watershed), data = d2)

summary(NO3.D.mod1)
hist(residuals(NO3.D.mod1))
r.squaredGLMM(NO3.D.mod1)

AIC(NO3.A.mod1, NO3.B.mod1, NO3.C.mod1, NO3.D.mod1)
# model D or C is best



#######
# PO4 #
#######
hist((d2$PO4_mg_L)) #normal

# Model A:
#   A. lmer(y ~ time + elevation + SA + depth + sample depth + (1|site))
PO4.A.mod1 <- lmer(PO4_mg_L ~ scale(Elevation) + Fish + scale(log10(SA_m2)) 
                   + scale(log10(max_depth)) + sample_depth +  sample_num
                   + (1 | Site), data = d2)
summary(PO4.A.mod1) # elevation, fish, sample depth and sample number
hist(residuals(PO4.A.mod1))
r.squaredGLMM(PO4.A.mod1)

# Model B:
#   B. lmer(y ~ time + elevation + SA + depth + (1|site))
PO4.B.mod1 <- lmer(PO4_mg_L ~ scale(Elevation) + Fish + scale(log10(SA_m2)) 
                   + scale(log10(max_depth)) + sample_num
                   + (1 | Site), data = d2)
summary(PO4.B.mod1)
hist(residuals(PO4.B.mod1))
r.squaredGLMM(PO4.B.mod1)

# Model C:
#   C. lmer(y ~ time + elevation + SA + depth + (1|site) + (1|watershed))
PO4.C.mod1 <- lmer(PO4_mg_L ~ scale(Elevation) + Fish + scale(log10(SA_m2)) 
                   + scale(log10(max_depth)) + sample_depth +  sample_num
                   + (1 | Site) + (1| watershed), data = d2)

summary(PO4.C.mod1)
hist(residuals(PO4.C.mod1))
r.squaredGLMM(PO4.C.mod1)

# Model D:
#   D. lmer(y ~ time + elevation + SA + depth + sample depth + (1|site) + (1|watershed))
PO4.D.mod1 <- lmer(PO4_mg_L ~ scale(Elevation) + Fish + scale(log10(SA_m2)) 
                   + scale(log10(max_depth)) + sample_num
                   + (1 | Site) + (1| watershed), data = d2)

summary(PO4.D.mod1)
hist(residuals(PO4.D.mod1))
r.squaredGLMM(PO4.D.mod1)

AIC(PO4.A.mod1, PO4.B.mod1, PO4.C.mod1, PO4.D.mod1)
# model B is best

######
# FI #
######
hist((d2$FI)) #normal
hist(log10(d2$FI)) #... better, not by loads 

# Model A:
#   A. lmer(y ~ time + elevation + SA + depth + sample depth + (1|site))
FI.A.mod1 <- lmer(log10(FI) ~ scale(Elevation) + Fish + scale(log10(SA_m2)) 
                  + scale(log10(max_depth)) + sample_depth +  sample_num
                  + (1 | Site), data = d2)
summary(FI.A.mod1) # elevation, fish, sample depth and sample number
hist(residuals(FI.A.mod1))
r.squaredGLMM(FI.A.mod1)

# Model B:
#   B. lmer(y ~ time + elevation + SA + depth + (1|site))
FI.B.mod1 <- lmer(log10(FI) ~ scale(Elevation) + Fish + scale(log10(SA_m2)) 
                  + scale(log10(max_depth)) + sample_num
                  + (1 | Site), data = d2)
summary(FI.B.mod1)
hist(residuals(FI.B.mod1))
r.squaredGLMM(FI.B.mod1)

# Model C:
#   C. lmer(y ~ time + elevation + SA + depth + (1|site) + (1|watershed))
FI.C.mod1 <- lmer(log10(FI) ~ scale(Elevation) + Fish + scale(log10(SA_m2)) 
                  + scale(log10(max_depth)) + sample_depth +  sample_num
                  + (1 | Site) + (1| watershed), data = d2)

summary(FI.C.mod1)
hist(residuals(FI.C.mod1))
r.squaredGLMM(FI.C.mod1)

# Model D:
#   D. lmer(y ~ time + elevation + SA + depth + sample depth + (1|site) + (1|watershed))
FI.D.mod1 <- lmer(log10(FI) ~ scale(Elevation) + Fish + scale(log10(SA_m2)) 
                  + scale(log10(max_depth)) + sample_num
                  + (1 | Site) + (1| watershed), data = d2)

summary(FI.D.mod1)
hist(residuals(FI.D.mod1))
r.squaredGLMM(FI.D.mod1)

AIC(FI.A.mod1, FI.B.mod1, FI.C.mod1, FI.D.mod1)
# B is still best model 



########
# SUVA #
########

hist((d2$SUVA))
# Model A:
#   A. lmer(y ~ time + elevation + SA + depth + sample depth + (1|site))
SUVA.A.mod1 <- lmer(SUVA ~ scale(Elevation) + Fish + scale(log10(SA_m2)) 
                    + scale(log10(max_depth)) + sample_depth +  sample_num
                    + (1 | Site), data = d2)
summary(SUVA.A.mod1) # elevation, fish, sample depth and sample number
hist(residuals(SUVA.A.mod1))
r.squaredGLMM(SUVA.A.mod1)

# Model B:
#   B. lmer(y ~ time + elevation + SA + depth + (1|site))
SUVA.B.mod1 <- lmer(SUVA ~ scale(Elevation) + Fish + scale(log10(SA_m2)) 
                    + scale(log10(max_depth)) + sample_num
                    + (1 | Site), data = d2)
summary(SUVA.B.mod1)
hist(residuals(SUVA.B.mod1))
r.squaredGLMM(SUVA.B.mod1)

# Model C:
#   C. lmer(y ~ time + elevation + SA + depth + (1|site) + (1|watershed))
SUVA.C.mod1 <- lmer(SUVA ~ scale(Elevation) + Fish + scale(log10(SA_m2)) 
                    + scale(log10(max_depth)) + sample_depth +  sample_num
                    + (1 | Site) + (1| watershed), data = d2)

summary(SUVA.C.mod1)
hist(residuals(SUVA.C.mod1))
r.squaredGLMM(SUVA.C.mod1)

# Model D:
#   D. lmer(y ~ time + elevation + SA + depth + sample depth + (1|site) + (1|watershed))
SUVA.D.mod1 <- lmer(SUVA ~ scale(Elevation) + Fish + scale(log10(SA_m2)) 
                    + scale(log10(max_depth)) + sample_num
                    + (1 | Site) + (1| watershed), data = d2)

summary(SUVA.D.mod1)
hist(residuals(SUVA.D.mod1))
r.squaredGLMM(SUVA.D.mod1)

AIC(SUVA.A.mod1, SUVA.B.mod1, SUVA.C.mod1, SUVA.D.mod1)
# B is best model 


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






##############
# Water Temp #
##############
hist((w2$Temp))

# Model A:
#   A. lmer(y ~ time + elevation + SA + depth + sample depth + (1|site))
TEMP.A.mod1 <- lmer(Temp ~ scale(Elevation) + fish + scale(log10(SA_m2)) 
                    + scale(log10(max_depth)) + sample_depth +  sample_num
                    + (1 | Site), data = w2)
summary(TEMP.A.mod1) 
hist(residuals(TEMP.A.mod1))
r.squaredGLMM(TEMP.A.mod1)

# Model B:
#   B. lmer(y ~ time + elevation + SA + depth + (1|site))
TEMP.B.mod1 <- lmer(Temp ~ scale(Elevation) + fish + scale(log10(SA_m2)) 
                    + scale(log10(max_depth)) + sample_num
                    + (1 | Site), data = w2)
summary(TEMP.B.mod1)
hist(residuals(TEMP.B.mod1))
r.squaredGLMM(TEMP.B.mod1)

# Model C:
#   C. lmer(y ~ time + elevation + SA + depth + (1|site) + (1|watershed))
TEMP.C.mod1 <- lmer(Temp ~ scale(Elevation) + fish + scale(log10(SA_m2)) 
                    + scale(log10(max_depth)) + sample_depth +  sample_num
                    + (1 | Site) + (1| watershed), data = w2)

summary(TEMP.C.mod1)
hist(residuals(TEMP.C.mod1))
r.squaredGLMM(TEMP.C.mod1)

# Model D:
#   D. lmer(y ~ time + elevation + SA + depth + sample depth + (1|site) + (1|watershed))
TEMP.D.mod1 <- lmer(Temp ~ scale(Elevation) + fish + scale(log10(SA_m2)) 
                    + scale(log10(max_depth)) + sample_num
                    + (1 | Site) + (1| watershed), data = w2)

summary(TEMP.D.mod1)
hist(residuals(TEMP.D.mod1)) # not that normal
r.squaredGLMM(TEMP.D.mod1)

AIC(TEMP.A.mod1, TEMP.B.mod1, TEMP.C.mod1, TEMP.D.mod1)
# model A is best



################
# Conductivity #
################
hist((w2$Conductivity))
range(w2$Conductivity)
hist(log10(w2$Conductivity)) # better

# Model A:
#   A. lmer(y ~ time + elevation + SA + depth + sample depth + (1|site))
COND.A.mod1 <- lmer(log10(Conductivity) ~ scale(Elevation) + fish + scale(log10(SA_m2)) 
                    + scale(log10(max_depth)) + sample_depth +  sample_num
                    + (1 | Site), data = w2)
summary(COND.A.mod1) 
hist(residuals(COND.A.mod1)) # not that normal
r.squaredGLMM(COND.A.mod1)

# Model B:
#   B. lmer(y ~ time + elevation + SA + depth + (1|site))
COND.B.mod1 <- lmer(log10(Conductivity) ~ scale(Elevation) + fish + scale(log10(SA_m2)) 
                    + scale(log10(max_depth)) + sample_num
                    + (1 | Site), data = w2)
summary(COND.B.mod1)
hist(residuals(COND.B.mod1)) # not that normal
r.squaredGLMM(COND.B.mod1)

# Model C:
#   C. lmer(y ~ time + elevation + SA + depth + (1|site) + (1|watershed))
COND.C.mod1 <- lmer(log10(Conductivity) ~ scale(Elevation) + fish + scale(log10(SA_m2)) 
                    + scale(log10(max_depth)) + sample_depth +  sample_num
                    + (1 | Site) + (1| watershed), data = w2)

summary(COND.C.mod1)
hist(residuals(COND.C.mod1))
r.squaredGLMM(COND.C.mod1)

# Model D:
#   D. lmer(y ~ time + elevation + SA + depth + sample depth + (1|site) + (1|watershed))
COND.D.mod1 <- lmer(log10(Conductivity) ~ scale(Elevation) + fish + scale(log10(SA_m2)) 
                    + scale(log10(max_depth)) + sample_num
                    + (1 | Site) + (1| watershed), data = w2)

summary(COND.D.mod1)
hist(residuals(COND.D.mod1)) # not that normal
r.squaredGLMM(COND.D.mod1)

AIC(COND.A.mod1, COND.B.mod1, COND.C.mod1, COND.D.mod1)
# model b is best 


##########
# SECCHI #
##########
#   * only model B and D work here as there is no sample/ measurement depth
#   * DF = w3 as it only has 1 row with secchi 
hist(w3$Secchi)

# Model B:
#   B. lmer(y ~ time + elevation + SA + depth + (1|site))
SECCHI.B.mod1 <- lmer(Secchi ~ scale(Elevation) + fish + scale(log10(SA_m2)) 
                      + scale(log10(max_depth)) + sample_num
                      + (1 | Site), data = w3)
summary(SECCHI.B.mod1) 
hist(residuals(SECCHI.B.mod1)) # not that normal
r.squaredGLMM(SECCHI.B.mod1)


# Model D:
#   D. lmer(y ~ time + elevation + SA + depth + sample depth + (1|site) + (1|watershed))
SECCHI.D.mod1 <- lmer(Secchi ~ scale(Elevation) + fish + scale(log10(SA_m2)) 
                      + scale(log10(max_depth)) + sample_num
                      + (1 | Site) + (1| watershed), data = w3)

summary(SECCHI.D.mod1)
hist(residuals(SECCHI.D.mod1)) # not that normal
r.squaredGLMM(SECCHI.D.mod1)

AIC(SECCHI.B.mod1, SECCHI.D.mod1)

###################
# PAR Attenuation #
###################
#   * only model B and D work here as there is no sample/ measurement depth
#   * DF = w3 as it only has 1 row with attenuation of PAR  
hist(w3$Light.at.surface)

# Model B:
#   B. lmer(y ~ time + elevation + SA + depth + (1|site))
PAR_AT.B.mod1 <- lmer(Light.at.surface ~ scale(Elevation) + fish + scale(log10(SA_m2)) 
                      + scale(log10(max_depth)) + sample_num
                      + (1 | Site), data = w3)
summary(PAR_AT.B.mod1) 
hist(residuals(PAR_AT.B.mod1)) # not that normal
r.squaredGLMM(PAR_AT.B.mod1)


# Model D:
#   D. lmer(y ~ time + elevation + SA + depth + sample depth + (1|site) + (1|watershed))
PAR_AT.D.mod1 <- lmer(Light.at.surface ~ scale(Elevation) + fish + scale(log10(SA_m2)) 
                      + scale(log10(max_depth)) + sample_num
                      + (1 | Site) + (1| watershed), data = w3)

summary(PAR_AT.D.mod1)
hist(residuals(PAR_AT.D.mod1)) # not that normal
r.squaredGLMM(PAR_AT.D.mod1)

AIC(PAR_AT.B.mod1, PAR_AT.D.mod1)


#################
# Change in RHO #
#################
#   * only model B and D work here as there is no sample/ measurement depth
#   * DF = CP as it only has 1 row with RHO 
hist(cp$rho_change.hypo_sur.)


# Model B:
#   B. lmer(y ~ time + elevation + SA + depth + (1|site))
RHO.B.mod1 <- lmer(rho_change.hypo_sur. ~ scale(Elevation) + Fish + scale(log10(SA_m2)) 
                   + scale(log10(max_depth)) + sample_num
                   + (1 | Site), data = cp)
summary(RHO.B.mod1) 
hist(residuals(RHO.B.mod1)) # not that normal
r.squaredGLMM(RHO.B.mod1)

plot(RHO.B.mod1)


# Model D:
#   D. lmer(y ~ time + elevation + SA + depth + sample depth + (1|site) + (1|watershed))
RHO.D.mod1 <- lmer(rho_change.hypo_sur. ~ scale(Elevation) + Fish + scale(log10(SA_m2)) 
                   + scale(log10(max_depth)) + sample_num
                   + (1 | Site) + (1| watershed), data = cp)

summary(RHO.D.mod1)
hist(residuals(RHO.D.mod1)) # not that normal
r.squaredGLMM(RHO.D.mod1)

AIC(RHO.B.mod1, RHO.D.mod1)


##################################################
##################################################
# ABIOTIC FACTORS : doc, TDN, TDP, NO3, PO4, FI, SUVA, Temp, chlor-a, conductivity, 
#    secchi, light at surface, rho change

##################################################
# What relevant factors of abiotic qualities of lakes are related to ##################################################
##################################################
# ABIOTIC FACTORS : doc, TDN, TDP, NO3, PO4, FI, SUVA, Temp, chlor-a, conductivity, 
#    secchi, light at surface, rho change

##################################################
# What relevant factors of abiotic qualities of lakes are related to elevation?: 
#    1. doc.D.mod1 (positive trend for fish)
#    2. TDN.B.mod1 (only a trend)
#    3. NO3.D.mod1 (only a trend, sig. for sample number)
#    4. PO4.A.mod1 (trend of sample location, sig. for fish, and sample number)
#    5. FI.B.mod1 (trend of elevation, sig. for surface area, max depth and sample number)
#    6. TEMP.A.mod1 (sig. for sample depth and sample number)
#    7. COND.B.mod1 (sig. for sample number)
#    8. PAR_AT.B.mod1 
#    9. RHO.B.mod1 (sig. for surface area, max depth and sample number)

# Other sig drivers not related to ##################################################
##################################################
# ABIOTIC FACTORS : doc, TDN, TDP, NO3, PO4, FI, SUVA, Temp, chlor-a, conductivity, 
#    secchi, light at surface, rho change

##################################################
# What relevant factors of abiotic qualities of lakes are related to elevation?: 
#    1. doc.D.mod1 (positive trend for fish)
#    2. TDN.B.mod1 (only a trend)
#    3. NO3.D.mod1 (only a trend, sig. for sample number)
#    4. PO4.A.mod1 (trend of sample location, sig. for fish, and sample number)
#    5. FI.B.mod1 (trend of elevation, sig. for surface area, max depth and sample number)
#    6. TEMP.A.mod1 (sig. for sample depth and sample number)
#    7. COND.B.mod1 (sig. for sample number)
#    8. PAR_AT.B.mod1 
#    9. RHO.B.mod1 (sig. for surface area, max depth and sample number)

# Other sig drivers not related to elevation?
#    1. SUVA.B.mod1 (trend of fish, sig. for sample number)
#    2. SECCHI.B.mod1 (sig for max depth)

?
#    1. SUVA.B.mod1 (trend of fish, sig. for sample number)
#    2. SECCHI.B.mod1 (sig for max depth)


#    1. doc.D.mod1 (positive trend for fish)
#    2. TDN.B.mod1 (only a trend)
#    3. NO3.D.mod1 (only a trend, sig. for sample number)
#    4. PO4.A.mod1 (trend of sample location, sig. for fish, and sample number)
#    5. FI.B.mod1 (trend of elevation, sig. for surface area, max depth and sample number)
#    6. TEMP.A.mod1 (sig. for sample depth and sample number)
#    7. COND.B.mod1 (sig. for sample number)
#    8. PAR_AT.B.mod1 
#    9. RHO.B.mod1 (sig. for surface area, max depth and sample number)

# Other sig drivers not related to elvation?
#    1. SUVA.B.mod1 (trend of fish, sig. for sample number)
#    2. SECCHI.B.mod1 (sig for max depth)



table(d$sample_num)

hist(d$sample_num)
d3$FI

## plot pannel of lakes by visit number and by depth
# ABIOTIC FACTORS : doc, TDN, TDP, NO3, PO4, FI, SUVA, Temp, chlor-a, conductivity, 
#    secchi, light at surface, rho change
ggplot(data = d3, # data
         aes(x = sample_num, # aesthetics
             y = TDN_mg_L,
             color = sample_depth)) +
  geom_line() +  geom_point() + labs(y = "TDN", x = "Sample number", colour = "Depth") + 
  facet_grid(~Site) 

ggsave("TDN_time.pdf",
       scale = 2, width = 24, height = 5, units = c("cm"), dpi = 300) 

  
ggplot(data = d3, # data
       aes(x = sample_num, # aesthetics
           y = TDP_mg_L,
           color = sample_depth)) +
  geom_line() +  geom_point() + labs(y = "TDP", x = "Sample number", colour = "Depth") + 
  facet_grid(~Site) 

ggsave("TDP_time.pdf",
       scale = 2, width = 24, height = 5, units = c("cm"), dpi = 300) 
  
  
ggplot(data = d3, # data
       aes(x = sample_num, # aesthetics
           y = NO3_mg_L,
           color = sample_depth)) +
  geom_line() +  geom_point() + labs(y = "NO3", x = "Sample number", colour = "Depth") + 
  facet_grid(~Site) 

ggsave("NO3_time.pdf",
       scale = 2, width = 24, height = 5, units = c("cm"), dpi = 300) 


ggplot(data = d3, # data
       aes(x = sample_num, # aesthetics
           y = FI,
           color = sample_depth)) +
  geom_line() +  geom_point() + labs(y = "FI", x = "Sample number", colour = "Depth") + 
  facet_grid(~Site) 

ggsave("FI_time.pdf",
       scale = 2, width = 24, height = 5, units = c("cm"), dpi = 300) 

d3$SUVA
ggplot(data = d3, # data
       aes(x = sample_num, # aesthetics
           y = SUVA,
           color = sample_depth)) +
  geom_line() +  geom_point() + labs(y = "SUVA", x = "Sample number", colour = "Depth") + 
  facet_grid(~Site) 

ggsave("SUVA_time.pdf",
       scale = 2, width = 24, height = 5, units = c("cm"), dpi = 300) 

d3$DOC_mg_L
ggplot(data = d3, # data
       aes(x = sample_num, # aesthetics
           y = DOC_mg_L,
           color = sample_depth)) +
  geom_line() +  geom_point() + labs(y = "DOC", x = "Sample number", colour = "Depth") + 
  facet_grid(~Site) 

ggsave("DOC_time.pdf",
       scale = 2, width = 24, height = 5, units = c("cm"), dpi = 300) 


d3$chla
ggplot(data = d3, # data
       aes(x = sample_num, # aesthetics
           y = chla,
           color = sample_depth)) +
  geom_line() +  geom_point() + labs(y = "chla", x = "Sample number", colour = "Depth") + 
  facet_grid(~Site) 

ggsave("chla_time.pdf",
       scale = 2, width = 24, height = 5, units = c("cm"), dpi = 300) 



w4$Temp
ggplot(data = w4, # data
       aes(x = sample_num, # aesthetics
           y = Temp,
           color = sample_depth)) +
  geom_line() +  geom_point() + labs(y = "Water Temp", x = "Sample number", colour = "Depth") + 
  facet_grid(~Site) 

ggsave("Temp_time.pdf",
       scale = 2, width = 24, height = 5, units = c("cm"), dpi = 300) 

w4$Conductivity
ggplot(data = w4, # data
       aes(x = sample_num, # aesthetics
           y = Conductivity,
           color = sample_depth)) +
  geom_line() +  geom_point() + labs(y = "Conductivity", x = "Sample number", colour = "Depth") + 
  facet_grid(~Site) 

ggsave("CONDC_time.pdf",
       scale = 2, width = 24, height = 5, units = c("cm"), dpi = 300) 

w4$Light.at.surface
cp$Secchi
ggplot(data = cp2, # data
       aes(x = sample_num, # aesthetics
           y = Secchi)) +
  geom_line() +  geom_point() + labs(y = "Secchi", x = "Sample number") + 
  facet_grid(~Site) 

ggsave("Secchi_time.pdf",
       scale = 2, width = 24, height = 5, units = c("cm"), dpi = 300) 



cp$rho_change.hypo_sur.
ggplot(data = cp2, # data
       aes(x = sample_num, # aesthetics
           y = rho_change.hypo_sur.)) +
  geom_line() +  geom_point() + labs(y = "Stability", x = "Sample number") + 
  facet_grid(~Site) 

ggsave("Stability_time.pdf",
       scale = 2, width = 24, height = 5, units = c("cm"), dpi = 300) 


ggplot(data = cp2, # data
       aes(x = sample_num, # aesthetics
           y = rho_change.hypo_sur.)) +
  geom_line() +  geom_point() + labs(y = "Stability", x = "Sample number") + 
  facet_grid(~Site) 

ggsave("Stability_time.pdf",
       scale = 2, width = 24, height = 5, units = c("cm"), dpi = 300) 


cp2$attentuation
ggplot(data = cp2, # data
       aes(x = sample_num, # aesthetics
           y = attentuation)) +
  geom_line() +  geom_point() + labs(y = "PAR attenuation", x = "Sample number") + 
  facet_grid(~Site) 

ggsave("atten_time.pdf",
       scale = 2, width = 24, height = 5, units = c("cm"), dpi = 300) 









ggplot(data = d3, # data
       aes(x = sample_num, # aesthetics
           y = chla,
           color = sample_depth)) +
  geom_line() +  geom_point() + # geom 
  facet_grid(~Site) 




ggplot(data = d3, # data
       aes(x = sample_num, # aesthetics
           y = chla,
           color = sample_depth)) +
  geom_point() + # geom 
  facet_grid(~Site) + # scale
  stat_smooth(method = "lm", se = TRUE)

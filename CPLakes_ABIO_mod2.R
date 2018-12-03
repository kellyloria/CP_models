###################
# 2016 Comp Lakes 
# Univariate analysis of abiotic charactersitics of alpine lakes along an elevation gradient
###################
# load needed packages 

library(lme4)
library(lmerTest)# for p-value
library(MuMIn) # forr squared
library(nlme)
library(car)


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

head(w)
w2 <- subset(w, sample_depth=="sur" | sample_depth=="hypo",
             select= Date:Elevation)
summary(w2$sample_depth)
w2$Site


#######################################
# Relevant factors for abiotic qualities of lakes and their dataframes: 
# DF= w: Temp, chlor-a, conductivity, secchi, light at surface
# DF= d: doc, TDN, TDP, NO3, PO4, FI, SUVA
# DF= cp: rho change

###########
### DOC ###
###########

# interactions
DOC.i.mod1 <- lmer(log10(DOC_mg_L +1) ~ scale(SA_m2) + scale(max_depth) + scale(Elevation) + 
                    epi + Fish + scale(Elevation)*epi + scale(Elevation)*Fish + 
                    scale(Elevation)*scale(max_depth) + (1|sample_num) + (1|Site), data = d2)

summary(DOC.i.mod1)
hist(residuals(DOC.i.mod1))

DOC.i.mod2 <- lmer(DOC_mg_L ~ scale(SA_m2) + scale(max_depth) + scale(Elevation) + 
                   epi + Fish + scale(Elevation)*epi + scale(Elevation)*Fish + 
                   scale(Elevation)*scale(max_depth) + (1|sample_num) + (1|Site), data = d2)

summary(DOC.i.mod2)
hist(residuals(DOC.i.mod2))
# no significant effects of interactions 

DOC.mod3 <- lmer(DOC_mg_L ~ scale(SA_m2) + scale(max_depth) + scale(Elevation) + 
                   epi + Fish + (1|sample_num) + (1|Site), data = d2)

summary(DOC.mod3)


DOC.mod4 <- lmer(log10(DOC_mg_L +1) ~ scale(SA_m2) + scale(max_depth) + scale(Elevation) + 
                   epi + Fish + (1|sample_num) + (1|Site), data = d2)

summary(DOC.mod4)
hist(residuals(DOC.mod4))
vif(DOC.mod4)

AIC(DOC.mod3, DOC.mod4)

###########
### TDN ###
###########


TDN.i.mod1 <- lmer(TDN_mg_L ~ scale(SA_m2) + scale(max_depth) + scale(Elevation) + 
                   epi + Fish + scale(Elevation)*epi + scale(Elevation)*Fish + 
                   scale(Elevation)*scale(max_depth) + (1|sample_num) + (1|Site), data = d2)

summary(TDN.i.mod1)
hist(residuals(TDN.i.mod1))


TDN.mod2 <- lmer(TDN_mg_L ~ scale(SA_m2) + scale(max_depth) + scale(Elevation) + 
                   epi + Fish + (1|sample_num) + (1|Site), data = d2)

summary(TDN.mod2)
hist(residuals(TDN.mod2))
vif(TDN.mod2)

###########
### NO3 ###
###########
hist((d2$NO3_mg_L))
hist(log10(d2$NO3_mg_L +1)) # not better
hist(asin(sqrt(d2$NO3_mg_L)))

N03.i.mod1 <- lmer(NO3_mg_L ~ scale(SA_m2) + scale(max_depth) + scale(Elevation) + 
                   epi + Fish + scale(Elevation)*epi + scale(Elevation)*Fish + 
                   scale(Elevation)*scale(max_depth) + (1|sample_num) + (1|Site), data = d2)

summary(N03.i.mod1)
hist(residuals(N03.i.mod1))



NO3.mod2 <- lmer(NO3_mg_L ~ scale(SA_m2) + scale(max_depth) + scale(Elevation) + 
                   epi + Fish + (1|sample_num) + (1|Site), data = d2)

summary(NO3.mod2)
hist(residuals(NO3.mod2))
vif(NO3.mod2)

N03.i.mod3 <- lmer(asin(sqrt(NO3_mg_L)) ~ scale(SA_m2) + scale(max_depth) + scale(Elevation) + 
                   epi + Fish + scale(Elevation)*epi + scale(Elevation)*Fish + 
                   scale(Elevation)*scale(max_depth) + (1|sample_num) + (1|Site), data = d2)

summary(N03.i.mod3)
hist(residuals(N03.i.mod3))


N03.i.mod3a <- lmer(asin(sqrt(NO3_mg_L)) ~ scale(SA_m2) + scale(max_depth) + scale(Elevation) + 
                   epi + Fish + scale(Elevation)*Fish + 
                   scale(Elevation)*scale(max_depth) + (1|sample_num) + (1|Site), data = d2)

summary(N03.i.mod3a)
hist(residuals(N03.i.mod3a))

N03.i.mod3b <- lmer(asin(sqrt(NO3_mg_L)) ~ scale(SA_m2) + scale(max_depth) + scale(Elevation) + 
                    epi + Fish + scale(Elevation)*scale(max_depth) + (1|sample_num) + (1|Site), data = d2)

summary(N03.i.mod3b)
hist(residuals(N03.i.mod3b))

AIC(N03.i.mod3, N03.i.mod3a, N03.i.mod3b)

N03.mod4 <- lmer(asin(sqrt(NO3_mg_L)) ~ scale(SA_m2) + scale(max_depth) + scale(Elevation) + 
                   epi + Fish  + (1|sample_num) + (1|Site), data = d2)

summary(N03.mod4)
hist(residuals(N03.mod4))
vif(N03.mod4)

###########
### PO4 ###
###########
d$PO4_mg_L
P04.i.mod1 <- lmer(PO4_mg_L ~ scale(SA_m2) + scale(max_depth) + scale(Elevation) + 
                   epi + Fish + scale(Elevation)*epi + scale(Elevation)*Fish + 
                   scale(Elevation)*scale(max_depth) + (1|sample_num) + (1|Site), data = d2)

summary(P04.i.mod1)
hist(residuals(P04.i.mod1))

P04.i.mod1a <- lmer(PO4_mg_L ~ scale(SA_m2) + scale(max_depth) + scale(Elevation) + 
                   epi + Fish + scale(Elevation)*Fish + 
                   scale(Elevation)*scale(max_depth) + (1|sample_num) + (1|Site), data = d2)

summary(P04.i.mod1a)
hist(residuals(P04.i.mod1a))

AIC(P04.mod1, P04.mod1a)

PO4.mod2 <- lmer(PO4_mg_L ~ scale(SA_m2) + scale(max_depth) + scale(Elevation) + 
                   epi + Fish + (1|sample_num) + (1|Site), data = d2)

summary(PO4.mod2)
hist(residuals(PO4.mod2))
vif(PO4.mod2)

##########
### FI ###
##########
hist((d2$FI)) #normal
hist(log10(d2$FI))

FI.i.mod1 <- lmer(log10(FI) ~ scale(SA_m2) + scale(max_depth) + scale(Elevation) + 
                   epi + Fish + scale(Elevation)*epi + scale(Elevation)*Fish + 
                   scale(Elevation)*scale(max_depth) + (1|sample_num) + (1|Site), data = d2)

summary(FI.i.mod1)
hist(residuals(FI.i.mod1))
vif(FI.i.mod1)


FI.mod2 <- lmer(log10(FI) ~ scale(SA_m2) + scale(max_depth) + scale(Elevation) + 
                   epi + Fish + (1|sample_num) + (1|Site), data = d2)

summary(FI.mod2)
hist(residuals(FI.mod2))
vif(FI.mod2)

#############
### SUVA ###
############
hist((d2$SUVA)) #normal

SUVA.i.mod1 <- lmer(SUVA ~ scale(SA_m2) + scale(max_depth) + scale(Elevation) + 
                  epi + Fish + scale(Elevation)*epi + scale(Elevation)*Fish + 
                  scale(Elevation)*scale(max_depth) + (1|sample_num) + (1|Site), data = d2)

summary(SUVA.i.mod1)
hist(residuals(SUVA.i.mod1))

SUVA.mod2 <- lmer(SUVA ~ scale(SA_m2) + scale(max_depth) + scale(Elevation) + 
                  epi + Fish + (1|sample_num) + (1|Site), data = d2)

summary(SUVA.mod2)
hist(residuals(SUVA.mod2))
vif(SUVA.mod2)

###############
### Chlor-a ###
##############

hist((d2$chla))
hist(log10(d2$chla +1)) 

chla.i.mod1 <- lmer(log10(chla +1) ~ scale(SA_m2) + scale(max_depth) + scale(Elevation) + 
                    epi + Fish + scale(Elevation)*epi + scale(Elevation)*Fish + 
                    scale(Elevation)*scale(max_depth) + (1|sample_num) + (1|Site), data = d2)

summary(chla.i.mod1)
hist(residuals(chla.i.mod1))

chla.i.mod3 <- lmer(log10(chla +1) ~ scale(SA_m2) + scale(max_depth) + scale(Elevation) + 
                    epi + Fish + scale(Elevation)*epi  + 
                    scale(Elevation)*scale(max_depth) + (1|sample_num) + (1|Site), data = d2)

summary(chla.i.mod3)
hist(residuals(chla.i.mod3))

chla.mod3a <- lmer(log10(chla +1) ~ scale(SA_m2) + scale(max_depth) + scale(Elevation) + 
                    epi + Fish + scale(Elevation)*epi + (1|sample_num) + (1|Site), data = d2)

summary(chla.mod3a)
hist(residuals(chla.mod3a))


chla.mod2 <- lmer(log10(chla +1) ~ scale(SA_m2) + scale(max_depth) + scale(Elevation) + 
                    epi + Fish + (1|sample_num) + (1|Site), data = d2)

summary(chla.mod2)
hist(residuals(chla.mod2))
vif(chla.mod2)

###########
### SO4 ###
###########

hist((d2$SO4_mg_L))
hist(log10(d2$SO4_mg_L +1)) 

SO4.i.mod1 <- lmer(log10(SO4_mg_L +1) ~ scale(SA_m2) + scale(max_depth) + scale(Elevation) + 
                    epi + Fish + scale(Elevation)*epi + scale(Elevation)*Fish + 
                    scale(Elevation)*scale(max_depth) + (1|sample_num) + (1|Site), data = d2)

summary(SO4.i.mod1)
hist(residuals(SO4.i.mod1))

SO4.i.mod1a <- lmer(log10(SO4_mg_L +1) ~ scale(SA_m2) + scale(max_depth) + scale(Elevation) + 
                   epi + Fish + scale(Elevation)*Fish + 
                   scale(Elevation)*scale(max_depth) + (1|sample_num) + (1|Site), data = d2)

summary(SO4.i.mod1a)
hist(residuals(SO4.i.mod1a))


SO4.i.mod1b <- lmer(log10(SO4_mg_L +1) ~ scale(SA_m2) + scale(max_depth) + scale(Elevation) + 
                    epi + Fish + scale(Elevation)*Fish + (1|sample_num) + (1|Site), data = d2)

summary(SO4.i.mod1b)
hist(residuals(SO4.i.mod1b))

AIC(SO4.i.mod1, SO4.i.mod1a, SO4.i.mod1b)



SO4.mod2 <- lmer(log10(SO4_mg_L +1) ~ scale(SA_m2) + scale(max_depth) + scale(Elevation) + 
                   epi + Fish + (1|sample_num) + (1|Site), data = d2)

summary(SO4.mod2)
hist(residuals(SO4.mod2))


############
### Temp ###
############
# DF= w: Temp, conductivity, secchi, light at surface

hist((w2$Temp))

Temp.i.mod1 <- lmer(Temp ~ scale(SA_m2) + scale(max_depth) + scale(Elevation) + 
                   epi + Fish + scale(Elevation)*epi + scale(Elevation)*Fish + 
                   scale(Elevation)*scale(max_depth) + (1|sample_num) + (1|Site), data = d2)

summary(Temp.i.mod1)
hist(residuals(Temp.i.mod1))


Temp.mod2 <- lmer(Temp ~ scale(SA_m2) + scale(max_depth) + scale(Elevation) + 
                   epi + Fish + (1|sample_num) + (1|Site), data = w2)

summary(Temp.mod2)
hist(residuals(Temp.mod2))


############
### COND ###
############
# DF= w: Temp, conductivity, secchi, light at surface

hist((w2$fish))
hist(log10(w2$Conductivity))

Cond.i.mod1 <- lmer(log10(Conductivity) ~ scale(SA_m2) + scale(max_depth) + scale(Elevation) + 
                    epi + fish + scale(Elevation)*epi + scale(Elevation)*fish + 
                    scale(Elevation)*scale(max_depth) + (1|sample_num) + (1|Site), data = w2)

summary(Cond.i.mod1)
hist(residuals(Cond.i.mod1))

Cond.i.mod1a <- lmer(log10(Conductivity) ~ scale(SA_m2) + scale(max_depth) + scale(Elevation) + 
                    epi + fish + scale(Elevation)*epi + scale(Elevation)*fish
                   + (1|sample_num) + (1|Site), data = w2)

summary(Cond.i.mod1a)
hist(residuals(Cond.i.mod1a))

Cond.i.mod1b <- lmer(log10(Conductivity) ~ scale(SA_m2) + scale(max_depth) + scale(Elevation) + 
                     epi + fish + scale(Elevation)*fish
                   + (1|sample_num) + (1|Site), data = w2)

summary(Cond.i.mod1b)
hist(residuals(Cond.i.mod1b))

AIC(Cond.i.mod1, Cond.i.mod1a, Cond.i.mod1b)

Cond.mod2 <- lmer(log10(Conductivity) ~ scale(SA_m2) + scale(max_depth) + scale(Elevation) + 
                    epi + fish + (1|sample_num) + (1|Site), data = w2)

summary(Cond.mod2)
hist(residuals(Cond.mod2))

##############
### SECCHI ###
##############
# DF= w: Temp, conductivity, secchi, light at surface

hist((w2$Secchi))
hist(log10(w2$Conductivity))

Secchi.i.mod1 <- lmer(Secchi ~ scale(SA_m2) + scale(max_depth) + scale(Elevation) + 
                    epi + fish + scale(Elevation)*epi + scale(Elevation)*fish + 
                    scale(Elevation)*scale(max_depth) + (1|sample_num) + (1|Site), data = w2)

summary(Secchi.i.mod1)
hist(residuals(Secchi.i.mod1))

Secchi.i.mod1a <- lmer(Secchi ~ scale(SA_m2) + scale(max_depth) + scale(Elevation) + 
                      epi + fish + scale(Elevation)*fish + 
                      scale(Elevation)*scale(max_depth) + (1|sample_num) + (1|Site), data = w2)

summary(Secchi.i.mod1a)
hist(residuals(Secchi.i.mod1a))

AIC(Secchi.i.mod1, Secchi.i.mod1a)


Secchi.mod2 <- lmer(Secchi ~ scale(SA_m2) + scale(max_depth) + scale(Elevation) + 
                    epi + fish + (1|sample_num) + (1|Site), data = w2)

summary(Secchi.mod2)
hist(residuals(Secchi.mod2))

#########################
### Light attenuation ###
#########################

PAR_atten.mod1 <- lmer(Light.at.surface ~ scale(SA_m2) + scale(max_depth) + scale(Elevation) + 
                         fish + scale(Elevation)*fish + scale(Elevation)*scale(max_depth) + (1|sample_num) + (1|Site), data = w2)

summary(PAR_atten.mod1)
hist(residuals(PAR_atten.mod1))


PAR_atten.mod2 <- lmer(Light.at.surface ~ scale(SA_m2) + scale(max_depth) + scale(Elevation) + 
                         fish + (1|sample_num) + (1|Site), data = w2)

summary(PAR_atten.mod2)
hist(residuals(PAR_atten.mod2))



###########
### RHO ###
###########
hist(cp$rho_change.hypo_sur.)
hist(log10(cp$rho_change.hypo_sur. +1))
cp$f
RHO.i.mod1 <- lmer(rho_change.hypo_sur. ~ scale(SA_m2) + scale(max_depth) + scale(Elevation) + Fish + 
                   scale(Elevation)*Fish + scale(Elevation)*scale(max_depth) + (1|sample_num) + (1|Site), data = cp)

summary(RHO.i.mod1)
hist(residuals(RHO.i.mod1))

RHO.i.mod1a <- lmer(rho_change.hypo_sur. ~ scale(SA_m2) + scale(max_depth) + scale(Elevation) + Fish +
                    scale(Elevation)*scale(max_depth) + (1|sample_num) + (1|Site), data = cp)

summary(RHO.i.mod1a)
hist(residuals(RHO.i.mod1a))
AIC(RHO.i.mod1, RHO.i.mod1a)

RHO.mod2 <- lmer(rho_change.hypo_sur. ~ scale(SA_m2) + scale(max_depth) + scale(Elevation) + 
                   Fish + (1|sample_num) + (1|Site), data = cp)

summary(RHO.mod2)
hist(residuals(RHO.mod2))

# mer.tools for vif with mixed effect models 
# don't remove terms with VIF -- 
# VIF just tells us if  
#just tell us whether terms are just playing well together -- 
# do it with non interaction models

# kick Piet a paired down version of the code once finished 



############## SUMMARY ########################### 
##################################################
# ABIOTIC FACTORS : doc, TDN, TDP, NO3, PO4, SO4, FI, SUVA, Temp, chlor-a, conductivity, 
#    secchi, light at surface, rho change

##################################################
# What relevant factors of abiotic qualities of lakes are related to elevation?: 
#    1. DOC: DOC.mod4 (no sig. interactions)(model seems better when DOC log10+1 transformed)
#    2. TDN: TDN.mod2 (no sig. interactions)
#    3. NO3: N03.i.mod3b (trend: for interaction for max depth* elevation; fish)(sig. for elevation)
#    4. PO4: P04.i.mod1a (sig. for interaction for scale(Elevation):Fish; scale(max_depth):scale(Elevation))
#    5. SO4: SO4.i.mod1a (sig. for interaction for scale(Elevation):Fish)(sig. for surface area and elevation)
#    7. Water temperature: Temp.mod2 (no sig. interactions)(sig. for elevation, epi and max depth)
#    8. Conductivity: Cond.i.mod1b (sig. for interaction for scale(Elevation):Fish)(sig. for elevation and epi)
#    9. Secchi depth: Secchi.i.mod1a (sig. for interaction for scale(max_depth):scale(Elevation), trend for scale(Elevation):fish)
#   10. PAR attenuation: PAR_atten.mod2 (no sig. interactions)(sig. for Elevation)
#   11. Water stability: RHO.i.mod1a (sig. for interaction for scale(max_depth):scale(Elevation)(sig. for elevation and mex depth)

# Other sig drivers not related to elevation?
#    1. FI: FI.mod2 (no sig. interactions)(Sig. for SA_m2, max depth, epi)
#    2. Chlorophyll-a: chla.mod2 (no sig. interactions)(sig. for epi)

# not significant for anything
#    1. SUVA 

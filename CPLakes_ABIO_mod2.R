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
library(gridExtra)

library(ggplot2)
library(broom)
library(dplyr)
library(dotwhisker)
library(coefplot)



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


summary(d2$Site)
d3 <- subset(d2, Site=="Albion" | Site=="GL4" | Site=="Lost Lake" 
             | Site=="Blue Lake" | Site=="Isabelle Lake" | Site=="Snowbank Lake" 
             | Site=="Diamond Lake" | Site=="MudLake" | Site=="Jasper Lake" 
             | Site=="Upper Diamond Lake "| Site=="Forest Lake"| Site=="Lion Lake 2"
             | Site=="Pear Reservoir" | Site=="Yankee Doodle " | Site=="GL1" 
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
w3$Site

w4 <- subset(w2, Site=="Albion" | Site=="GL4" | Site=="Lost Lake" | Site=="RedRock Lake" 
             | Site=="Blue Lake" | Site=="Isabelle Lake" | Site=="Snowbank Lake" 
             | Site=="Diamond Lake" | Site=="MitchellLake" | Site=="Jasper Lake" 
             | Site=="Upper Diamond Lake "| Site=="Forest Lake"| Site=="Lion Lake 2"
             | Site=="Pear Reservoir" | Site=="Yankee Doodle " | Site=="GL1" 
             | Site=="Long Lake" | Site=="Red Deer Lake", 
             select= Date:Elevation)
summary(w4$Site)

w4m <- subset(w2, Site=="Albion" | Site=="GL4" | Site=="Lost Lake" | Site=="Red Rock Lake" 
             | Site=="BlueLake" | Site=="Isabelle Lake" | Site=="Snowbank Lake" 
             | Site=="Diamond Lake" | Site=="Mitchell Lake" | Site=="Jasper Lake" 
             | Site=="Upper Diamond Lake "| Site=="Forest Lake"| Site=="Lion Lake 2"
             | Site=="Pear Reservoir" | Site=="Yankee Doodle " | Site=="GL1" 
             | Site=="Long Lake" | Site=="Mud Lake" | Site=="Red Deer Lake", 
             select= Date:Elevation)
summary(w4m$Site)


cp2 <- subset(cp, Site=="Albion" | Site=="GL4" | Site=="Lost Lake" 
              | Site=="Blue Lake" | Site=="Isabelle Lake" | Site=="Snowbank Lake" 
              | Site=="Diamond Lake" | Site=="MudLake" | Site=="Jasper Lake" 
              | Site=="Upper Diamond Lake"| Site=="Forest Lake"| Site=="Lion Lake 2"
              | Site=="Pear Reservoir " | Site=="Yankee Doodle Lake" | Site=="GL1" 
              | Site=="Long Lake " | Site=="Red Deer Lake " | Site=="GL1", 
              select= Site:Sort)
summary(cp2$Site)

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


# no significant effects of interactions 


DOC.mod4 <- lmer(log10(DOC_mg_L +1) ~ scale(SA_m2) + scale(max_depth) + scale(Elevation) + 
                   epi + Fish + (1|sample_num) + (1|Site), data = d2)

summary(DOC.mod4) # best model
hist(residuals(DOC.mod4))


summary(d2$Site)
r.squaredGLMM(DOC.mod4)

###########
### TDN ###
###########


TDN.i.mod1 <- lmer(TDN_mg_L ~ scale(SA_m2) + scale(max_depth) + scale(Elevation) + 
                   epi + Fish + scale(Elevation)*epi + scale(Elevation)*Fish + 
                   scale(Elevation)*scale(max_depth) + (1|sample_num) + (1|Site), data = d3)

summary(TDN.i.mod1)
hist(residuals(TDN.i.mod1))
# no sig. interactions

TDN.mod2 <- lmer(TDN_mg_L ~ scale(SA_m2) + scale(max_depth) + scale(Elevation) + 
                   epi + Fish + (1|sample_num) + (1|Site), data = d2)

summary(TDN.mod2) # best model
hist(residuals(TDN.mod2))
vif(TDN.mod2)
r.squaredGLMM(TDN.mod2)

###########
### NO3 ###
###########
hist((d2$NO3_mg_L))
hist(log10(d2$NO3_mg_L +1)) # not better

hist((d3$NO3_mg_L))
hist(log10(d3$NO3_mg_L +1))


N03.i.mod1 <- lmer(log10(NO3_mg_L +1) ~ scale(SA_m2) + scale(max_depth) + scale(Elevation) + 
                   epi + Fish + scale(Elevation)*epi + scale(Elevation)*Fish + 
                   scale(Elevation)*scale(max_depth) + (1|sample_num) + (1|Site), data = d3)

summary(N03.i.mod1)
hist(residuals(N03.i.mod1))


N03.i.mod2 <- lmer(log10(NO3_mg_L +1) ~ scale(SA_m2) + scale(max_depth) + scale(Elevation) + 
                     epi + Fish + scale(Elevation)*epi + scale(Elevation)*Fish + 
                     scale(Elevation)*scale(max_depth) + (1|sample_num) + (1|Site), data = d2)

summary(N03.i.mod2)
hist(residuals(N03.i.mod2))


N03.i.mod2a <- lmer(log10(NO3_mg_L +1) ~ scale(SA_m2) + scale(max_depth) + scale(Elevation) + 
                     epi + Fish + scale(Elevation)*scale(max_depth) + (1|sample_num) + (1|Site), data = d2)

summary(N03.i.mod2a) # best model 
hist(residuals(N03.i.mod2a))
r.squaredGLMM(N03.i.mod2a)

AIC(N03.i.mod3, N03.i.mod3a, N03.i.mod3b)

N03.mod4 <- lmer(log10(NO3_mg_L +1) ~ scale(SA_m2) + scale(max_depth) + scale(Elevation) +
                   epi + (1|sample_num) + (1|Site), data = d2)

summary(N03.mod4)
hist(residuals(N03.mod4))
vif(N03.mod4)

###########
### PO4 ###
###########
hist(d2$PO4_mg_L)
P04.i.mod1 <- lmer(PO4_mg_L ~ scale(SA_m2) + scale(max_depth) + scale(Elevation) + 
                   epi + Fish + scale(Elevation)*epi + scale(Elevation)*Fish + 
                   scale(Elevation)*scale(max_depth) + (1|sample_num) + (1|Site), data = d2)

summary(P04.i.mod1)
hist(residuals(P04.i.mod1))
?vif()
vif(P04.i.mod1)

P04.i.mod1a <- lmer(PO4_mg_L ~ scale(SA_m2) + scale(max_depth) + scale(Elevation) + 
                   epi + Fish + scale(Elevation)*Fish + 
                   scale(Elevation)*scale(max_depth) + (1|sample_num) + (1|Site), data = d2)

summary(P04.i.mod1a)
hist(residuals(P04.i.mod1a))

P04.i.mod2 <- lmer(PO4_mg_L ~ scale(SA_m2) + scale(max_depth) + scale(Elevation) + 
                     epi + Fish + scale(Elevation)*epi + scale(Elevation)*Fish + 
                     scale(Elevation)*scale(max_depth) + (1|sample_num) + (1|Site), data = d3)

summary(P04.i.mod2)
hist(residuals(P04.i.mod2))

AIC(P04.mod1, P04.mod1a)

PO4.mod2 <- lmer(PO4_mg_L ~ scale(SA_m2) + scale(max_depth) + scale(Elevation) + 
                   epi + Fish + (1|sample_num) + (1|Site), data = d2)

summary(PO4.mod2) # best model 
hist(residuals(PO4.mod2))
vif(PO4.mod2)
r.squaredGLMM(PO4.mod2)


summary(d2$Site)

PO4.mod3 <- lmer(PO4_mg_L ~ scale(SA_m2) + scale(max_depth) + scale(Elevation) + 
                   epi + Fish + (1|sample_num) + (1|Site), data = d3)

summary(PO4.mod3)
hist(residuals(PO4.mod3))
vif(PO4.mod3)

#### TDP ####
hist(d2$TDP_mg_L)
hist(log10(d2$TDP_mg_L + 1))
hist(d3$TDP_mg_L) # D3 looks best
hist(log10(d3$TDP_mg_L + 1))

TDP.i.mod1a <- lmer(TDP_mg_L ~ scale(SA_m2) + scale(max_depth) + scale(Elevation) + 
                      epi + Fish + scale(Elevation)*epi + scale(Elevation)*Fish + 
                      scale(Elevation)*scale(max_depth) + (1|sample_num) + (1|Site), data = d2)

summary(TDP.i.mod1a)
hist(residuals(TDP.i.mod1a))

TDP.i.mod1b <- lmer(TDP_mg_L ~ scale(SA_m2) + scale(max_depth) + scale(Elevation) + 
                      epi + Fish + scale(Elevation)*Fish + 
                      scale(Elevation)*scale(max_depth) + (1|sample_num) + (1|Site), data = d2)

summary(TDP.i.mod1b)  
hist(residuals(TDP.i.mod1b))

TDP.i.mod1c <- lmer(TDP_mg_L ~ scale(SA_m2) + scale(max_depth) + scale(Elevation) + 
                      epi + Fish + scale(Elevation)*scale(max_depth) + (1|sample_num) + 
                      (1|Site), data = d3)

summary(TDP.i.mod1c)
hist(residuals(TDP.i.mod1c))

TDP.mod2 <- lmer(TDP_mg_L ~ scale(SA_m2) + scale(max_depth) + scale(Elevation) + 
                      epi + (1|sample_num) + (1|Site), data = d2)

summary(TDP.mod2)
hist(residuals(TDP.mod2))



##########
### FI ###
##########
hist((d2$FI)) #normal
hist(log10(d2$FI))
max(na.omit(d2$FI))
FI.i.mod1 <- lmer(log10(FI) ~ scale(SA_m2) + scale(max_depth) + scale(Elevation) + 
                   epi + Fish + scale(Elevation)*epi + scale(Elevation)*Fish + 
                   scale(Elevation)*scale(max_depth) + (1|sample_num) + (1|Site), data = d2)

summary(FI.i.mod1)
hist(residuals(FI.i.mod1))
vif(FI.i.mod1)

FI.mod2 <- lmer(log10(FI) ~ scale(SA_m2) + scale(max_depth) + scale(Elevation) + 
                   epi + Fish + (1|sample_num) + (1|Site), data = d2)

summary(FI.mod2) # best model 
hist(residuals(FI.mod2))
vif(FI.mod2)
r.squaredGLMM(FI.mod2)


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

hist((d3$chla))
hist(log10(d3$chla +1)) 

chla.i.mod1 <- lmer(log10(chla +1) ~ scale(SA_m2) + scale(max_depth) + scale(Elevation) + 
                    epi + Fish + scale(Elevation)*epi + scale(Elevation)*Fish + 
                    scale(Elevation)*scale(max_depth) + (1|sample_num) + (1|Site), data = d2)

summary(chla.i.mod1)
hist(residuals(chla.i.mod1))


chla.i.mod1a <- lmer(log10(chla +1) ~ scale(SA_m2) + scale(max_depth) + scale(Elevation) + 
                      epi + Fish + scale(Elevation)*epi + scale(Elevation)*Fish + 
                      scale(Elevation)*scale(max_depth) + (1|sample_num) + (1|Site), data = d3)

summary(chla.i.mod1a)
hist(residuals(chla.i.mod1a))

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

summary(chla.mod2) #best model
hist(residuals(chla.mod2))
r.squaredGLMM(chla.mod2)

###########
### SO4 ###
###########

hist((d2$SO4_mg_L))
hist(log10(d2$SO4_mg_L +1)) 

hist((d3$SO4_mg_L))
hist(log10(d3$SO4_mg_L +1)) 

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
                   epi +  scale(Elevation)*epi  + 
                   scale(Elevation)*scale(max_depth) + (1|sample_num) + (1|Site), data = d2)

summary(Temp.i.mod1)
hist(residuals(Temp.i.mod1))


Temp.mod2 <- lmer(Temp ~ scale(SA_m2) + scale(max_depth) + scale(Elevation) + 
                   epi + (1|sample_num) + (1|Site), data = w2)

summary(Temp.mod2)
hist(residuals(Temp.mod2))
r.squaredGLMM(Temp.mod2)
vif(Temp.mod2)

############
### COND ###
############
# DF= w: Temp, conductivity, secchi, light at surface

hist((w2$fish))
hist(log10(w2$Conductivity))

Cond.i.mod1 <- lmer(log10(Conductivity) ~ scale(SA_m2) + scale(max_depth) + scale(Elevation) + 
                    epi + fish + scale(Elevation)*epi + scale(Elevation)*fish + 
                    scale(Elevation)*scale(max_depth) + (1|sample_num) + (1|Site), data = w4m)

summary(Cond.i.mod1)
hist(residuals(Cond.i.mod1))

Cond.i.mod1a <- lmer(log10(Conductivity) ~ scale(SA_m2) + scale(max_depth) + scale(Elevation) + 
                      epi + fish+ scale(Elevation)*fish + 
                      scale(Elevation)*scale(max_depth) + (1|sample_num) + (1|Site), data = w4m)

summary(Cond.i.mod1a)
hist(residuals(Cond.i.mod1a))


Cond.i.mod1b <- lmer(log10(Conductivity) ~ scale(SA_m2) + scale(max_depth) + scale(Elevation) + 
                       epi + fish+ scale(Elevation)*fish + (1|sample_num) + (1|Site), data = w4m)

summary(Cond.i.mod1b) 
hist(residuals(Cond.i.mod1b))



Cond.mod2 <- lmer(log10(Conductivity) ~ scale(SA_m2) + scale(max_depth) + scale(Elevation) + 
                    epi + fish + (1|sample_num) + (1|Site), data = w4m)

summary(Cond.mod2) # best model
hist(residuals(Cond.mod2)) ##cut mud lake out
r.squaredGLMM(Cond.mod2)


##############
### SECCHI ###
##############
# DF= w: Temp, conductivity, secchi, light at surface

hist((w2$Secchi))
hist(log10(w2$Secchi))
hist(log10(w4$Secchi))

Secchi.i.mod1 <- lmer(log10(Secchi) ~ scale(SA_m2) + scale(max_depth) + scale(Elevation) + 
                        Fish + scale(Elevation)*Fish + scale(Elevation)*scale(max_depth) + 
                        (1|sample_num) + (1|Site), data = cp)

summary(Secchi.i.mod1)
hist(residuals(Secchi.i.mod1))

Secchi.i.mod1a <- lmer(log10(Secchi) ~ scale(SA_m2) + scale(max_depth) + scale(Elevation) + 
                        Fish + scale(Elevation)*scale(max_depth) + 
                        (1|sample_num) + (1|Site), data = cp)

summary(Secchi.i.mod1a)
hist(residuals(Secchi.i.mod1a))


Secchi.mod2 <- lmer(log10(Secchi) ~ scale(SA_m2) + scale(max_depth) + scale(Elevation) + Fish +
                  (1|sample_num) + (1|Site), data = cp)

summary(Secchi.mod2)
hist(residuals(Secchi.mod2))
r.squaredGLMM(Secchi.mod2)

#########################
### Light attenuation ###
#########################

PAR_atten.mod1 <- lmer(Light.at.surface ~ scale(SA_m2) + scale(max_depth) + scale(Elevation) + 
                         fish + scale(Elevation)*fish + scale(Elevation)*scale(max_depth) + (1|sample_num) + (1|Site), data = w4)

summary(PAR_atten.mod1)
hist(residuals(PAR_atten.mod1))


PAR_atten.mod2 <- lmer(Light.at.surface ~ scale(SA_m2) + scale(max_depth) + scale(Elevation) + 
                         fish + (1|sample_num) + (1|Site), data = w2)

summary(PAR_atten.mod2)
hist(residuals(PAR_atten.mod2))

hist(cp$visit_re)
PAR_atten.mod3 <- lmer(attentuation ~ scale(SA_m2) + scale(max_depth) + scale(Elevation) + 
                         Fish + as.factor(1|sample_num) + (1|Site), data = cp)

summary(PAR_atten.mod3)
hist(residuals(PAR_atten.mod3))
r.squaredGLMM(PAR_atten.mod3)


###########
### RHO ###
########### 
hist(cp$rho_change.hypo_sur.)
hist(log10(cp$rho_change.hypo_sur. +1))
min(cp$rho_change.hypo_sur.) #one negative value check to see if that makes sense maybe just force to zero 
# maybe subset without blue lake 
RHO.i.mod1 <- lmer(rho_change.hypo_sur. ~ scale(SA_m2) + scale(max_depth) + scale(Elevation)+ 
                   + scale(Elevation)*scale(max_depth) + (1|sample_num) + (1|Site), data = cp2)

summary(RHO.i.mod1)
hist(residuals(RHO.i.mod1))



RHO.mod2 <- lmer(rho_change.hypo_sur. ~ scale(SA_m2) + scale(max_depth) + scale(Elevation) + 
                    (1|sample_num) + (1|Site), data = cp2)

hist(cp2$delt_rho)
summary(RHO.mod2)
hist(residuals(RHO.mod2))
r.squaredGLMM(RHO.mod2)
cp2$delt_rho

RHO.mod2 <- lmer(delt_rho ~ scale(SA_m2) + scale(max_depth) + scale(Elevation) + 
                   (1|sample_num) + (1|Site), data = cp2)

summary(RHO.mod2)
hist(residuals(RHO.mod2))
r.squaredGLMM(RHO.mod2)
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

## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
# Split by Fish, Epi and Median depth 
d2$epi 
w4

# for fish interactions
d2.if0 <- subset(d2, Fish==0,
             select= Site:phyto_den)
summary(d2.if0$Fish)

d2.if1 <- subset(d2, Fish==1,
                select= Site:phyto_den)
summary(d2.if1$Fish)

w4$fish
head(w4)
w4.if0 <- subset(w4, fish==0,
                 select= Date:Elevation)
summary(w4.if0$fish)

w4.if1 <- subset(w4, fish==1,
                 select= Date:Elevation)
summary(w4.if1$fish)


w4.imd0 <- subset(w4, max_depth<9,
                  select=Date:Elevation)
summary(w4.imd0$max_depth)

w4.imd1 <- subset(w4, max_depth>9,
                  select=Date:Elevation)
summary(w4.imd1$max_depth)



# for epi interactions
d2.ie0 <- subset(d2, epi==0,
                 select= Site:phyto_den)
summary(d2.ie0$epi)

d2.ie1 <- subset(d2, epi==1,
                 select= Site:phyto_den)
summary(d2.ie1$Fish)


summary(d2$max_depth)
d2.imd0 <- subset(d2, max_depth<9,
                 select= Site:phyto_den)
summary(d2.imd0$max_depth)

d2.imd1 <- subset(d2, max_depth>9,
                  select= Site:phyto_den)
summary(d2.imd1$max_depth)



summary(d2$max_depth)
cp2.imd0 <- subset(cp2, max_depth<9,
                  select= Site:Sort)
summary(cp2.imd0$max_depth)

cp2.imd1 <- subset(cp2, max_depth>9,
                  select= Site:Sort)
summary(cp2.imd1$max_depth)


#######
# N03 #
N03.trend0.mod2a <- lmer(log10(NO3_mg_L +1) ~ scale(SA_m2) + scale(max_depth) + scale(Elevation) + 
                      epi + Fish + scale(Elevation)*scale(max_depth) + (1|sample_num) + (1|Site), data = d2)

summary(N03.trend0.mod2a) # best model 
hist(residuals(N03.trend0.mod2a))

# Nitrate interation and max depth
d2$NO3_mg_L
ggplot(d2.imd0,
       aes(x = Elevation,
           y = NO3_mg_L)) + geom_point(color='red') + 
  geom_smooth(data=d2.imd0, aes(x = Elevation, y = NO3_mg_L), fill= "red", color="grey25", 
              se = TRUE, stat = "smooth", method = "lm") + 
  labs(y = "NO3", x = "Elevation", title="Shallow (red) and deep (blue)") + 
  geom_point(data= d2.imd1, aes(x = Elevation,
                                y = NO3_mg_L), color= 'blue') + theme(plot.margin = unit(c(0.5, 1,0.5,0.5), "cm")) +
  geom_smooth(data=d2.imd1, aes(x = Elevation, y = NO3_mg_L), fill= 'blue', color= 'grey25', 
              se = TRUE, stat = "smooth", method = "lm")  
ggsave("Interaction_NO3_maxdepth.pdf",
       scale = 2, width = 7, height = 5, units = c("cm"), dpi = 300) 


#######
# PO4 #

P04.i.mod1a <- lmer(PO4_mg_L ~ scale(SA_m2) + scale(max_depth) + scale(Elevation) + 
                      epi + Fish + scale(Elevation)*Fish + 
                      scale(Elevation)*scale(max_depth) + (1|sample_num) + (1|Site), data = d2)

summary(P04.i.mod1a)
hist(residuals(P04.i.mod1a))

# PO4 fish interaction
d2$PO4_mg_L
ggplot(data = d2,
            aes(x = Elevation, 
                y = PO4_mg_L, color=as.factor(Fish))) + labs(y = "PO4", x = "Elevation") +
  geom_point() +  theme(plot.margin = unit(c(0.5, 1,0.5,0.5), "cm")) +
  geom_smooth(method = "lm") 
ggsave("Interaction_PO4_fish.pdf",
       scale = 2, width = 7, height = 5, units = c("cm"), dpi = 300) 

ggplot(data = d2,
       aes(x = Elevation, 
           y = PO4_mg_L, color=as.factor(epi))) + labs(y = "PO4", x = "Elevation") +
  geom_point() +  theme(plot.margin = unit(c(0.5, 1,0.5,0.5), "cm")) +
  geom_smooth(method = "lm") 

ggsave("Interaction_PO4_epi.pdf",
       scale = 2, width = 7, height = 5, units = c("cm"), dpi = 300) 

# PO4 interation and max depth
ggplot(d2.imd0,
             aes(x = Elevation,
                 y = PO4_mg_L)) + geom_point(color='red') + 
  geom_smooth(data=d2.imd0, aes(x = Elevation, y = PO4_mg_L), fill= "red", color="grey25", 
              se = TRUE, stat = "smooth", method = "lm") + 
  labs(y = "PO4", x = "Elevation", title="Shallow (red) and deep (blue)") + 
  geom_point(data= d2.imd1, aes(x = Elevation,
              y = PO4_mg_L), color= 'blue') + theme(plot.margin = unit(c(0.5, 1,0.5,0.5), "cm")) +
  geom_smooth(data=d2.imd1, aes(x = Elevation, y = PO4_mg_L), fill= 'blue', color= 'grey25', 
              se = TRUE, stat = "smooth", method = "lm")  
ggsave("Interaction_PO4_maxdepth.pdf",
       scale = 2, width = 7, height = 5, units = c("cm"), dpi = 300) 




# TDP #
TDP.i.mod1b <- lmer(TDP_mg_L ~ scale(SA_m2) + scale(max_depth) + scale(Elevation) + 
                      epi + Fish + scale(Elevation)*Fish + 
                      scale(Elevation)*scale(max_depth) + (1|sample_num) + (1|Site), data = d2)

summary(TDP.i.mod1b) # best model 
hist(residuals(TDP.i.mod1b))



# TDP and max depth
d2$TDP_mg_L
ggplot(d2.imd0,
       aes(x = Elevation,
           y = TDP_mg_L)) + geom_point(color='red') + 
  geom_smooth(data=d2.imd0, aes(x = Elevation, y = TDP_mg_L), fill= "red", color="grey25", 
              se = TRUE, stat = "smooth", method = "lm") + 
  labs(y = "TDP", x = "Elevation", title="Shallow (red) and deep (blue)") + 
  geom_point(data= d2.imd1, aes(x = Elevation,
                                y = TDP_mg_L), color= 'blue') + theme(plot.margin = unit(c(0.5, 1,0.5,0.5), "cm")) +
  geom_smooth(data=d2.imd1, aes(x = Elevation, y = TDP_mg_L), fill= 'blue', color= 'grey25', 
              se = TRUE, stat = "smooth", method = "lm")  
ggsave("Interaction_TDP_maxdepth.pdf",
       scale = 2, width = 7, height = 5, units = c("cm"), dpi = 300) 

d2$TDN_mg_L
ggplot(data = d2, # data
       aes(x = Elevation, # aesthetics
           y = TDN_mg_L, colour=as.factor(Fish))) + labs(y = "TDP", x = "Elevation") +
  geom_point() + theme(plot.margin = unit(c(0.5, 1,0.5,0.5), "cm"))+
  geom_smooth(method = "lm") 
ggsave("Interaction_TDP_fish.pdf",
       scale = 2, width = 7, height = 5, units = c("cm"), dpi = 300) 

################
# Conductivity #
hist(residuals(Cond.i.mod1b))
w4$Conductivity
ggplot(data = w4m, # data
       aes(x = Elevation, # aesthetics
           y = Conductivity, colour=as.factor(fish))) + labs(y = "Conductivity", x = "Elevation") +
  geom_point() + theme(plot.margin = unit(c(0.5, 1,0.5,0.5), "cm"))+
  geom_smooth(method = "lm") 
ggsave("Interaction_conduct_fish.pdf",
       scale = 2, width = 7, height = 5, units = c("cm"), dpi = 300) 


ggplot(data = w4m, # data
       aes(x = Elevation, # aesthetics
           y = Conductivity, colour=as.factor(epi))) + labs(y = "Conductivity", x = "Elevation") +
  geom_point() + theme(plot.margin = unit(c(0.5, 1,0.5,0.5), "cm"))+
  geom_smooth(method = "lm") 
ggsave("Interaction_conduct_epi.pdf",
       scale = 2, width = 7, height = 5, units = c("cm"), dpi = 300) 

##########
# Secchi #

ggplot(data = cp2, # data
       aes(x = Elevation, # aesthetics
           y = Secchi, colour=as.factor(Fish))) + labs(y = "Secchi", x = "Elevation") +
  geom_point() + theme(plot.margin = unit(c(0.5, 1,0.5,0.5), "cm"))+
  geom_smooth(method = "lm") 
ggsave("Interaction_secchi_fish.pdf",
       scale = 2, width = 7, height = 5, units = c("cm"), dpi = 300) 

# Secchi and max depth
w4$Secchi
ggplot(cp2.imd0,
       aes(x = Elevation,
           y = Secchi)) + geom_point(color='red') + 
  geom_smooth(data=cp2.imd0, aes(x = Elevation, y = Secchi), fill= "red", color="grey25", 
              se = TRUE, stat = "smooth", method = "lm") + 
  labs(y = "Secchi", x = "Elevation", title="Shallow (red) and deep (blue)") + 
  geom_point(data= cp2.imd1, aes(x = Elevation,
                                y = Secchi), color= 'blue') + theme(plot.margin = unit(c(0.5, 1,0.5,0.5), "cm")) +
  geom_smooth(data=cp2.imd1, aes(x = Elevation, y = Secchi), fill= 'blue', color= 'grey25', 
              se = TRUE, stat = "smooth", method = "lm")  
ggsave("Interaction_secchi_maxdepth.pdf",
       scale = 2, width = 7, height = 5, units = c("cm"), dpi = 300) 

summary(d2$SA_m2)


### Coeficient plot
?coefplot

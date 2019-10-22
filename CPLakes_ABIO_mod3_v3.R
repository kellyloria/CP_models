###################
# 2016 Comp Lakes 
# Univariate analysis of abiotic charactersitics of alpine lakes along an elevation gradient
###################
# Post revisions # 


d2 <- read.csv("2016_Comp_Lakes_water_chem.csv", header=T)
names(d)
w2 <- read.csv("COMPLAKES_2016data_a.csv", header=T)
names(w2)
cp <- read.csv("comp_lakes_master_data_2.csv", header=T)
names(cp)


summary(d2$Site)
d_v1 <- subset(d2, Site=="Albion" | Site=="GL4" | Site=="Lost Lake" 
               | Site=="Blue Lake" | Site=="Isabelle Lake" | Site=="Snowbank Lake" 
               | Site=="Diamond Lake" | Site=="Mud Lake" | Site=="Jasper Lake" 
               | Site=="Upper Diamond Lake "| Site=="Forest Lake"| Site=="Lion Lake 2"
               | Site=="Pear Reservoir" | Site=="Yankee Doodle " | Site=="GL1" 
               | Site=="Long Lake" | Site=="Red Deer Lake" | Site=="GL1"
               | Site=="Mitchell Lake" | Site=="Red Rock Lake",
               select= Site:NDVI_QAQC)
summary(d_v1$Site)


w_v1 <- subset(w2, Site=="Albion" | Site=="GL4" | Site=="Lost Lake" 
               | Site=="Blue Lake" | Site=="Isabelle Lake" | Site=="Snowbank Lake" 
               | Site=="Diamond Lake" | Site=="Mud Lake" | Site=="Jasper Lake" 
               | Site=="Upper Diamond Lake "| Site=="Forest Lake"| Site=="Lion Lake 2"
               | Site=="Pear Reservoir" | Site=="Yankee Doodle " | Site=="GL1" 
               | Site=="Long Lake" | Site=="Red Deer Lake" | Site=="GL1"
               | Site=="Mitchell Lake" | Site=="Red Rock Lake", select= Date:Elevation)
summary(w_v1$Site)



# water chem so d df:
# interactions
hist(d_v1$DOC_mg_L)
hist(log10(d_v1$DOC_mg_L +1))
DOC.i.mod1 <- lmer(log10(DOC_mg_L +1) ~ scale(SA_m2) + scale(max_depth) + scale(Elevation) + 
                     epi + Fish + scale(Elevation)*epi + scale(Elevation)*Fish + 
                     scale(Elevation)*scale(max_depth) + (1|sample_num) + (1|Site), data = d_v1)

summary(DOC.i.mod1)
hist(residuals(DOC.i.mod1))


# w/o interactions
DOC.i.mod1 <- lmer(log10(DOC_mg_L +1) ~ scale(SA_m2) + scale(max_depth) + scale(Elevation) + 
                     epi + Fish + (1|sample_num) + (1|Site), data = d_v1)

summary(DOC.i.mod1)  #better model 
hist(residuals(DOC.i.mod1))
r.squaredGLMM(DOC.i.mod1)

#####
# TDN 
hist(d_v1$TDN_uMOL_L)
hist(log10(d_v1$TDN_uMOL_L))
TDN.i.mod1 <- lmer(TDN_uMOL_L ~ scale(SA_m2) + scale(max_depth) + scale(Elevation) + 
                     epi + Fish + scale(Elevation)*epi + scale(Elevation)*Fish + 
                     scale(Elevation)*scale(max_depth) + (1|sample_num) + (1|Site), data = d_v1)

summary(TDN.i.mod1) # significant
hist(residuals(TDN.i.mod1))

# w/o interactions
TDN.i.mod1 <- lmer((TDN_uMOL_L) ~ scale(SA_m2) + scale(max_depth) + scale(Elevation) + 
                     epi + Fish + (1|sample_num) + (1|Site), data = d_v1)

summary(TDN.i.mod1)
hist(residuals(TDN.i.mod1))
r.squaredGLMM(TDN.i.mod1)

TDN.i.mod1 <- lmer((log10(TDN_uMOL_L +1)) ~ scale(SA_m2) + scale(max_depth) + scale(Elevation) + 
                     epi + Fish + (1|sample_num) + (1|Site), data = d_v1)

summary(TDN.i.mod1)

#####
# NO3 
hist(d_v1$NO3_uMOL)
hist(log10(d_v1$NO3_uMOL +1))
NO3.i.mod1 <- lmer(log10(NO3_uMOL +1) ~ scale(SA_m2) + scale(max_depth) + scale(Elevation) + 
                     epi + Fish + scale(Elevation)*epi + scale(Elevation)*Fish + 
                     scale(Elevation)*scale(max_depth) + (1|sample_num) + (1|Site), data = d_v1)

summary(NO3.i.mod1) 
hist(residuals(NO3.i.mod1))

# w/o interactions
NO3.i.mod1 <- lmer(log10(NO3_uMOL +1) ~ scale(SA_m2) + scale(max_depth) + scale(Elevation) + 
                     epi + Fish + (1|sample_num) + (1|Site), data = d_v1)

summary(NO3.i.mod1) # significant
hist(residuals(NO3.i.mod1))
r.squaredGLMM(NO3.i.mod1)

#####
# TDP 
hist(d_v1$TDP_uMOL)
hist(log10(d_v1$TDP_uMOL +1))
TDP.i.mod1 <- lmer(log10(TDP_uMOL +1) ~ scale(SA_m2) + scale(max_depth) + scale(Elevation) + 
                     epi + Fish + scale(Elevation)*epi + scale(Elevation)*Fish + 
                     scale(Elevation)*scale(max_depth) + (1|sample_num) + (1|Site), data = d_v1)

summary(TDP.i.mod1)
hist(residuals(TDP.i.mod1))

# w/o interactions
TDP.i.mod1 <- lmer(log10(TDP_uMOL +1) ~ scale(SA_m2) + scale(max_depth) + scale(Elevation) + 
                     epi + Fish + (1|sample_num) + (1|Site), data = d_v1)

summary(TDP.i.mod1) # not significant
hist(residuals(TDP.i.mod1)) 
r.squaredGLMM(TDP.i.mod1)

#####
# PO4 
hist(d_v1$PO4_uMOL)
PO4.i.mod1 <- lmer((PO4_uMOL) ~ scale(SA_m2) + scale(max_depth) + scale(Elevation) + 
                     epi + Fish + scale(Elevation)*epi + scale(Elevation)*Fish + 
                     scale(Elevation)*scale(max_depth) + (1|sample_num) + (1|Site), data = d_v1)

summary(PO4.i.mod1)
hist(residuals(PO4.i.mod1))


# w/o interactions
PO4.i.mod1 <- lmer((PO4_uMOL) ~ scale(SA_m2) + scale(max_depth) + scale(Elevation) + 
                     epi + Fish + (1|sample_num) + (1|Site), data = d_v1)

summary(PO4.i.mod1)
hist(residuals(PO4.i.mod1)) # best fit model
r.squaredGLMM(PO4.i.mod1)

PO4.i.mod2 <- lmer((PO4_uMOL) ~ scale(SA_m2) + scale(max_depth) + scale(Elevation) + 
                     epi + (1|sample_num) + (1|Site), data = d_v1)

summary(PO4.i.mod2)
hist(residuals(PO4.i.mod2))
plot(residuals(PO4.i.mod2))
AIC(PO4.i.mod1,PO4.i.mod2)

#####
# SO4 
hist(d_v1$SO4_uMOL)
hist(log10(d_v1$SO4_uMOL))
SO4.i.mod1 <- lmer((SO4_uMOL) ~ scale(SA_m2) + scale(max_depth) + scale(Elevation) + 
                     epi + Fish + scale(Elevation)*epi + scale(Elevation)*Fish + 
                     scale(Elevation)*scale(max_depth) + (1|sample_num) + (1|Site), data = d_v1)

summary(SO4.i.mod1)
hist(residuals(SO4.i.mod1))

# w/o interactions
SO4.i.mod1 <- lmer(log10(SO4_uMOL) ~ scale(SA_m2) + scale(max_depth) + scale(Elevation) + 
                     epi + Fish + (1|sample_num) + (1|Site), data = d_v1)

summary(SO4.i.mod1)
hist(residuals(SO4.i.mod1)) # still not pattern
r.squaredGLMM(SO4.i.mod1)

#####
# FI 
hist(d_v1$FI)
hist(log10(d_v1$SO4_uMOL))
FI.i.mod1 <- lmer((FI) ~ scale(SA_m2) + scale(max_depth) + scale(Elevation) + 
                     epi + Fish + scale(Elevation)*epi + scale(Elevation)*Fish + 
                     scale(Elevation)*scale(max_depth) + (1|sample_num) + (1|Site), data = d_v1)

summary(FI.i.mod1)
hist(residuals(FI.i.mod1))

# w/o interactions
FI.i.mod1 <- lmer((FI) ~ scale(SA_m2) + scale(max_depth) + scale(Elevation) + 
                     epi + Fish + (1|sample_num) + (1|Site), data = d_v1)

summary(FI.i.mod1)
hist(residuals(FI.i.mod1)) 
r.squaredGLMM(FI.i.mod1)

# water quality
# interactions
hist(w_v1$Chl.a)
hist(log10(w_v1$Chl.a +1))
CHLA.i.mod1 <- lmer(log10(Chl.a +1) ~ scale(SA_m2) + scale(max_depth) + scale(Elevation) + 
                     epi + fish + scale(Elevation)*epi + scale(Elevation)*fish + 
                     scale(Elevation)*scale(max_depth) + (1|sample_num) + (1|Site), data = w_v1)

summary(CHLA.i.mod1)
hist(residuals(CHLA.i.mod1))

# w/o interactions
CHLA.mod1 <- lmer(log10(Chl.a +1) ~ scale(SA_m2) + scale(max_depth) + scale(Elevation) + 
                      epi + fish + (1|sample_num) + (1|Site), data = w_v1)

summary(CHLA.mod1)
hist(residuals(CHLA.mod1))

r.squaredGLMM(CHLA.mod1)


CHLA.mod2 <- lmer(log10(Chl.a +1) ~ scale(SA_m2) + scale(max_depth) + scale(Elevation) + 
                    epi +scale(sample_num) + (1|Site), data = w_v1)

summary(CHLA.mod2)
hist(residuals(CHLA.mod2))


######
# Temp
hist(w_v1$Temp)
hist(log10(w_v1$Temp))
Temp.i.mod1 <- lmer(Temp ~ scale(SA_m2) + scale(max_depth) + scale(Elevation) + scale(Elevation)*epi  + 
                      scale(Elevation)*scale(max_depth) + (1|sample_num) + (1|Site), data = w_v1)

summary(Temp.i.mod1)
hist(residuals(Temp.i.mod1))

# w/o interactions
Temp.mod1 <- lmer((Temp) ~ scale(SA_m2) + scale(max_depth) + scale(Elevation) + epi+ (1|sample_num) + (1|Site), data = w_v1)

summary(Temp.mod1)
hist(residuals(Temp.mod1)) # significant
r.squaredGLMM(Temp.mod1)


####
# pH
# interactions
hist(w_v1$pH)
hist(log10(w_v1$Chl.a +1))
pH.i.mod1 <- lmer(pH ~ scale(SA_m2) + scale(max_depth) + scale(Elevation) + 
                      epi + fish + scale(Elevation)*epi + scale(Elevation)*fish + 
                      scale(Elevation)*scale(max_depth) + (1|sample_num) + (1|Site), data = w_v1)

summary(pH.i.mod1)
hist(residuals(pH.i.mod1))

# w/o interactions
pH.mod1 <- lmer(pH ~ scale(SA_m2) + scale(max_depth) + scale(Elevation) + 
                    epi + fish + (1|sample_num) + (1|Site), data = w_v1)

summary(pH.mod1)
hist(residuals(pH.mod1))
r.squaredGLMM(pH.mod1)

######
# SP.C
# interactions
hist(w_v1$St.C)
hist(log10(w_v1$St.C))
STC.i.mod1 <- lmer(log10(St.C) ~ scale(SA_m2) + scale(max_depth) + scale(Elevation) + 
                    epi + fish + scale(Elevation)*epi + scale(Elevation)*fish + 
                    scale(Elevation)*scale(max_depth) + (1|sample_num) + (1|Site), data = w_v1)

summary(STC.i.mod1)
hist(residuals(STC.i.mod1))

# w/o interactions
STC.mod1 <- lmer(log10(St.C) ~ scale(SA_m2) + scale(max_depth) + scale(Elevation) + 
                  epi + fish + (1|sample_num) + (1|Site), data = w_v1)

summary(STC.mod1) #significant
hist(residuals(STC.mod1))
r.squaredGLMM(STC.mod1)


Secchi.i.mod1 <- lmer(log10(Secchi) ~ scale(SA_m2) + scale(max_depth) + scale(Elevation) + 
                        Fish + scale(Elevation)*Fish + scale(Elevation)*scale(max_depth) + 
                        (1|sample_num) + (1|Site), data = cp)

summary(Secchi.i.mod1)
hist(residuals(Secchi.i.mod1))

# FI data 
precip <- read.csv("lake_weather_stations_dist_3.csv", header=T)
summary(precip)


precip2 <- read.csv("lake_weather_stations_dist_4.csv", header=T)
summary(precip2)

# anverage stratification across all lakes from 1st to 3rd visit
(0.0259049-0.0163927)/0.0259049 # 36.71% increase

# anverage chlor across all lakes from 1st to 3rd visit
(5.110-5.0992)/5.110 # 2.4 %

#Precip plot 1

precip1_plot <- ggplot(precip, aes(x = precip_1, y = FI_1)) + ylab("FI") +
  geom_point(aes(colour = as.factor(sample), shape=as.factor(Fish)), alpha=0.75) + #geom_smooth(method = "glm", colour = 'grey35', fill='grey85') +
  theme_bw() + theme(text = element_text(size=14),
                          legend.position="none", 
                          plot.margin = unit(c(0.15, 0.1,0.15,1.2), "cm")) + 
  scale_x_continuous(limits=c(0,20), breaks=seq(0,20,4))+
  scale_y_continuous(limits=c(1.25,1.50),breaks=seq(1.25,1.5,0.05)) + 
  scale_colour_manual(values = c("dodgerblue4", "sky blue")) +
  scale_shape_manual(values = c(15, 17))


precip2_plot <- ggplot(precip, aes(x = precip_2, y = FI_2)) + ylab("FI") +
  geom_point(aes(colour = as.factor(sample), shape=as.factor(Fish)), alpha=0.75) + #geom_smooth(method = "glm", colour = 'grey35', fill='grey85') +
  theme_bw() + theme(text = element_text(size=14),
                          legend.position="none", 
                          plot.margin = unit(c(0.15, 0.1,0.15,1.2), "cm")) + 
  scale_x_continuous(limits=c(0,20), breaks=seq(0,20,4)) + 
  scale_y_continuous(limits=c(1.25,1.50),breaks=seq(1.25,1.5,0.05)) + 
  scale_colour_manual(values = c("dodgerblue4", "sky blue")) +
  scale_shape_manual(values = c(15, 17))


precip3_plot <- ggplot(precip, aes(x = precip_3, y = FI_3)) + ylab("FI") +
  geom_point(aes(colour = as.factor(sample), shape=as.factor(Fish)), alpha=0.75) + #geom_smooth(method = "glm", colour = 'grey35', fill='grey85') +
  theme_bw() + theme(text = element_text(size=14),
                          legend.position="none", 
                          plot.margin = unit(c(0.15, 0.1,0.15,1.2), "cm")) + 
  scale_x_continuous(limits=c(0,20), breaks=seq(0,20,4)) + 
  scale_y_continuous(limits=c(1.25,1.50),breaks=seq(1.25,1.5,0.05)) + 
  scale_colour_manual(values = c("dodgerblue4", "sky blue")) +
  scale_shape_manual(values = c(15, 17))

precip_pannel <- grid.arrange(precip1_plot, precip2_plot, precip3_plot, nrow = 1, ncol=3)


precip4_plot.A <- ggplot(precip2, aes(x = visit, y = FI)) + ylab("FI") +
  geom_point(aes(colour = as.factor(sample), shape=as.factor(Fish))) #geom_smooth(method = "glm", colour = 'grey35', fill='grey85') +
  #theme_bw() + theme(text = element_text(size=14), 
                    # plot.margin = unit(c(0.15, 0.1,0.15,1.2), "cm")) +  
  #scale_colour_manual(values = c("dodgerblue4", "sky blue")) +
  #scale_shape_manual(values = c(15, 17)) + 
  #scale_y_continuous(limits=c(1.25,1.50),breaks=seq(1.25,1.5,0.05)) +
  #facet_grid(.~Site)

precip4_plot.b <- ggplot(precip2, aes(x = visit, y = precip)) + ylab("FI") 
  geom_line() + #geom_smooth(method = "glm", colour = 'grey35', fill='grey85') +
  #theme_bw() + theme(text = element_text(size=14), 
                    # plot.margin = unit(c(0.15, 0.1,0.15,1.2), "cm")) +  
  #scale_colour_manual(values = c("dodgerblue4", "sky blue")) +
  #scale_shape_manual(values = c(15, 17)) + 
  #scale_y_continuous(limits=c(0, 20),breaks=seq(0,20,4)) +
  #facet_grid(.~Site)

#scale_x_continuous(limits=c(0,20), breaks=seq(0,20,4)) 
#scale_y_continuous(limits=c(1.25,1.50),breaks=seq(1.25,1.5,0.05)) + 



precip4_plot <- ggplot(precip2) +
  geom_point(aes(y =FI, x = visit, colour = as.factor(sample), shape=as.factor(Fish))) +
  geom_line(aes(y =precip, x = visit)) +
  theme_bw() + theme(text = element_text(size=14), 
                     plot.margin = unit(c(0.15, 0.1,0.15,1.2), "cm")) +  
  #scale_colour_manual(values = c("dodgerblue4", "sky blue", "grey30")) +
  #scale_y_continuous(limits=c(1.25,1.50),breaks=seq(1.25,1.5,0.05), 
                     #sec.axis = sec_axis(~.*7.5, name = "Precip")) +
  #scale_shape_manual(values = c(15, 17)) + facet_grid(.~Site) 

library(grid)
grid.newpage()
grid.draw(rbind(ggplotGrob(precip4_plot.A), ggplotGrob(precip4_plot.b), size = "last"))






# sube set out deep lakes 
# blue pear and red dear
blue <- precip2 %>% 
  subset(Site =="Blue Lake")





pdf("blue1plot.pdf", width = 4.5, height = 4) 

#precip and FI
par(mar = c(5,5,2,5))
with(blue, plot(visit, precip, type="l", col="dodgerblue4", 
             ylab="Precip (mm)",
             ylim=c(0,20)))

par(new = T)
with(blue, plot(visit, FI, pch=16, axes=F, xlab=NA, ylab=NA, cex=1.2, 
                col = blue$sample))
axis(side = 4)
mtext(side = 4, line = 3, 'FI')
legend("topleft",
       legend=c("Precip.", "FI hyp.", "FI sur."),
       lty=c(1,0,0), pch=c(NA, 16, 16), col=c("dodgerblue4", "black", "red"), cex=0.8)

dev.off()
#precip and chlora
pdf("blue12plot.pdf", width = 4.5, height = 4) 

par(mar = c(5,5,2,5))
with(blue, plot(visit, precip, type="l", col="dodgerblue4", 
                ylab="Precip (mm)",
                ylim=c(0,20)))

par(new = T)
with(blue, plot(visit, Chl.a, pch=17, axes=F, xlab=NA, ylab=NA, cex=1.2, 
                col = blue$sample))
axis(side = 4)
mtext(side = 4, line = 3, 'Chlor-a')
legend("topleft",
       legend=c("Precip.", "Chla hyp.", "Chla sur."),
       lty=c(1,0, 0), pch=c(NA, 17, 17), col=c("dodgerblue4", "black", "red"),cex=0.8)
dev.off()


# precip and zooplankton 

pdf("blue3plot.pdf", width = 4.5, height = 4) 
par(mar = c(5,5,2,5))
with(blue, plot(visit, precip, type="l", col="dodgerblue4", 
                ylab="Precip (mm)",
                ylim=c(0,20)))

par(new = T)
with(blue, plot(visit, zoop_den, pch=18, axes=F, xlab=NA, ylab=NA, cex=1.2))
axis(side = 4)
mtext(side = 4, line = 3, 'Zooplankton density')
legend("topleft",
       legend=c("Precip.", "Zoo."),
       lty=c(1,0), pch=c(NA, 18), col=c("dodgerblue4", "black"),cex=0.8)

dev.off()



pear <- precip2 %>% 
  subset(Site =="Pear Reservoir")

pdf("pear1plot.pdf", width = 4.5, height = 4) 

#precip and FI
par(mar = c(5,5,2,5))
with(pear, plot(visit, precip, type="l", col="dodgerblue4", 
                ylab="Precip (mm)",
                ylim=c(0,20)))

par(new = T)
with(pear, plot(visit, FI, pch=16, axes=F, xlab=NA, ylab=NA, cex=1.2, 
                col = blue$sample))
axis(side = 4)
mtext(side = 4, line = 3, 'FI')

dev.off()
#precip and chlora
pdf("pear2plot.pdf", width = 4.5, height = 4) 

par(mar = c(5,5,2,5))
with(pear, plot(visit, precip, type="l", col="dodgerblue4", 
                ylab="Precip (mm)",
                ylim=c(0,20)))

par(new = T)
with(pear, plot(visit, Chl.a, pch=17, axes=F, xlab=NA, ylab=NA, cex=1.2, 
                col = blue$sample))
axis(side = 4)
mtext(side = 4, line = 3, 'Chlor-a')

dev.off()


# precip and zooplankton 

pdf("pear3plot.pdf", width = 4.5, height = 4) 
par(mar = c(5,5,2,5))
with(pear, plot(visit, precip, type="l", col="dodgerblue4", 
                ylab="Precip (mm)",
                ylim=c(0,20)))

par(new = T)
with(pear, plot(visit, zoop_den, pch=18, axes=F, xlab=NA, ylab=NA, cex=1.2))
axis(side = 4)
mtext(side = 4, line = 3, 'Zooplankton density')

dev.off()


reddeer <- precip2 %>% 
  subset(Site =="Red Deer Lake")

pdf("reddeer1plot.pdf", width = 4.5, height = 4) 

#precip and FI
par(mar = c(5,5,2,5))
with(reddeer, plot(visit, precip, type="l", col="dodgerblue4", 
                ylab="Precip (mm)",
                ylim=c(0,20)))

par(new = T)
with(reddeer, plot(visit, FI, pch=16, axes=F, xlab=NA, ylab=NA, cex=1.2, 
                col = blue$sample))
axis(side = 4)
mtext(side = 4, line = 3, 'FI')

dev.off()
#precip and chlora
pdf("reddeer2plot.pdf", width = 4.5, height = 4) 

par(mar = c(5,5,2,5))
with(reddeer, plot(visit, precip, type="l", col="dodgerblue4", 
                ylab="Precip (mm)",
                ylim=c(0,20)))

par(new = T)
with(reddeer, plot(visit, Chl.a, pch=17, axes=F, xlab=NA, ylab=NA, cex=1.2, 
                col = blue$sample))
axis(side = 4)
mtext(side = 4, line = 3, 'Chlor-a')

dev.off()


# precip and zooplankton 

pdf("reddeer3plot.pdf", width = 4.5, height = 4) 
par(mar = c(5,5,2,5))
with(reddeer, plot(visit, precip, type="l", col="dodgerblue4", 
                ylab="Precip (mm)",
                ylim=c(0,20)))

par(new = T)
with(reddeer, plot(visit, zoop_den, pch=18, axes=F, xlab=NA, ylab=NA, cex=1.2))
axis(side = 4)
mtext(side = 4, line = 3, 'Zooplankton density')

dev.off()









blue_plot <- ggplot(blue, aes(x = visit, y = FI)) + ylab("FI") +
  geom_point(aes(colour = as.factor(sample), shape=as.factor(Fish)), alpha=0.75) + 
  geom_line() +
  #geom_smooth(method = "glm", colour = 'grey35', fill='grey85') +
  theme_bw() + theme(text = element_text(size=14),
                     legend.position="none", 
                     plot.margin = unit(c(0.15, 0.1,0.15,1.2), "cm")) + 
  scale_x_continuous(limits=c(1,3), breaks=seq(1,3,1)) + 
  scale_y_continuous(limits=c(1.25,1.50),breaks=seq(1.25,1.5,0.05)) + 
  scale_colour_manual(values = c("dodgerblue4", "sky blue")) +
  scale_shape_manual(values = c(15, 17))


#Ratio of zooplankton to chlor-a over time
hist(precip2$zoop_den)
hist(log10(precip2$zoop_den +1))


DOC.i.mod1 <- lmer(log10(Chl.a +1) ~ scale(precip) + (1|visit) + (1|Site), data = precip2)

summary(DOC.i.mod1)
hist(residuals(DOC.i.mod1))
r.squaredGLMM(DOC.i.mod1)

DOC.i.mod1 <- lmer(FI ~ scale(precip) + (1|sample) + (1|visit) + (1|Site), data = precip2)

summary(DOC.i.mod1)
hist(residuals(DOC.i.mod1))

DOC.i.mod1 <- lmer(log10(Strat +1) ~ scale(precip) + (1|visit) + (1|Site), data = precip2)

summary(DOC.i.mod1)
hist(residuals(DOC.i.mod1))

DOC.i.mod1 <- lmer(log10(zoop_den +1) ~ scale(precip) + (1|visit) + (1|Site), data = precip2)

summary(DOC.i.mod1)
hist(residuals(DOC.i.mod1))

# double check stratification 
w2 <- read.csv("COMPLAKES_2016data_a.csv", header=T)
names(w2)

w_sur <- w2 %>% 
  subset(rho_deph == "sur")
summary(w_sur)

w_hyp <- w2 %>% 
  subset(rho_deph == "bottom")
summary(w_hyp)

hist(w2$delta_rho)
hist(log10(w2$delta_rho +1))

RHO.mod <- lmer(log10(delta_rho+1) ~ scale(SA_m2) + scale(max_depth) + scale(Elevation) + 
                   (1|sample_num) + (1|Site), data = w2)

summary(RHO.mod)
hist(residuals(RHO.mod))
r.squaredGLMM(RHO.i.mod2)
cp2$delt_rho


RHO.i.mod <- lmer(log10(delta_rho+1) ~ scale(SA_m2) + scale(max_depth) + scale(Elevation)+ 
                     + scale(Elevation)*scale(max_depth) + 
                    (1|sample_num) + (1|Site), data = w2)

summary(RHO.i.mod)
hist(residuals(RHO.i.mod))
r.squaredGLMM(RHO.i.mod)


RHO.i.mod2 <- lmer(log10(delta_rho+1) ~ scale(sample_num) + (1|Site), data = w2)

summary(RHO.i.mod2)
hist(residuals(RHO.i.mod2))


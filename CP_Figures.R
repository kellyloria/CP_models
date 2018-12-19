library(lme4)
library(lmerTest)# for p-value
library(MuMIn) # forr squared
library(nlme)
library(car)
library(ggplot2)
library(gridExtra)

d <- read.csv("2016_Comp_Lakes_water_chem.csv", header=T)
names(d)
w <- read.csv("COMPLAKES_2016data_a.csv", header=T)
names(w)
cp <- read.csv("comp_lakes_master_data_2.csv", header=T)
names(cp)

d2 <- subset(d, sample_depth=="sur" | sample_depth=="hypo",
             select= Site:phyto_den)
summary(d2$sample_depth)

head(w)
w2 <- subset(w, sample_depth=="sur" | sample_depth=="hypo",
             select= Date:Elevation)
summary(w2$sample_depth)
w2$Site



### water temp ###
CP_temp_plot <- ggplot(w2, aes(x = Elevation, y = Temp)) + ylab("Water temperature") + 
  geom_point(aes(colour = as.factor(epi), shape=as.factor(fish))) + geom_smooth(method = "glm", colour = 'grey30') + 
  theme_classic() + theme(text = element_text(size=14), 
                          axis.title.x=element_blank(),
                          axis.text.x=element_blank(), legend.position="none", 
                          plot.margin = unit(c(0.15, 0,0.5,0.75), "cm")) + 
  scale_x_continuous(breaks=seq(2400,3600,150)) + 
  scale_y_continuous(limits=c(2,19), breaks=seq(0,21,3)) + 
  scale_colour_manual(values = c("dodgerblue4", "sky blue"))


#ggsave("CP_temp_plot.pdf",
#       scale = 2, width = 6, height = 5.25, units = c("cm"), dpi = 300)


### conductivity ###
CP_cond_plot <- ggplot(w4m, aes(x = Elevation, y = log10(Conductivity))) + ylab("log10(Conductivity)") +
  geom_point(aes(colour = as.factor(epi), shape=as.factor(fish))) + geom_smooth(method = "glm", colour = 'grey30') + 
  theme_classic() + theme(text = element_text(size=14),
                          axis.title.x=element_blank(),
                          axis.text.x=element_blank(), legend.position="none", 
                          plot.margin = unit(c(0.15, 0,0.5,0.5), "cm")) + 
  scale_x_continuous(breaks=seq(2400,3600,150)) + 
  scale_y_continuous(limits=c(0.6,2.1), breaks=seq(0,3.25,0.25)) + 
  scale_colour_manual(values = c("dodgerblue4", "sky blue"))

#ggsave("CP_cond_plot.pdf",
#       scale = 2, width = 6, height = 5.25, units = c("cm"), dpi = 300)


### rho ##
CP_rho_plot <- ggplot(cp, aes(x = Elevation, y = rho_change.hypo_sur.)) + ylab("Δρ") +
  geom_point(aes(shape=as.factor(Fish))) + geom_smooth(method = "glm", colour = 'grey30') + 
  theme_classic() + theme(text = element_text(size=14), 
                          axis.text.x = element_text(angle=45, hjust=1.2),
                          axis.title.x=element_blank(),
                          legend.position="none", 
                          plot.margin = unit(c(0.15, 0,0.5,0.5), "cm")) + 
  scale_x_continuous(breaks=seq(2400,3600,150)) + 
  scale_y_continuous(limits=c(0,0.78), breaks=seq(0,1,0.15))


#ggsave("CP_rho_plot.pdf",
#       scale = 2, width = 6.25, height = 5.5, units = c("cm"), dpi = 300)

CP_abio_pannel1 <- grid.arrange(CP_temp_plot, CP_cond_plot, CP_rho_plot, nrow = 3)

ggsave("CP_abio_pannel.pdf", CP_abio_pannel1, scale = 2, width = 6, height = 10, units = c("cm"), dpi = 300)


### TDN ###
CP_TDN_plot <- ggplot(d2, aes(x = Elevation, y = TDN_mg_L)) + ylab("TDN") +
  geom_point(aes(colour = as.factor(epi), shape=as.factor(Fish))) + geom_smooth(method = "glm", colour = 'grey30') + 
  theme_classic() + theme(text = element_text(size=14),
                          axis.title.x=element_blank(),
                          axis.text.x=element_blank(), legend.position="none", 
                          plot.margin = unit(c(0.15, 0,0.5,0.85), "cm")) + 
  scale_x_continuous(breaks=seq(2400,3600,150)) + 
  scale_y_continuous(limits=c(0,0.45), breaks=seq(0,0.5,0.1)) + 
  scale_colour_manual(values = c("dodgerblue4", "sky blue"))


#ggsave("CP_temp_plot.pdf",
#       scale = 2, width = 6, height = 5.25, units = c("cm"), dpi = 300)


### TDP ###
range(na.omit(d2$TDP_mg_L))
CP_TDP_plot <- ggplot(d2, aes(x = Elevation, y = TDP_mg_L)) + ylab("TDP") +
  geom_point(aes(colour = as.factor(epi), shape=as.factor(Fish)))  + 
  theme_classic() + theme(text = element_text(size=14),
                          axis.title.x=element_blank(),
                          axis.text.x=element_blank(), legend.position="none", 
                          plot.margin = unit(c(0.15, 0,0.5,0.5), "cm")) + 
  scale_x_continuous(breaks=seq(2400,3600,150)) + 
  scale_y_continuous(breaks=seq(0,0.02,0.003)) + 
  scale_colour_manual(values = c("dodgerblue4", "sky blue"))  


### DOC ###
range(na.omit(d2$DOC_mg_L))
CP_DOC_plot <- ggplot(d2, aes(x = Elevation, y = DOC_mg_L)) + ylab("DOC") +
  geom_point(aes(colour = as.factor(epi), shape=as.factor(Fish)))  + geom_smooth(method = "glm", colour = 'grey30') +
  theme_classic() + theme(text = element_text(size=14), 
                          axis.text.x = element_text(angle=45, hjust=1.2),
                          axis.title.x=element_blank(),
                          legend.position="none", 
                          plot.margin = unit(c(0.15, 0,0.5,1), "cm")) + 
  scale_x_continuous(breaks=seq(2400,3600,150)) + 
  scale_y_continuous(limits=c(0,9), breaks=seq(0,9,1.5)) + 
  scale_colour_manual(values = c("dodgerblue4", "sky blue"))



CP_abio_pannel2 <- grid.arrange(CP_temp_plot, CP_TDN_plot, CP_cond_plot, CP_TDP_plot, CP_rho_plot, CP_DOC_plot, nrow = 3, ncol=2)

ggsave("CP_abio_pannel2.pdf", CP_abio_pannel2, scale = 2, width = 11, height = 10, units = c("cm"), dpi = 300)

### NO3 ###
range(na.omit(d2$NO3_mg_L))
CP_NO3_plot <- ggplot(d2, aes(x = Elevation, y = NO3_mg_L)) + ylab("NO3") +
  geom_point(aes(colour = as.factor(epi), shape=as.factor(Fish)))  + geom_smooth(method = "glm", colour = 'grey30') +
  theme_classic() + theme(text = element_text(size=14),
                          axis.title.x=element_blank(),
                          axis.text.x=element_blank(), legend.position="none", 
                          plot.margin = unit(c(0.15, 0,0.5,1.01), "cm")) + 
  scale_x_continuous(breaks=seq(2400,3600,150)) + 
  scale_y_continuous(limits=c(-0.02,1), breaks=seq(0.000,1, 0.25)) + 
  scale_colour_manual(values = c("dodgerblue4", "sky blue"))


#ggsave("CP_temp_plot.pdf",
#       scale = 2, width = 6, height = 5.25, units = c("cm"), dpi = 300)


### P04 ###
CP_PO4_plot <- ggplot(d2, aes(x = Elevation, y = PO4_mg_L)) + ylab("PO4") +
  geom_point(aes(colour = as.factor(epi), shape=as.factor(Fish)))  + geom_smooth(method = "glm", colour = 'grey30') +
  theme_classic() + theme(text = element_text(size=14),
                          axis.title.x=element_blank(),
                          axis.text.x=element_blank(), legend.position="none", 
                          plot.margin = unit(c(0.15, 0,0.5,0.5), "cm")) + 
  scale_x_continuous(breaks=seq(2400,3600,150)) + 
  scale_y_continuous(breaks=seq(0,0.004,0.00064)) + 
  scale_colour_manual(values = c("dodgerblue4", "sky blue"))

#ggsave("CP_cond_plot.pdf",
#       scale = 2, width = 6, height = 5.25, units = c("cm"), dpi = 300)


### FI ###
CP_FI_plot <- ggplot(d2, aes(x = Elevation, y = FI)) + ylab("FI") +
  geom_point(aes(colour = as.factor(epi), shape=as.factor(Fish))) +
  theme_classic() + theme(text = element_text(size=14), 
                          axis.text.x = element_text(angle=45, hjust=1.2),
                          axis.title.x=element_blank(),
                          legend.position="none", 
                          plot.margin = unit(c(0.15, 0,0.5,1.2), "cm")) + 
  scale_x_continuous(breaks=seq(2400,3600,150)) + 
  scale_y_continuous(limits=c(1.25,1.55),breaks=seq(0,2,0.07)) + 
  scale_colour_manual(values = c("dodgerblue4", "sky blue"))



CP_abio_pannel3 <- grid.arrange(CP_temp_plot, CP_TDN_plot, CP_NO3_plot, CP_cond_plot, CP_TDP_plot, 
                                CP_P04_plot, CP_rho_plot, CP_DOC_plot, CP_FI_plot, nrow = 3, ncol=3)

ggsave("CP_abio_pannel3.pdf", CP_abio_pannel3, scale = 2, width = 16, height = 10, units = c("cm"), dpi = 300)


### chlora ###
CP_chlor_plot <- ggplot(d2, aes(x = Elevation, y = log10(d2$chla +1))) + ylab("log10(Chla+1)") +
  geom_point(aes(colour = as.factor(epi), shape=as.factor(Fish))) +
  theme_classic() + theme(text = element_text(size=14),
                          axis.title.x=element_blank(),
                          axis.text.x=element_blank(), legend.position="none", 
                          plot.margin = unit(c(0.15, 0.15,0.5,0.85), "cm")) + 
  scale_x_continuous(breaks=seq(2400,3600,150)) + 
  scale_y_continuous(limits=c(0.1,1.42), breaks=seq(0.000,1.4, 0.25)) + 
  scale_colour_manual(values = c("dodgerblue4", "sky blue"))


#ggsave("CP_temp_plot.pdf",
#       scale = 2, width = 6, height = 5.25, units = c("cm"), dpi = 300)


### secchi ###
range(na.omit(cp$Secchi))

CP_secchi_plot <- ggplot(cp, aes(x = Elevation, y = Secchi)) + ylab("Secchi depth") +
  geom_point(aes(shape=as.factor(Fish)))  + geom_smooth(method = "glm", colour = 'grey30') +
  theme_classic() + theme(text = element_text(size=14),
                          axis.title.x=element_blank(),
                          axis.text.x=element_blank(), legend.position="none", 
                          plot.margin = unit(c(0.15, 0.15,0.5,1.1), "cm")) + 
  scale_x_continuous(breaks=seq(2400,3600,150)) + 
  scale_y_continuous(breaks=seq(0,8,1.5)) 


#ggsave("CP_cond_plot.pdf",
#       scale = 2, width = 6, height = 5.25, units = c("cm"), dpi = 300)


### Kpar ###
range(na.omit(cp$attentuation))

CP_kpar_plot <- ggplot(cp, aes(x = Elevation, y = attentuation)) + ylab("Kpar") +
  geom_point(aes(shape=as.factor(Fish)))  + geom_smooth(method = "glm", colour = 'grey30') +
  theme_classic() + theme(text = element_text(size=14), 
                          axis.text.x = element_text(angle=45, hjust=1.2),
                          axis.title.x=element_blank(),
                          legend.position="none", 
                          plot.margin = unit(c(0.15, 0.15,0.5,1), "cm")) + 
  scale_x_continuous(breaks=seq(2400,3600,150)) + 
  scale_y_continuous(limits=c(0,3.5),breaks=seq(0,3.5,0.75)) 



CP_abio_pannel4 <- grid.arrange(CP_temp_plot, CP_TDN_plot, CP_NO3_plot, CP_chlor_plot, CP_cond_plot, CP_TDP_plot, 
                                CP_PO4_plot, CP_secchi_plot, CP_rho_plot, CP_DOC_plot, CP_FI_plot, CP_kpar_plot, nrow = 3, ncol=4)

ggsave("CP_abio_pannel4.pdf", CP_abio_pannel4, scale = 2, width = 17, height = 10, units = c("cm"), dpi = 300)

########### panel of bio figures #############

PHYTD.mod1 <- lmer(log10(phyto_den) ~ scale(SA_m2) + scale(max_depth) + scale(Elevation) + 
                     epi + Fish + (1|sample_num) + (1|Site), data = d2)

summary(PHYTD.mod1)
hist(residuals(PHYTD.mod1))
range(na.omit(log10(d2$phyto_den)))
CP_phyden_plot <- ggplot(d2, aes(x = Elevation, y = log10(phyto_den))) + ylab("log10(phyto density)") +
  geom_point(aes(colour = as.factor(epi), shape=as.factor(Fish))) + geom_smooth(method = "glm", colour = 'grey30') +
  theme_classic() + theme(text = element_text(size=14),                          
                          axis.title.x=element_blank(),
                          axis.text.x=element_blank(),
                          legend.position="none", 
                          plot.margin = unit(c(0.15, 0,0.5,1.2), "cm")) + 
  scale_x_continuous(breaks=seq(2400,3600,150)) + 
  scale_y_continuous(limits=c(1.5,4.7),breaks=seq(0,5,0.75)) + 
  scale_colour_manual(values = c("dodgerblue4", "sky blue"))

range(na.omit(log10(cp$zoop_den2 +1)))
CP_zooden_plot <- ggplot(cp, aes(x = Elevation, y = log10(zoop_den2 +1))) + ylab("log10(zoo density+1)") +
  geom_point(aes(shape=as.factor(Fish))) + geom_smooth(method = "glm", colour = 'grey30') +
  theme_classic() + theme(text = element_text(size=14), 
                          axis.text.x = element_text(angle=45, hjust=1.2),
                          axis.title.x=element_blank(),
                          legend.position="none", 
                          plot.margin = unit(c(0.15, 0,0.5,1.2), "cm")) + 
  scale_x_continuous(breaks=seq(2400,3600,150)) + 
  scale_y_continuous(limits=c(0,2),breaks=seq(0,2.25,0.45)) + 
  scale_colour_manual(values = c("dodgerblue4", "sky blue"))

range(na.omit(cp$Ave_size))
CP_avesize_plot <- ggplot(cp, aes(x = Elevation, y = Ave_size)) + ylab("Average zoo size") +
  geom_point(aes(shape=as.factor(Fish))) + geom_smooth(method = "glm", colour = 'grey30') +
  theme_classic() + theme(text = element_text(size=14),                          
                          axis.title.x=element_blank(),
                          axis.text.x=element_blank(),
                          legend.position="none", 
                          plot.margin = unit(c(0.15, 0,0.5,1.2), "cm")) + 
  scale_x_continuous(breaks=seq(2400,3600,150)) + 
  scale_y_continuous(limits=c(0,2.25),breaks=seq(0,5,0.5)) + 
  scale_colour_manual(values = c("dodgerblue4", "sky blue"))

cp$prop_gravid
CP_fecund_plot <- ggplot(cp, aes(x = Elevation, y = prop_gravid)) + ylab("Proportion gravid zoo") +
  geom_point(aes(shape=as.factor(Fish))) +
  theme_classic() + theme(text = element_text(size=14), 
                          axis.text.x = element_text(angle=45, hjust=1.2),
                          axis.title.x=element_blank(),
                          legend.position="none", 
                          plot.margin = unit(c(0.15, 0.15,0.5,1.2), "cm")) + 
  scale_x_continuous(breaks=seq(2400,3600,150)) + 
  scale_y_continuous(limits=c(0,0.8),breaks=seq(0,1,0.25)) + 
  scale_colour_manual(values = c("dodgerblue4", "sky blue"))


CP_bio_pannel <- grid.arrange(CP_phyden_plot, CP_avesize_plot, CP_zooden_plot, CP_fecund_plot, nrow = 2, ncol=2)

ggsave("CP_bio_pannel.pdf", CP_bio_pannel, scale = 2, width = 8.5, height = 6, units = c("cm"), dpi = 300)



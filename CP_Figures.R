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
             select= Site:NDVI_QAQC)
summary(d2$sample_depth)
summary(d2$Site)

head(w)
w2 <- subset(w, sample_depth=="sur" | sample_depth=="hypo",
             select= Date:Elevation)
summary(w2$sample_depth)
w2$Site

d2$NDVI_QAQC

d3 <- subset(d2, Site=="Albion" | Site=="GL4" | Site=="Lost Lake" 
             | Site=="Blue Lake" | Site=="Isabelle Lake" | Site=="Snowbank Lake" 
             | Site=="Diamond Lake" | Site=="MudLake" | Site=="Jasper Lake" 
             | Site=="Upper Diamond Lake "| Site=="Forest Lake"| Site=="Lion Lake 2"
             | Site=="Pear Reservoir" | Site=="Yankee Doodle " | Site=="GL1" 
             | Site=="Long Lake" | Site=="Red Deer Lake" | Site=="Mitchell Lake" 
             | Site=="Red Rock Lake", 
             select= Site:NDVI_QAQC)
summary(d3$Site)

cp2 <- subset(cp, Site=="Albion" | Site=="GL4" | Site=="Lost Lake" 
             | Site=="Blue Lake" | Site=="Isabelle Lake" | Site=="Snowbank Lake" 
             | Site=="Diamond Lake" | Site=="Mud Lake" | Site=="Jasper Lake" 
             | Site=="Upper Diamond Lake "| Site=="Forest Lake"| Site=="Lion Lake 2"
             | Site=="Pear Reservoir" | Site=="Yankee Doodle " | Site=="GL1" 
             | Site=="Long Lake" | Site=="Red Deer Lake" | Site=="Mitchell Lake" | Site=="Red Rock Lake", 
             select= Site:Sort)
summary(cp2$Site)


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



#CP_abio_pannel4 <- grid.arrange(CP_temp_plot, CP_TDN_plot, CP_NO3_plot, CP_chlor_plot, CP_cond_plot, CP_TDP_plot, 
#                                CP_PO4_plot, CP_secchi_plot, CP_rho_plot, CP_DOC_plot, CP_FI_plot, CP_kpar_plot, nrow = 3, ncol=4)

#ggsave("CP_abio_pannel4.pdf", CP_abio_pannel4, scale = 2, width = 17, height = 10, units = c("cm"), dpi = 300)












##############################################
########### panel of bio figures #############
CP_chlor_plot <- ggplot(w2, aes(x = Elevation, y = log10(Chl.a +1))) + ylab("Transformed Chl-a") +
  geom_point(aes(colour = as.factor(epi), shape=as.factor(fish))) +
  theme_classic() + theme(text = element_text(size=14), 
                          axis.text.x = element_text(angle=45, hjust=1.2),
                          axis.title.x=element_blank(),
                          legend.position="none", 
                          plot.margin = unit(c(0.15, 0.1,0.15,1.2), "cm")) + 
  scale_x_continuous(breaks=seq(2400,3650,200)) + 
  scale_y_continuous(limits=c(0,1.5), breaks=seq(0.000,1.5, 0.25)) + 
  scale_colour_manual(values = alpha(c("dodgerblue4", "sky blue"), 0.85)) +
  scale_shape_manual(values =c(15, 17))

#YOU ARE HERE #

range(na.omit(log10(d2$phyto_den)))
CP_phyden_plot <- ggplot(d, aes(x = Elevation, y = log10(phyto_den))) + ylab("Transformed density") +
  geom_point(aes(colour = as.factor(epi), shape=as.factor(Fish))) + geom_smooth(method = "glm", colour = 'grey35', fill='grey85') +
  theme_classic() + theme(text = element_text(size=14), 
                          axis.text.x = element_text(angle=45, hjust=1.2),
                          axis.title.x=element_blank(),
                          legend.position="none", 
                          plot.margin = unit(c(0.15, 0.1,0.15,1.2), "cm")) + 
  scale_x_continuous(breaks=seq(2400,3650,200)) + 
  scale_y_continuous(limits=c(1.5,5),breaks=seq(0,5.0,0.55)) + 
  scale_colour_manual(values = alpha(c("dodgerblue4", "sky blue"),0.85)) +
  scale_shape_manual(values = c(15, 17))



range(na.omit(log10(cp$zoop_den2 +1)))
CP_zooden_plot <- ggplot(cp2, aes(x = Elevation, y = log10(zoop_den2 +1))) + ylab("Transformed density") +
  geom_point(aes(shape=as.factor(Fish)), alpha=0.75) + geom_smooth(method = "glm", colour = 'grey35', fill='grey85') +
  theme_classic() + theme(text = element_text(size=14), 
                          axis.text.x = element_text(angle=45, hjust=1.2),
                          axis.title.x=element_blank(),
                          legend.position="none", 
                          plot.margin = unit(c(0.15, 0.1,0.15,1.2), "cm")) + 
  scale_x_continuous(breaks=seq(2400,3650,200)) + 
  scale_y_continuous(limits=c(0,2.25),breaks=seq(0,2.25,0.45))+
  scale_shape_manual(values = c(15, 17))


range(d2$NDVI_QAQC)
CP_NDVI_plot <- ggplot(d, aes(x = Elevation, y = NDVI)) + ylab("NDVI") +
  geom_point(aes(shape=as.factor(Fish)), alpha=0.75)  + geom_smooth(method = "glm", colour = 'grey35', fill='grey85') +
  theme_classic() + theme(text = element_text(size=14), 
                          axis.text.x = element_text(angle=45, hjust=1.2),
                          axis.title.x=element_blank(),
                          legend.position="none", 
                          plot.margin = unit(c(0.15, 0.1,0.15,1.2), "cm")) + 
  scale_x_continuous(breaks=seq(2400,3650,200)) + 
  scale_y_continuous(limits=c(0.12,0.75), breaks=seq(0,0.75, 0.15)) +
  scale_shape_manual(values = c(15, 17))




range(na.omit(cp$Ave_size))
CP_avesize_plot <- ggplot(cp2, aes(x = Elevation, y = Ave_size)) + ylab("Average size") +
  geom_point(aes(shape=as.factor(Fish)), alpha=0.75) + geom_smooth(method = "glm", colour = 'grey35', fill='grey85') +
  theme_classic() + theme(text = element_text(size=14), 
                          axis.text.x = element_text(angle=45, hjust=1.2),
                          axis.title.x=element_blank(),
                          legend.position="none", 
                          plot.margin = unit(c(0.15, 0.1,0.15,1.2), "cm")) + 
  scale_x_continuous(breaks=seq(2400,3650,200)) + 
  scale_y_continuous(limits=c(0.00,2.26),breaks=seq(0.00,2.25,0.45)) + 
  scale_colour_manual(values = c("dodgerblue4", "sky blue")) +
  scale_shape_manual(values = c(15, 17))

cp$prop_gravid
CP_fecund_plot <- ggplot(cp2, aes(x = Elevation, y = prop_gravid)) + ylab("Proportion gravid") +
  geom_point(aes(shape=as.factor(Fish)), alpha=0.75) +
  theme_classic() + theme(text = element_text(size=14), 
                          axis.text.x = element_text(angle=45, hjust=1.2),
                          axis.title.x=element_blank(),
                          legend.position="none", 
                          plot.margin = unit(c(0.15, 0.1,0.15,1.2), "cm")) + 
  scale_x_continuous(breaks=seq(2400,3650,200)) + 
  scale_y_continuous(limits=c(0,1.0),breaks=seq(0,1,0.25)) + 
  scale_colour_manual(values = c("dodgerblue4", "sky blue")) +
  scale_shape_manual(values = c(15, 17))


CP_bio_pannel <- grid.arrange(CP_chlor_plot, CP_NDVI_plot, CP_phyden_plot, CP_avesize_plot, CP_zooden_plot, CP_fecund_plot, nrow = 3, ncol=2)

ggsave("CP_bio_panel3.pdf", CP_bio_pannel, scale = 1.95, width = 89, height = 100, units = c("mm"), dpi = 300)













##########
# coef_plot.csv

co <- read.csv("coef_plot2.csv", header=T)
names(co)
summary(co)

abio_coefplot <- ggplot(co, aes(y = response, x = elev_coef)) + xlab("Coefficient estimate") + ylab("") + 
  geom_vline(xintercept=0, color = "grey35", size=0.75) +
  geom_errorbarh(aes(xmin=(elev_coef-st_error), xmax=(elev_coef+st_error)), height = 0, size= 1) + 
  geom_point(aes(shape=as.character.factor(Sig), color=as.character.factor(Sig)), size=4) +
  theme_linedraw() + theme(axis.title.x = element_text(size=19, margin = margin(t = 10, r = 20, b = 0, l = 0)), 
                     axis.text.x  = element_text(vjust=0.5, size=16), 
                     axis.text.y  = element_text(vjust=0.5, size=16),
                          text = element_text(size=16), 
                         plot.margin = unit(c(0.15, 0.15,0.15,0.15), "cm"),
                     legend.justification=c(-0.1,-0.1), legend.position=c(0,0.65),
                     legend.background = element_rect(size=.5, linetype="dotted"))  + 
  scale_x_continuous(limits=c(-0.4,0.309), breaks = c(-0.4, -0.3, -0.2, -0.1, 0, 0.1, 0.2)) +
  scale_colour_manual(values = c("dodgerblue2", "tomato3", "palegreen4"))

ggsave("abio_coefplot.pdf", abio_coefplot, scale = 5, width = 38, height = 30, units = c("mm"), dpi = 1500)

names(d2)
cp$Ave_size
d3 <- subset(d2, Elevation > 3450,
             select= Site:phyto_den)
summary(d3$Elevation)
mean(na.omit(d3$FI))

d4 <- subset(d2, Elevation < 3000,
             select= Site:phyto_den)
summary(na.omit(d4$FI))

mean(na.omit(d4$FI))/(mean(na.omit(d3$Ave_size)))
1-0.5032808

range(na.omit(d2$FI))

CP_phyden_plot <- ggplot(d2, aes(x = Elevation, y = FI)) + ylab("log10(phyto density)") +
  geom_point(aes(colour = as.factor(epi), shape=as.factor(Fish))) + geom_smooth(method = "glm", colour = 'grey30') +
  theme_classic() + theme(text = element_text(size=14),                          
                          axis.title.x=element_blank(),
                          plot.margin = unit(c(0.15, 0,0.5,1.2), "cm")) + 
  scale_x_continuous(breaks=seq(2400,3600,150))  + 
  scale_colour_manual(values = c("dodgerblue4", "sky blue"))

mean(na.omit(d3$FI))


## distrophic plots
range(d2$DOC_mg_L)
NDVI_chlor_plot <- ggplot(d2, aes(x = NDVI_QAQC, y =chla)) + ylab("Chl-a") + xlab("NDVI") +
  geom_point(aes(colour = as.factor(label), shape=as.factor(Fish))) +
  theme_bw() + theme(text = element_text(size=14),
                          legend.position="none", 
                          plot.margin = unit(c(0.15, 0.1,0.15,1.2), "cm")) + 
  scale_x_continuous(breaks=seq(0,0.75,0.05)) + 
  scale_y_continuous(breaks=seq(0, 35, 5)) + 
  scale_colour_manual(values = c("darkolivegreen4","darkolivegreen3","dodgerblue4", "sky blue")) +
  scale_shape_manual(values = c(15, 17))

DOC_chlor_plot <- ggplot(d2, aes(x =NDVI_QAQC, y =DOC_mg_L)) + ylab("DOC") + xlab("NDVI") +
  geom_point(aes(colour = as.factor(label), shape=as.factor(Fish))) +
  theme_bw() + theme(text = element_text(size=14),
                          legend.position="none", 
                          plot.margin = unit(c(0.15, 0.1,0.15,1.2), "cm")) + 
  scale_x_continuous(breaks=seq(0,0.75,0.05)) + 
  scale_y_continuous(breaks=seq(0, 10, 1.5)) + 
  scale_colour_manual(values = c("darkolivegreen4","darkolivegreen3","dodgerblue4", "sky blue")) +
  scale_shape_manual(values = c(15, 17))

dis_pannel <- grid.arrange(NDVI_chlor_plot, DOC_chlor_plot, nrow = 2, ncol=1)

ggsave("dis_pannel.pdf", dis_pannel, scale = 4.05, width = 30, height = 35, units = c("mm"), dpi = 1500)







co2 <- read.csv("coef_plot_v3.csv", header=T)
names(co2)
summary(co2)

co1 <- read.csv("coeff_plot3.csv", header=T)
names(co1)
summary(co1)


bio_coefplot <- ggplot(co1, aes(y = response, x = elev_coef)) + xlab("Elevation coefficiente") + ylab("") + 
  geom_vline(xintercept=0, color = "grey35", size=0.75) +
  geom_errorbarh(aes(xmin=(elev_coef-st_error), xmax=(elev_coef+st_error)), height = 0, size= 1) + 
  geom_point(aes(shape=as.character.factor(Sig), color=as.character.factor(Sig)), size=2) +
  theme_linedraw() + 
  scale_x_continuous(limits=c(-0.45,0.2), breaks = c(-0.4, -0.3, -0.2, -0.1, 0, 0.1, 0.2, 0.3)) +
  scale_colour_manual(values = c("dodgerblue2", "tomato3", "palegreen4"))

ggsave("NSFbio_coefplot.pdf", bio_coefplot, scale = 3.8, width = 45, height = 10, units = c("mm"), dpi = 1500)


abio_coefplot <- ggplot(co2, aes(y = response, x = elev_coef)) + xlab("Elevation coefficiente") + ylab("") + 
  geom_vline(xintercept=0, color = "grey40", size=0.75) +
  geom_errorbarh(aes(xmin=(elev_coef-st_error), xmax=(elev_coef+st_error)), height = 0, size= 0.5) + 
  geom_point(aes(shape=as.character.factor(Sig), color=as.character.factor(Sig)), size=2) +
  theme_linedraw() + 
  scale_x_continuous(limits=c(-0.22,0.259), breaks=seq(-0.2,0.25, 0.05)) +
  scale_colour_manual(values = c("palegreen4", "dodgerblue2", "tomato3")) +
  scale_shape_manual(values=c(19,17,15))

ggsave("CPabio_coefplot.jpeg", abio_coefplot, scale = 2.75, width = 75, height = 35, units = c("mm"), dpi = 1500)



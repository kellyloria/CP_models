# Comp lakes figures

# load needed packages 
library(vegan)
library(sp) #space
library(ggplot2)
library(lme4)
library(lmerTest)# for p-value
library(MuMIn) # forr squared
library(PerformanceAnalytics) # correlation analysis 
library(TTR)
library(gridExtra)
library(Rmisc)


# read in data file 
d <- read.csv("comp_lakes_master_data_2.csv", header=T)
k <- read.csv("COMPLAKES_2016data_a.csv", header=T) # no inlet outlet data
l <- read.csv("2016_Comp_Lakes_nutrients.csv", header=T)
cl <- read.csv("cl_averages.csv", header=T)


# subset sites for D
d$Site
cl$uv
# exclude mud lake because it is shallow seepage lake

gl4 <- d[c(1:6),]
blue <- d[c(7:9),]
snow <- d[c(10:12),]
lion <- d[c(13:15),]
updia <- d[c(16:19),]
alb <- d[c(20:22),]
gl1 <- d[c(23:26),]
isa <- d[c(27:29),]
jas <- d[c(30:31),]
mith <- d[c(32),]
forest <- d[c(33:34),]
yank <- d[c(35:36),]
pear <- d[c(37:39),]
dia <- d[c(40:42),]
mud <- d[c(43:45),]
red <- d[c(46:48),]
lost <- d[c(49:51),]
long <- d[c(52:54),]
redr <- d[c(55),]
d2 <- d[c(1:37, 41:55),]
d2$Site
summary(d2$Site)
k$Site

k2 <- k[c(1:207, 223:567),]
summary(k2$Site)
summary(k2$elevation)
w_sur <- subset(k2, Depth=="0",
                select= Date:visit)
summary(w_sur$Depth)
head(w_sur)

d2 <- d[c(1:37, 41:55),]
summary(d2$Site)

#not clear why I subset so for now
l2 <- l[c(1:37, 41:107, 111:125),]

summary(l2$SAMPLE)

se <- function(x) {sd(x,na.rm=TRUE)/sqrt(length(x))}

?sd



d$Elevation
DOC_chl <- aggregate(LN_DOC_ChlA ~ Elevation, data=d2, FUN=mean) 
DOC_chl$SE <- aggregate(LN_DOC_ChlA ~ Elevation, data=d2, FUN=se)[,2]

d$DOC_ave
DOC <- aggregate(DOC_ave ~ Elevation, data=d2, FUN=mean) 
DOC$SE <- aggregate(DOC_ave ~ Elevation, data=d2, FUN=se)[,2]

d$Cl_ave
Cl <- aggregate(Cl_ave ~ Elevation, data=d2, FUN=mean) 
Cl$SE <- aggregate(Cl_ave ~ Elevation, data=d2, FUN=se)[,2]

d$UV_abs_ave
UVab <- aggregate(UV_abs_ave ~ Elevation, data=d2, FUN=mean) 
UVab$SE <- aggregate(UV_abs_ave ~ Elevation, data=d2, FUN=se)[,2]

d$TDN_ave
TDN <- aggregate(TDN_ave ~ Elevation, data=d2, FUN=mean) 
TDN$SE <- aggregate(TDN_ave ~ Elevation, data=d2, FUN=se)[,2]

d$IP_ave
IP <- aggregate(IP_ave ~ Elevation, data=d2, FUN=mean) 
IP$SE <- aggregate(IP_ave ~ Elevation, data=d2, FUN=se)[,2]

d$TDP_ave
TDP <- aggregate(TDP_ave ~ Elevation, data=d2, FUN=mean) 
TDP$SE <- aggregate(TDP_ave ~ Elevation, data=d2, FUN=se)[,2]

d$DOP_ave
DOP <- aggregate(DOP_ave ~ Elevation, data=d2, FUN=mean) 
DOP$SE <- aggregate(DOP_ave ~ Elevation, data=d2, FUN=se)[,2]

d$NO3_ave
NO3 <- aggregate(NO3_ave ~ Elevation, data=d2, FUN=mean) 
NO3$SE <- aggregate(NO3_ave ~ Elevation, data=d2, FUN=se)[,2]
range(na.omit(l2$TDN_mg_L))


d$SO4_ave
SO4 <- aggregate(SO4_ave ~ Elevation, data=d2, FUN=mean) 
SO4$SE <- aggregate(SO4_ave ~ Elevation, data=d2, FUN=se)[,2]

d$chal_ave
chla <- aggregate(chal_ave ~ Elevation, data=d2, FUN=mean) 
chla$SE <- aggregate(chal_ave ~ Elevation, data=d2, FUN=se)[,2]


######## with l2
######### Subset for surface samples 
l$SAMPLE
l$Depth

#remove mud
l_nm <- l[c(1:89, 96:125),]
summary(l_nm$SAMPLE)
l_nm$Depth

#select only surface
l2 <- subset(l_nm, Depth==0.0,
              select=SAMPLE:LN_DOC_Chla)
head(l2)
(l2$Depth)

names(l)

l2$elevation
DOC_chl <- aggregate(LN_DOC_Chla ~ elevation, data=l2, FUN=mean) 
DOC_chl$SE <- aggregate(LN_DOC_Chla ~ elevation, data=l2, FUN=se)[,2]
DOC_chl$treeline[DOC_chl$elevation > 3420] <- "16"
DOC_chl$treeline[DOC_chl$elevation < 3425] <- "17"

l2$DOC_mg_L
DOC <- aggregate(DOC_mg_L ~ elevation, data=l2, FUN=mean) 
DOC$SE <- aggregate(DOC_mg_L ~ elevation, data=l2, FUN=se)[,2]
DOC$treeline[DOC$elevation > 3420] <- "16"
DOC$treeline[DOC$elevation < 3425] <- "17"

l2$CL_mg_L
Cl <- aggregate(CL_mg_L  ~ elevation, data=l2, FUN=mean) 
Cl$SE <- aggregate(CL_mg_L  ~ elevation, data=l2, FUN=se)[,2]
Cl$treeline[Cl$elevation > 3420] <- "16"
Cl$treeline[Cl$elevation < 3425] <- "17"

w_sur$elevation
UVab <- aggregate(Light.at.surface ~ elevation, data=w_sur, FUN=mean) 
UVab$SE <- aggregate(Light.at.surface ~ elevation, data=w_sur, FUN=se)[,2]
UVab$treeline[UVab$elevation > 3420] <- "16"
UVab$treeline[UVab$elevation < 3425] <- "17"

l2$TDN_mg_L
TDN <- aggregate(TDN_mg_L ~ elevation, data=l2, FUN=mean) 
TDN$SE <- aggregate(TDN_mg_L ~ elevation, data=l2, FUN=se)[,2]
TDN$treeline[TDN$elevation > 3420] <- "16"
TDN$treeline[TDN$elevation < 3425] <- "17"

l2$TDP_mg_L
TDP <- aggregate(TDP_mg_L ~ elevation, data=l2, FUN=mean) 
TDP$SE <- aggregate(TDP_mg_L ~ elevation, data=l2, FUN=se)[,2]
TDP$treeline[TDP$elevation > 3420] <- "16"
TDP$treeline[TDP$elevation < 3425] <- "17"

NO3 <- aggregate(NO3_mg_L ~ elevation, data=l2, FUN=mean) 
NO3$SE <- aggregate(NO3_mg_L ~ elevation, data=l2, FUN=se)[,2]
NO3$treeline[NO3$elevation > 3420] <- "16"
NO3$treeline[NO3$elevation < 3425] <- "17"

l2$SO4_mg_L
SO4 <- aggregate(SO4_mg_L ~ elevation, data=l2, FUN=mean) 
SO4$SE <- aggregate(SO4_mg_L ~ elevation, data=l2, FUN=se)[,2]
SO4$treeline[SO4$elevation > 3420] <- "16"
SO4$treeline[SO4$elevation < 3425] <- "17"


l2$chla
chla <- aggregate(chla ~ elevation, data=l2, FUN=mean) 
chla$SE <- aggregate(chla ~ elevation, data=l2, FUN=se)[,2]
chla$treeline[chla$elevation > 3420] <- "16"
chla$treeline[chla$elevation < 3425] <- "17"

l2$PO4_mg_L
PO4 <- aggregate(PO4_mg_L ~ elevation, data=l2, FUN=mean) 
PO4$SE <- aggregate(PO4_mg_L ~ elevation, data=l2, FUN=se)[,2]
PO4$treeline[PO4$elevation > 3420] <- "16"
PO4$treeline[PO4$elevation < 3425] <- "17"

l2$FI
FI <- aggregate(FI ~ elevation, data=l2, FUN=mean) 
FI$SE <- aggregate(FI ~ elevation, data=l2, FUN=se)[,2]
FI$treeline[FI$elevation > 3420] <- "16"
FI$treeline[FI$elevation < 3425] <- "17"

Suva <- aggregate(SUVA ~ elevation, data=l2, FUN=mean) 
Suva$SE <- aggregate(SUVA ~ elevation, data=l2, FUN=se)[,2]
Suva$treeline[Suva$elevation > 3420] <- "16"
Suva$treeline[Suva$elevation < 3425] <- "17"

#surface
Temp <- aggregate(Temp ~ elevation, data=w_sur, FUN=mean) 
Temp$SE <- aggregate(Temp ~ elevation, data=w_sur, FUN=se)[,2]
Temp$treeline[Temp$elevation > 3420] <- "16"
Temp$treeline[Temp$elevation < 3425] <- "17"

#surface
pH <- aggregate(pH ~ elevation, data=w_sur, FUN=mean) 
pH$SE <- aggregate(pH ~ elevation, data=w_sur, FUN=se)[,2]
pH$treeline[pH$elevation > 3420] <- "16"
pH$treeline[pH$elevation < 3425] <- "17"


d2$attentuation
atten <- aggregate(attentuation ~ Elevation, data=d2, FUN=mean) 
atten$SE <- aggregate(attentuation ~ Elevation, data=d2, FUN=se)[,2]
atten$treeline[atten$Elevation > 3420] <- "16"
atten$treeline[atten$Elevation < 3425] <- "17"


k2$elevation
Conduct <- aggregate(log10(Conductivity) ~ elevation, data=k2, FUN=mean) 
Conduct$SE <- aggregate(log10(Conductivity) ~ elevation, data=k2, FUN=se)[,2]
Conduct$treeline[Conduct$elevation > 3420] <- "16"
Conduct$treeline[Conduct$elevation < 3425] <- "17"

#surface mwasures
Conduct <- aggregate(Conductivity ~ elevation, data=w_sur, FUN=mean) 
Conduct$SE <- aggregate(Conductivity ~ elevation, data=w_sur, FUN=se)[,2]
Conduct$treeline[Conduct$elevation > 3420] <- "16"
Conduct$treeline[Conduct$elevation < 3425] <- "17"

st_cond <- aggregate(St.C ~ elevation, data=w_sur, FUN=mean) 
st_cond$SE <- aggregate(St.C ~ elevation, data=w_sur, FUN=se)[,2]
st_cond$treeline[st_cond$elevation > 3420] <- "16"
st_cond$treeline[st_cond$elevation < 3425] <- "17"

sur_PAR <- aggregate(Light.at.surface ~ elevation, data=w_sur, FUN=mean) 
sur_PAR$SE <- aggregate(Light.at.surface ~ elevation, data=w_sur, FUN=se)[,2]
sur_PAR$treeline[sur_PAR$elevation > 3420] <- "16"
sur_PAR$treeline[sur_PAR$elevation < 3425] <- "17"

d2$rho_change.hypo_sur.
rho <- aggregate(rho_change.hypo_sur. ~ Elevation, data=d2, FUN=mean) 
rho$SE <- aggregate(rho_change.hypo_sur. ~ Elevation, data=d2, FUN=se)[,2]
rho$treeline[rho$Elevation > 3420] <- "16"
rho$treeline[rho$Elevation < 3425] <- "17"

d2$Secchi
Secchi <- aggregate(Secchi ~ Elevation, data=d2, FUN=mean) 
Secchi$SE <- aggregate(Secchi ~ Elevation, data=d2, FUN=se)[,2]
Secchi$treeline[Secchi$Elevation > 3420] <- "16"
Secchi$treeline[Secchi$Elevation < 3425] <- "17"

d2$zoop_rich
zoop.rich <- aggregate(zoop_rich ~ Elevation, data=d2, FUN=mean) 
zoop.rich$SE <- aggregate(zoop_rich ~ Elevation, data=d2, FUN=se)[,2]
zoop.rich$treeline[zoop.rich$Elevation > 3420] <- "16"
zoop.rich$treeline[zoop.rich$Elevation < 3425] <- "17"

d2$phy_gen_rich_ave
p.rich <- aggregate(phy_gen_rich_ave ~ Elevation, data=d2, FUN=mean) 
p.rich$SE <- aggregate(phy_gen_rich_ave ~ Elevation, data=d2, FUN=se)[,2]
p.rich$treeline[p.rich$Elevation > 3420] <- "16"
p.rich$treeline[p.rich$Elevation < 3425] <- "17"

d2$total_rich
t.rich <- aggregate(total_rich ~ Elevation, data=d2, FUN=mean) 
t.rich$SE <- aggregate(total_rich ~ Elevation, data=d2, FUN=se)[,2]
t.rich$treeline[t.rich$Elevation > 3420] <- "16"
t.rich$treeline[t.rich$Elevation < 3425] <- "17"




summary(mith$Site)

# Water quality
# Temp, Conductivity, ph, rho


#png("Comp lakes water quality.png",  
#    width = 1,
#    height = 2.25,
#    units = "in",
#    res = 1200,
#    pointsize = 1.5 )


png("Comp_lakes_sur_temp.png",  
    width = 2,
    height = 1.2,
    units = "in",
    res = 1200,
    pointsize = 2 )

pch.list <- as.numeric(Temp$treeline) #Using the 2nd example
pch.list

par(mar=c(6, 7, 3, 2),mgp=c(4.5,1.4,0))
plot(x= Temp$elevation, y=Temp$Temp, xlab="Elevation (m)", ylab="Surface temperature", 
     ylim=c(2.5, 15), pch=c(pch.list), cex=2.5, cex.axis=2, cex.lab=2)
#lablist.x<-as.vector(c(2983:3550))
#axis(1, at=seq(2950, 3550, by=100), labels = F)
arrows(Temp$elevation, (Temp$Temp + Temp$SE), 
       Temp$elevation, (Temp$Temp - Temp$SE),
       length=0, code=3, lwd=0.25)
abline(glm(Temp$Temp ~ Temp$elevation), col="gray 40", lty=6, lwd=0.35)

dev.off()

png("Comp_lakes_sur_conduct.png",  
    width = 2,
    height = 1.2,
    units = "in",
    res = 1200,
    pointsize = 2 )

pch.list1 <- as.numeric(Conduct$treeline)
pch.list1

par(mar=c(6, 7, 3, 2),mgp=c(4.5,1.4,0))
plot(x= Conduct$elevation, y=Conduct$Conductivity, xlab="Elevation (m)", ylab="Conductivity", 
     ylim=c(5, 39), pch=c(pch.list), cex=2.5, cex.axis=2, cex.lab=2)
arrows(Conduct$elevation, (Conduct$Conductivity + Conduct$SE), 
       Conduct$elevation, (Conduct$Conductivity - Conduct$SE),
       length=0, code=3, lwd=0.25)
abline(glm(Conduct$Conductivity ~ Conduct$elevation), col="gray 40", lty=6, lwd=0.35)

dev.off()


png("Comp_lakes_sur_pH.png",  
    width = 2,
    height = 1.2,
    units = "in",
    res = 1200,
    pointsize = 2)

pch.list2 <- as.numeric(pH$treeline)
pch.list2

par(mar=c(6, 7, 3, 2),mgp=c(4.5,1.4,0))
plot(x= pH$elevation, y=pH$pH, xlab="Elevation (m)", ylab="pH", 
     ylim=c(6, 9), pch=c(pch.list), cex=2.5, cex.axis=2, cex.lab=2)
arrows(pH$elevation, (pH$pH + pH$SE), 
       pH$elevation, (pH$pH - pH$SE),
       length=0, code=3, lwd=0.25)
#abline(glm(pH$pH ~ pH$elevation), col="gray 40", lty=2)
dev.off()


png("Comp_lakes_rho.png",  
    width = 2,
    height = 1.2,
    units = "in",
    res = 1200,
    pointsize = 2)

pch.list3 <- as.numeric(rho$treeline)
pch.list3

par(mar=c(6, 7, 3, 2),mgp=c(4.5,1.4,0))
plot(x=rho$Elevation, y=rho$rho_change.hypo_sur., xlab="Elevation (m)", ylab="Δρ", 
     ylim=c(0, 0.7), pch=c(pch.list3), cex=2.5, cex.axis=2, cex.lab=2)
arrows(rho$Elevation, (rho$rho_change.hypo_sur. + rho$SE), 
       rho$Elevation, (rho$rho_change.hypo_sur. - rho$SE),
       length=0, code=3, lwd=0.25)
#abline(glm(rho$`log10(rho_change.hypo_sur. + 1)` ~ rho$Elevation), col="gray 40", lty=2)

dev.off()


# Water clarity
# Temp, Conductivity, ph, rho

#png("Comp lakes water clarity.png",  
#    width = 1,
#    height = 2.25,
#    units = "in",
#    res = 1200,
#    pointsize = 1.5 )
#nf <- layout(matrix(c(1, 2, 3, 4, 5), nrow = 4, ncol = 1))
#layout.show(nf)

png("Comp_lakes_secchi.png",  
    width = 2,
    height = 1.2,
    units = "in",
    res = 1200,
    pointsize = 2)

pch.list <- as.numeric(Secchi$treeline)
pch.list

par(mar=c(6, 7, 3, 2),mgp=c(4.5,1.4,0))
plot(x= Secchi$Elevation, y=Secchi$Secchi, xlab="Elevation (m)", ylab="Secchi (m)", 
     ylim=c(0, 7), pch=c(pch.list), cex=2.5, cex.axis=2, cex.lab=2)
arrows(Secchi$Elevation, (Secchi$Secchi + Secchi$SE), 
       Secchi$Elevation, (Secchi$Secchi - Secchi$SE),
       length=0, code=3, lwd=0.25)
abline(glm(Secchi$Secchi ~ Secchi$Elevation), col="gray 40", lty=6, lwd=0.35)

dev.off()

png("Comp_lakes_surPAR.png",  
    width = 2,
    height = 1.2,
    units = "in",
    res = 1200,
    pointsize = 2)

pch.list1 <- as.numeric(sur_PAR$treeline)
pch.list1

par(mar=c(6, 7, 3, 2),mgp=c(4.5,1.4,0))
plot(x= sur_PAR$elevation, y=sur_PAR$Light.at.surface, xlab="Elevation (m)", ylab="Sur PAR", 
     ylim=c(0.25, 1), pch=c(pch.list1), cex=2.5, cex.axis=2, cex.lab=2)
arrows(sur_PAR$elevation, (sur_PAR$Light.at.surface + sur_PAR$SE), 
       sur_PAR$elevation, (sur_PAR$Light.at.surface - sur_PAR$SE),
       length=0, code=3, lwd=0.25)
abline(glm(sur_PAR$Light.at.surface ~ sur_PAR$elevation), col="gray 40", lty=6, lwd=0.35)

dev.off()


png("Comp_lakes_atten.png",  
    width = 2,
    height = 1.2,
    units = "in",
    res = 1200,
    pointsize = 2)
pch.list2 <- as.numeric(atten$treeline)
pch.list2

par(mar=c(6, 7, 3, 2),mgp=c(4.5,1.4,0))
plot(x= atten$Elevation, y=atten$attentuation, xlab="Elevation (m)", ylab="Attenuation", 
     ylim=c(0, 1.75), pch=c(pch.list2), cex=2.5, cex.axis=2, cex.lab=2)
arrows(atten$Elevation, (atten$attentuation + atten$SE), 
       atten$Elevation, (atten$attentuation - atten$SE),
       length=0, code=3, lwd=0.25)
abline(glm(atten$attentuation ~ atten$Elevation), col="gray 40", lty=6, lwd=0.35)
dev.off()

png("Comp_lakes_chlor.png",  
    width = 2,
    height = 1.2,
    units = "in",
    res = 1200,
    pointsize = 2)

pch.list3 <- as.numeric(chla$treeline)
pch.list3
par(mar=c(6, 7, 3, 2),mgp=c(4.5,1.4,0))
plot(x= chla$elevation, y=chla$chla, xlab="Elevation (m)", ylab="Chl-a", 
     ylim=c(0, 17), pch=c(pch.list3), cex=2.5, cex.axis=2, cex.lab=2)
arrows(chla$elevation, (chla$chla + chla$SE), 
       chla$elevation, (chla$chla - chla$SE),
       length=0, code=3, lwd=0.25)
#abline(glm(chla$chla ~ chla$elevation), col="gray 40", lty=2)
dev.off()


png("Comp_lakes_DOC.png",  
    width = 2,
    height = 1.2,
    units = "in",
    res = 1200,
    pointsize = 2)

pch.list4 <- as.numeric(DOC$treeline)
pch.list4
par(mar=c(6, 7, 3, 2),mgp=c(4.5,1.4,0))
plot(x= DOC$elevation, y=DOC$DOC_mg_L, xlab="Elevation (m)", ylab="DOC", 
     ylim=c(0, 6.5), pch=c(pch.list4), cex=2.5, cex.axis=2, cex.lab=2)
arrows(DOC$elevation, (DOC$DOC_mg_L + DOC$SE), 
       DOC$elevation, (DOC$DOC_mg_L - DOC$SE),
       length=0, code=3, lwd=0.25)
abline(glm(DOC$DOC_mg_L ~ DOC$elevation), col="gray 40", lty=6, lwd=0.35)
dev.off()




## Nutrients:
#TDN, No3, TDP, PO4

png("Comp_lakes_TDN.png",  
    width = 2,
    height = 1.2,
    units = "in",
    res = 1200,
    pointsize = 2 )

pch.list <- as.numeric(TDN$treeline)
pch.list

par(mar=c(6, 7, 3, 2),mgp=c(4.5,1.4,0))
plot(x= TDN$elevation, y=TDN$TDN_mg_L, xlab="Elevation (m)", ylab="TDN", 
     ylim=c(0, 0.4), pch=c(pch.list), cex=2.5, cex.axis=2, cex.lab=2)
arrows(TDN$elevation, (TDN$TDN_mg_L + TDN$SE), 
       TDN$elevation, (TDN$TDN_mg_L - TDN$SE),
       length=0, code=3, lwd=0.25)
abline(glm(TDN$TDN_mg_L ~ TDN$elevation), col="gray 40", lty=6, lwd=0.35)

dev.off()


png("Comp_lakes_NO3.png",  
    width = 2,
    height = 1.2,
    units = "in",
    res = 1200,
    pointsize = 2 )
par(mar=c(6, 7, 3, 2),mgp=c(4.5,1.4,0))

pch.list <- as.numeric(NO3$treeline)
pch.list
plot(x= NO3$elevation, y=NO3$NO3_mg_L, xlab="Elevation (m)", ylab="NO3", 
     ylim=c(0, 0.8), pch=c(pch.list), cex=2.5, cex.axis=2, cex.lab=2)
arrows(NO3$elevation, (NO3$NO3_mg_L + NO3$SE), 
       NO3$elevation, (NO3$NO3_mg_L - NO3$SE),
       length=0, code=3, lwd=0.25)
abline(glm(NO3$NO3_mg_L ~ NO3$elevation), col="gray 40", lty=6, lwd=0.35)
dev.off()

png("Comp_lakes_NO3.png",  
    width = 2,
    height = 1.2,
    units = "in",
    res = 1200,
    pointsize = 2 )
par(mar=c(6, 7, 3, 2),mgp=c(4.5,1.4,0))

pch.list <- as.numeric(TDP$treeline)
pch.list

par(mar=c(6, 7, 3, 2),mgp=c(4.5,1.4,0))
plot(x= TDP$elevation, y=TDP$TDP_mg_L, xlab="Elevation (m)", ylab="TDP", 
     ylim=c(0, 0.016), pch=c(pch.list), cex=2.5, cex.axis=2, cex.lab=2)
arrows(TDP$elevation, (TDP$TDP_mg_L + TDP$SE), 
       TDP$elevation, (TDP$TDP_mg_L - TDP$SE),
       length=0, code=3, lwd=0.25)
abline(glm(TDP$TDP_mg_L ~ TDP$elevation), col="gray 40", lty=6, lwd=0.35)

dev.off()


png("Comp_lakes_PO4.png",  
    width = 2,
    height = 1.2,
    units = "in",
    res = 1200,
    pointsize = 2 )

pch.list <- as.numeric(PO4$treeline)
pch.list

par(mar=c(6, 7, 3, 2),mgp=c(4.5,1.4,0))
plot(x= PO4$elevation, y=PO4$PO4_mg_L, xlab="Elevation (m)", ylab="PO4", 
     ylim=c(0.00, 0.003), pch=c(pch.list), cex=2, cex.axis=2, cex.lab=2)
arrows(PO4$elevation, (PO4$PO4_mg_L + PO4$SE), 
       PO4$elevation, (PO4$PO4_mg_L - PO4$SE),
       length=0, code=3, lwd=0.25)
#abline(glm(PO4$PO4_mg_L ~ PO4$elevation), col="gray 40", lty=6, lwd=0.35)

dev.off()





png("Comp_lakes_Cl.png",  
    width = 2,
    height = 1.2,
    units = "in",
    res = 1200,
    pointsize = 2 )

pch.list <- as.numeric(Cl$treeline)
pch.list

par(mar=c(6, 7, 3, 2),mgp=c(4.5,1.4,0))
plot(x= Cl$elevation, y=Cl$CL_mg_L, xlab="Elevation (m)", ylab="Cl", 
     ylim=c(0, 0.5), pch=c(pch.list), cex=2.5, cex.axis=2, cex.lab=2)
arrows(Cl$elevation, (Cl$CL_mg_L + Cl$SE), 
       Cl$elevation, (Cl$CL_mg_L - Cl$SE),
       length=0, code=3, lwd=0.25)
abline(glm(Cl$CL_mg_L ~ Cl$elevation), col="gray 40", lty=6, lwd=0.35)

dev.off()

png("Comp_lakes_SO4.png",  
    width = 2,
    height = 1.2,
    units = "in",
    res = 1200,
    pointsize = 2 )

pch.list <- as.numeric(SO4$treeline)
pch.list

par(mar=c(6, 7, 3, 2),mgp=c(4.5,1.4,0))
plot(x= SO4$elevation, y=SO4$SO4_mg_L, xlab="Elevation (m)", ylab="SO4", 
     ylim=c(0, 8), pch=c(pch.list), cex=2.5, cex.axis=2, cex.lab=2)
arrows(SO4$elevation, (SO4$SO4_mg_L + SO4$SE), 
       SO4$elevation, (SO4$SO4_mg_L - SO4$SE),
       length=0, code=3, lwd=0.25)
#abline(glm(SO4$`log10(SO4_mg_L + 1)` ~ SO4$elevation), col="gray 40", lty=2)




dev.off()


###############

png("Comp_lakes_FI.png",  
    width = 2,
    height = 1.2,
    units = "in",
    res = 1200,
    pointsize = 2)

pch.list <- as.numeric(FI$treeline)
pch.list

par(mar=c(6, 7, 3, 2),mgp=c(4.5,1.4,0))
plot(x= FI$elevation, y=FI$FI, xlab="Elevation (m)", ylab="FI", 
     ylim=c(1.15, 1.6), pch=c(pch.list), cex=2.5, cex.axis=2, cex.lab=2)
arrows(FI$elevation, (FI$FI + FI$SE), 
       FI$elevation, (FI$FI - FI$SE),
       length=0, code=3, lwd=0.25)
#abline(glm(FI$FI ~ FI$elevation), col="gray 40", lty=2)
dev.off()


####### Biotic variables 
png("Comp_lakes_zoop_rich.png",  
    width = 2,
    height = 1.2,
    units = "in",
    res = 1200,
    pointsize = 2)

pch.list <- as.numeric(zoop.rich$treeline)
pch.list

par(mar=c(6, 7, 3, 2),mgp=c(4.5,1.4,0))
plot(x= zoop.rich$Elevation, y=zoop.rich$zoop_rich, xlab="Elevation (m)", ylab="Zooplankton richness", 
     ylim=c(0, 15), pch=c(pch.list), cex=2.5, cex.axis=2, cex.lab=2)
arrows(zoop.rich$Elevation, (zoop.rich$zoop_rich + zoop.rich$SE), 
       zoop.rich$Elevation, (zoop.rich$zoop_rich - zoop.rich$SE),
       length=0, code=3, lwd=0.25)
#abline(glm(zoop.rich$zoop_rich ~ zoop.rich$Elevation), col="gray 40", lty=6, lwd=0.35)

dev.off()


png("Comp_lakes_p_rich.png",  
    width = 2,
    height = 1.2,
    units = "in",
    res = 1200,
    pointsize = 2)

pch.list <- as.numeric(p.rich$treeline)
pch.list

par(mar=c(6, 7, 3, 2),mgp=c(4.5,1.4,0))
p.rich$phy_gen_rich_ave
plot(x= p.rich$Elevation, y=p.rich$phy_gen_rich_ave, xlab="Elevation (m)", ylab="Phytoplankton richness", 
     ylim=c(10, 33), pch=c(pch.list), cex=2.5, cex.axis=2, cex.lab=2)
arrows(p.rich$Elevation, (p.rich$phy_gen_rich_ave + p.rich$SE), 
       p.rich$Elevation, (p.rich$phy_gen_rich_ave - p.rich$SE),
       length=0, code=3,  lwd=0.25)
abline(glm(p.rich$phy_gen_rich_ave ~ p.rich$Elevation), col="gray 40", lty=6, lwd=0.35)

dev.off()

png("Comp_lakes_t_rich.png",  
    width = 2,
    height = 1.2,
    units = "in",
    res = 1200,
    pointsize = 2)

pch.list <- as.numeric(t.rich$treeline)
pch.list
par(mar=c(6, 7, 3, 2),mgp=c(4.5,1.4,0))
t.rich$total_rich
plot(x= t.rich$Elevation, y=t.rich$total_rich, xlab="Elevation (m)", ylab="Phytoplankton richness", 
     ylim=c(20, 60), pch=c(pch.list), cex=2.5, cex.axis=2, cex.lab=2)
arrows(t.rich$Elevation, (t.rich$total_rich + t.rich$SE), 
       t.rich$Elevation, (t.rich$total_rich - t.rich$SE),
       length=0, code=3,  lwd=0.25)
abline(glm(t.rich$total_rich ~ t.rich$Elevation), col="gray 40", lty=6, lwd=0.35)

dev.off()

###############
d2$zoop_den
zoop_den <- aggregate(zoop_den ~ Elevation, data=d2, FUN=mean) 
zoop_den$SE <- aggregate(zoop_den ~ Elevation, data=d2, FUN=se)[,2]
zoop_den$treeline[zoop_den$Elevation > 3420] <- "16"
zoop_den$treeline[zoop_den$Elevation < 3425] <- "17"

d2$BV_zoop
zoop_BV <- aggregate(BV_zoop ~ Elevation, data=d2, FUN=mean) 
zoop_BV$SE <- aggregate(BV_zoop ~ Elevation, data=d2, FUN=se)[,2]
zoop_BV$treeline[zoop_BV$Elevation > 3420] <- "16"
zoop_BV$treeline[zoop_BV$Elevation < 3425] <- "17"

d2$BV_phyto_ave
phyto_BV <- aggregate(BV_phyto_ave ~ Elevation, data=d2, FUN=mean) 
phyto_BV$SE <- aggregate(BV_phyto_ave ~ Elevation, data=d2, FUN=se)[,2]
phyto_BV$treeline[phyto_BV$Elevation > 3420] <- "16"
phyto_BV$treeline[phyto_BV$Elevation < 3425] <- "17"


d2$ratio_ZBV_PBV
z_ph_BV <- aggregate(ratio_ZBV_PBV ~ Elevation, data=d2, FUN=mean) 
z_ph_BV$SE <- aggregate(ratio_ZBV_PBV ~ Elevation, data=d2, FUN=se)[,2]
z_ph_BV$treeline[z_ph_BV$Elevation > 3420] <- "16"
z_ph_BV$treeline[z_ph_BV$Elevation < 3425] <- "17"


png("Comp_lakes_den.png",  
    width = 2,
    height = 1.2,
    units = "in",
    res = 1200,
    pointsize = 2)

pch.list <- as.numeric(zoop_den$treeline)
pch.list

par(mar=c(6, 7, 3, 2),mgp=c(4.5,1.4,0))
plot(x= zoop_den$Elevation, y=zoop_den$zoop_den, xlab="Elevation (m)", ylab="Zooplankton density", 
     ylim=c(0, 90),  pch=c(pch.list), cex=2.5, cex.axis=2, cex.lab=2)
arrows(zoop_den$Elevation, (zoop_den$zoop_den + zoop_den$SE), 
       zoop_den$Elevation, (zoop_den$zoop_den- zoop_den$SE),
       length=0, code=3,  lwd=0.25)
abline(glm(zoop_den$zoop_den ~ zoop_den$Elevation), col="gray 40", lty=6, lwd=0.35)

dev.off()

png("Comp_lakes_BV.png",  
    width = 2,
    height = 1.2,
    units = "in",
    res = 1200,
    pointsize = 2)

pch.list <- as.numeric(zoop_BV$treeline)
pch.list

par(mar=c(6, 7, 3, 2),mgp=c(4.5,1.4,0))
plot(x= zoop_BV$Elevation, y=zoop_BV$BV_zoop, xlab="Elevation (m)", ylab="Zooplankton biovolume", 
     ylim=c(0, 65), pch=c(pch.list), cex=2.5, cex.axis=2, cex.lab=2)
arrows(zoop_BV$Elevation, (zoop_BV$BV_zoop + zoop_BV$SE), 
       zoop_BV$Elevation, (zoop_BV$BV_zoop - zoop_BV$SE),
       length=0, code=3,  lwd=0.25)
#abline(glm(zoop_BV$BV_zoop ~ zoop_BV$Elevation), col="gray 40", lty=6, lwd=0.35)

dev.off()

png("Comp_lakes_BV_phyto.png",  
    width = 2,
    height = 1.2,
    units = "in",
    res = 1200,
    pointsize = 2)

pch.list <- as.numeric(phyto_BV$treeline)
pch.list

par(mar=c(6, 7, 3, 2),mgp=c(4.5,1.4,0))
plot(x= phyto_BV$Elevation, y=phyto_BV$BV_phyto_ave, xlab="Elevation (m)", ylab="Phytoplankton biovolume", 
     ylim=c(250, 5500), pch=c(pch.list), cex=2.5, cex.axis=2, cex.lab=2)
arrows(phyto_BV$Elevation, (phyto_BV$BV_phyto_ave + phyto_BV$SE), 
       phyto_BV$Elevation, (phyto_BV$BV_phyto_ave - phyto_BV$SE),
       length=0, code=3,  lwd=0.25)
abline(glm(phyto_BV$BV_phyto_ave ~ phyto_BV$Elevation), col="gray 40", lty=6, lwd=0.35)

dev.off()

png("Comp_lakes_BV_Z_phyto.png",  
    width = 2,
    height = 1.2,
    units = "in",
    res = 1200,
    pointsize = 2)

pch.list <- as.numeric(z_ph_BV$treeline)
pch.list
par(mar=c(6, 7, 3, 2),mgp=c(4.5,1.4,0))
plot(x= z_ph_BV$Elevation, y=z_ph_BV$ratio_ZBV_PBV, xlab="Elevation (m)", ylab="Zooplankton:Phytoplankton\n biovolume", 
     ylim=c(0, 0.07), pch=c(pch.list), cex=2.5, cex.axis=2, cex.lab=2)
arrows(z_ph_BV$Elevation, (z_ph_BV$ratio_ZBV_PBV + z_ph_BV$SE), 
       z_ph_BV$Elevation, (z_ph_BV$ratio_ZBV_PBV - z_ph_BV$SE),
       length=0, code=3,  lwd=0.25)
#abline(glm(z_ph_BV$ratio_ZBV_PBV ~ z_ph_BV$Elevation), col="gray 40", lty=6, lwd=0.35)

dev.off()
d2$Ave_size

z_size <- aggregate(Ave_size ~ Elevation, data=d2, FUN=mean) 
z_size$SE <- aggregate(Ave_size ~ Elevation, data=d2, FUN=se)[,2]
z_size$treeline[z_size$Elevation > 3420] <- "16"
z_size$treeline[z_size$Elevation < 3425] <- "17"

d2$prop_gravid
zoop_gravid <- aggregate(prop_gravid ~ Elevation, data=d2, FUN=mean) 
zoop_gravid$SE <- aggregate(prop_gravid ~ Elevation, data=d2, FUN=se)[,2]
zoop_gravid$treeline[zoop_gravid$Elevation > 3420] <- "16"
zoop_gravid$treeline[zoop_gravid$Elevation < 3425] <- "17"


png("Comp_lakes_BV_phyto.png",  
    width = 2,
    height = 1.2,
    units = "in",
    res = 1200,
    pointsize = 2)

pch.list <- as.numeric(z_size$treeline)
pch.list

par(mar=c(6, 7, 3, 2),mgp=c(4.5,1.4,0))
plot(x= z_size$Elevation, y=z_size$Ave_size, xlab="Elevation (m)", ylab="Zooplankton average\n size (mm)", 
     ylim=c(0.35, 2.1), pch=c(pch.list), cex=2.5, cex.axis=2, cex.lab=2)
arrows(z_size$Elevation, (z_size$Ave_size + z_size$SE), 
       z_size$Elevation, (z_size$Ave_size - z_size$SE),
       length=0, code=3,  lwd=0.25)
abline(glm(z_size$Ave_size ~ z_size$Elevation), col="gray 40", lty=6, lwd=0.35)

dev.off()

png("Comp_lakes_gravid.png",  
    width = 2,
    height = 1.2,
    units = "in",
    res = 1200,
    pointsize = 2)

pch.list <- as.numeric(zoop_gravid$treeline)
pch.list

par(mar=c(6, 7, 3, 2),mgp=c(4.5,1.4,0))
plot(x= zoop_gravid$Elevation, y=zoop_gravid$prop_gravid, xlab="Elevation (m)", ylab="Proportion of gravid zooplankton", 
     ylim=c(0, 0.52), pch=c(pch.list), cex=2.5, cex.axis=2, cex.lab=2)
arrows(zoop_gravid$Elevation, (zoop_gravid$prop_gravid + zoop_gravid$SE), 
       zoop_gravid$Elevation, (zoop_gravid$prop_gravid - zoop_gravid$SE),
       length=0, code=3,  lwd=0.25)
abline(glm(zoop_gravid$prop_gravid ~ zoop_gravid$Elevation), col="gray 40", lty=6, lwd=0.35)

dev.off()



rel_ab <- read.csv("Comp_lakes_taxa_plot.csv", header=T)
se <- function(x) {sd(x,na.rm=TRUE)/sqrt(length(x))}


new_df <- read.csv("Book2.csv", header=T)

library(gridExtra)

rel_ab$RANaupPerL
cop_RA <- aggregate(RA_COP ~ treeline, data=rel_ab, FUN=mean) 
cop_RA$SE <- aggregate(RA_COP ~ treeline, data=rel_ab, FUN=se)[,2]

clad_RA <- aggregate(RA_CLAD ~ treeline, data=rel_ab, FUN=mean) 
clad_RA$SE <- aggregate(RA_CLAD ~ treeline, data=rel_ab, FUN=se)[,2]

dith_RA <- aggregate(RADithPerL ~ Treeline, data=rel_ab, FUN=mean) 
dith_RA$SE <- aggregate(RADithPerL ~ Treeline, data=rel_ab, FUN=se)[,2]

hesp_RA <- aggregate(RAHespPerL ~ Treeline, data=rel_ab, FUN=mean) 
hesp_RA$SE <- aggregate(RAHespPerL ~ Treeline, data=rel_ab, FUN=se)[,2]

naup_RA <- aggregate(RANaupPerL ~ Treeline, data=rel_ab, FUN=mean) 
naup_RA$SE <- aggregate(RANaupPerL ~ Treeline, data=rel_ab, FUN=se)[,2]

chyd_RA <- aggregate(RAChyPerL ~ Treeline, data=rel_ab, FUN=mean) 
chyd_RA$SE <- aggregate(RAChyPerL ~ Treeline, data=rel_ab, FUN=se)[,2]

holo_RA <- aggregate(RAHoloPerL ~ Treeline, data=rel_ab, FUN=mean) 
holo_RA$SE <- aggregate(RAHoloPerL ~ Treeline, data=rel_ab, FUN=se)[,2]

pul_RA <- aggregate(RAPulPerL ~ Treeline, data=rel_ab, FUN=mean) 
pul_RA$SE <- aggregate(RAPulPerL ~ Treeline, data=rel_ab, FUN=se)[,2]

dia_RA <- aggregate(RA_Diatom ~ Treeline, data=rel_ab, FUN=mean) 
dia_RA$SE <- aggregate(RA_Diatom ~ Treeline, data=rel_ab, FUN=se)[,2]

cyan_RA <- aggregate(RA_Cyan ~ Treeline, data=rel_ab, FUN=mean) 
cyan_RA$SE <- aggregate(RA_Cyan ~ Treeline, data=rel_ab, FUN=se)[,2]

chry_RA <- aggregate(RA_Chry ~ Treeline, data=rel_ab, FUN=mean) 
chry_RA$SE <- aggregate(RA_Chry ~ Treeline, data=rel_ab, FUN=se)[,2]

chlor_RA <- aggregate(RA_chlor ~ Treeline, data=rel_ab, FUN=mean) 
chlor_RA$SE <- aggregate(RA_chlor ~ Treeline, data=rel_ab, FUN=se)[,2]



new_df <- read.csv("Book3.csv", header=T)

png("Comp_lakes_rel_ab.png",  
    width = 5,
    height = 5,
    units = "in",
    res = 1200,
    pointsize = 1.5 )

ggplot(data = new_df, aes(x = treeline, y = Abund, fill=Group)) + 
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  scale_fill_manual(values=c( "coral2", "tomato", "lightsalmon", "paleturquoise3", "cyan4", "dodgerblue3")) + theme_classic() + ylab("Relative abundance") + xlab("Treeline") + geom_errorbar(aes(ymin=Abund, ymax=Abund+SE), width=.2,
                       position=position_dodge(.9))

dev.off()
new_df$SE
p <- ggplot(data = new_df, aes(x = treeline, y = Abund, fill=Group)) + 
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  scale_fill_manual(values=c( "paleturquoise3", "coral3", "dodgerblue3", "coral4", "cyan4", "coral"))
p <- p + theme_classic() + ylab("Relative abundance") + xlab("Treeline")
p <- p + geom_errorbar(aes(ymin=Abund, ymax=Abund+SE), width=.2,
                       position=position_dodge(.9))


rel_ab <- read.csv("Phyto_RA.csv", header=T)

library(gridExtra)

rel_ab$Dinophyta_RA
BAC_RA <- aggregate(Bacillariophyta_RA ~ Treeline, data=rel_ab, FUN=mean) 
BAC_RA$SE <- aggregate(Bacillariophyta_RA ~ Treeline, data=rel_ab, FUN=se)[,2]

chlor_RA <- aggregate(Chlorophyta_RA ~ Treeline, data=rel_ab, FUN=mean) 
chlor_RA$SE <- aggregate(Chlorophyta_RA ~ Treeline, data=rel_ab, FUN=se)[,2]

cyan_RA <- aggregate(Cyanophyta_RA ~ Treeline, data=rel_ab, FUN=mean) 
cyan_RA$SE <- aggregate(Cyanophyta_RA ~ Treeline, data=rel_ab, FUN=se)[,2]

ci_RA <- aggregate(Ciliata_RA ~ Treeline, data=rel_ab, FUN=mean) 
ci_RA$SE <- aggregate(Ciliata_RA ~ Treeline, data=rel_ab, FUN=se)[,2]

cry_RA <- aggregate(Cryptophyta_RA ~ Treeline, data=rel_ab, FUN=mean) 
cry_RA$SE <- aggregate(Cryptophyta_RA ~ Treeline, data=rel_ab, FUN=se)[,2]

chry_RA <- aggregate(Chrysophyta_RA ~ Treeline, data=rel_ab, FUN=mean) 
chry_RA$SE <- aggregate(Chrysophyta_RA ~ Treeline, data=rel_ab, FUN=se)[,2]

eug_RA <- aggregate(Euglenophyta_RA ~ Treeline, data=rel_ab, FUN=mean) 
eug_RA$SE <- aggregate(Euglenophyta_RA ~ Treeline, data=rel_ab, FUN=se)[,2]

dino_RA <- aggregate(Dinophyta_RA ~ Treeline, data=rel_ab, FUN=mean) 
dino_RA$SE <- aggregate(Dinophyta_RA ~ Treeline, data=rel_ab, FUN=se)[,2]

new_df <- read.csv("Book4.csv", header=T)


png("Comp_lakes_rel_ab_phyo.png",  
    width = 5,
    height = 5,
    units = "in",
    res = 1200,
    pointsize = 1.5 )

ggplot(data = new_df, aes(x = treeline, y = Abund, fill=Group)) + 
  geom_bar(stat="identity", position=position_dodge(), colour="black") +
  scale_fill_manual(values=c( "palegoldenrod", "seagreen4","lightcyan", "lightskyblue","limegreen","turquoise", "turquoise3", "steelblue3")) + theme_classic() + ylab("Relative abundance") + xlab("Treeline") + geom_errorbar(aes(ymin=Abund, ymax=Abund+SE), width=.2,
                       position=position_dodge(.9))

dev.off()

rel_ab <- read.csv("Comp_lakes_taxa_plot.csv", header=T)

d2 <- new_df[c(1:2),]
t.test(rel_ab$DithPerL ~ rel_ab$Treeline)
mod.1 <- glm(rel_ab$DithPerL ~ rel_ab$Treeline)
summary(mod.1)
?t.test

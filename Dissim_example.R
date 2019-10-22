#########################
# Compare distance models for zooplankton by building matrices from 2 different .csv:

library(vegan)

# File 1: spatial comparisons for 2 things:
# 1. Lat and Long only = dis_v1
# 2. Elevation distance only = elev_dis3

# lat long distance for visit 1
coords <- d_v[, c("N_coordiate", "W_coordinate", "Site")]
names(coords) <- c("lat", "long", "site")
s.coord <- coordinates(coords[, c("long", "lat")])
dis_v1 <- dist(s.coord[, c("long", "lat")], method = "euclidean") 
max(dis_v1)
min(dis_v1)

# Lat Long and Elevation distance for 1 visit to visualize 
xcood1 <- d_v$N_coordiate
ycood1 <- d_v$W_coordinate
elev1 <- d_v$Elevation
dat2 <- cbind(xcood1, ycood1, elev1)

elev_dis1 <- dist(dat2)
max(elev_dis1)
min(elev_dis1)


###################################
# File 2: Species dissimilarity distances #
#     Data file in a wide form with sp.s as columns 

# can add in species column ID for presence absence data
host_cols <- d[,64:100]

# OR can add in species by name
host_cols <- c("Midge", "Chaoborus..albicus..flavicans", "Bosmina.longirostrus", "Diaptomus.shoshone.","Epischura.sp.", 
               "Alonella.sp.", "Chydorus.sp.", "Microcyclops.sp.", "Daphnia.longispinus", "Daphnia.pulex.pulicaria", 
               "Daphnia.middendorfiana", "Daphnia..unidentifiable.", "Holopedium.gibberum", "Ascomorpha.sp.", "Asplanchna.sp.",
               "Brachionus.sp.", "Collotheca.pelagica", "Conochiloides.Conochilus.spp.", "Euchlanis.sp.", "Filinia.terminalis",
               "Gastropus.sp.", "Keratella.cochlearis", "Keratella.serrulata", "Kellicotia.longispinus", "Monostyla.sp...see.Lecane.",
               "Notholca.squamula", "Notholca.folicea", "Trichocerca.sp.", "Keratella.hiemalis", "Keratella.valga", "Tylotrocha.monopus", 
               "ostrocod", "Diaphanosoma", "Scapholeberis.mucronata", "Macrothrix.spp", "rot.unknown", "Mite")


hjac <- vegdist(d_v[host_cols], binary = T, method="jaccard", upper=T)

# plot the disctance scores of the two matrices 
par(mar=c(5, 5, 4, 2.7))
plot(c(dis_v3), c(hjac), # where dis_v3 = spatial matrix and hjac = species dissim matrix 
     xlab="Distance", ylab="Bray", main="Zoop dissimilarity decay (lat and long)",cex.lab=1.6, cex.axis=1.4, cex.main=1.8)
legend("bottomright", legend=c("Mantel r: 0.264", "Sig: 0.001"), cex=1.4, bty = "n")
abline(lm(hjac ~ dis_v3), lty=3)
mantel(dis_v3, hjac) 

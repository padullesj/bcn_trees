
####
# Title: Script to analyze data on tree species composition for the paper
#        "Socioeconomics predict tree richness, abundance, and composition in Barcelona, Spain"
# Date: 15/02/2022
####

# Load libraries:
library(data.table)
library(plyr)
library(openxlsx)
library(stringr)
library(spdep)
library(dplyr)
library(tidyr)
library(vegan)
library(raster)
library(ggplot2)
library(ggpubr)
library(grid)
library(FSA)
library(grid)
library(ggplotify)
library(cowplot)

###
# Clean up and prepare data:
###

# Load data. It can be downloaded from https://opendata-ajuntament.barcelona.cat/data/en/dataset?tags=Arbres:
parcs<-unique(fread("raw_data/2021_1T_OD_Arbrat_Parcs_BCN.csv", select=c(1,10,19),  data.table = T))
viari<-unique(fread("raw_data/2021_1T_OD_Arbrat_Viari_BCN.csv", select=c(1,10,19),  data.table = T))
zona<-unique(fread("raw_data/2021_1T_OD_Arbrat_Zona_BCN.csv", select=c(1,10,19),  data.table = T))

# Add land-use identifier to the tables:
parcs$type<-"parcs"
viari$type<-"viari"
zona$type<-"zona"
total<-unique(rbind(parcs, viari, zona))

# Fix factor categories:
total$codi_barri<-as.factor(total$codi_barri)
total$codi_barri<-revalue(total$codi_barri, c("01" = "1", "02" = "2", "03"="3", "04"="4", "05" = "5",
                              "06" = "6", "07" = "7", "08"= "8", "09"="9"))

# Standardize taxonomy:
total$nom_cientific<-gsub(" x ", " ", total$nom_cientific)
total$nom_cientific<-gsub(" × ", " ", total$nom_cientific)
total$nom_cientific<-word(total$nom_cientific, 1, 2)
total<-subset(total, nom_cientific !="Indeterminat Arbre")
total$nom_cientific <-gsub("\'", "", total$nom_cientific)
total$nom_cientific[total$nom_cientific=="Prunus Accolade"]<-"Prunus hybrid"
total$nom_cientific[total$nom_cientific=="Malus Evereste"]<-"Malus hybrid"
total$nom_cientific[total$nom_cientific=="Ulmus Dodoens"]<-"Ulmus hybrid"
total$nom_cientific[total$nom_cientific=="Ulmus New"]<-"Ulmus hybrid"
total$nom_cientific[total$nom_cientific=="Ulmus Sapporo"]<-"Ulmus hybrid"
total<-unique(total)
total$codi<-NULL

# Get community data for each land-use type:
cm<-subset(total, type=="parcs") # For parks
com_par<-as.data.frame(dcast(cm, codi_barri ~ nom_cientific)) # Convert table into wide format
rownames(com_par)<-com_par$codi_barri # Assign rownames
com_par$codi_barri<-NULL # Delete column with the codes of the neighborhoods
ncol(com_par) # Get number of species

cm<-subset(total, type=="viari") # For streets
com_via<-as.data.frame(dcast(cm, codi_barri ~ nom_cientific))
com_via<-com_via[-c(74),] 
rownames(com_via)<-com_via$codi_barri
com_via$codi_barri<-NULL
ncol(com_via)

cm<-subset(total, type=="zona") # For zonal areas
com_zon<-as.data.frame(dcast(cm, codi_barri ~ nom_cientific))
com_zon<-com_zon[-c(74),]
rownames(com_zon)<-com_zon$codi_barri
com_zon$codi_barri<-NULL
ncol(com_zon)

com<-as.data.frame(dcast(total, codi_barri ~ nom_cientific)) # Across land-use types
com<-com[-c(74),]
rownames(com)<-com$codi_barri
com$codi_barri<-NULL
ncol(com)

# Get relative abundances:
parcs<-as.data.frame(prop.table(as.matrix(com_par), margin = 1))
viari<-as.data.frame(prop.table(as.matrix(com_via), margin = 1))
zona<-as.data.frame(prop.table(as.matrix(com_zon), margin = 1))
total<-as.data.frame(prop.table(as.matrix(com), margin = 1))

# Load predictors:
bcn<-read.table("processed_data/bcn_div.csv")
bcn$log_RDLpc_2019<-log(bcn$RDLpc_2019)

# Get shapefile of neighborhoods in Barcelona:
nt <- st_read("raw_data/carto/unitats/SHP/unitats/0301040100_Barris_UNITATS_ADM.shp")

# Set theme options:
theme_opts<-list(theme(panel.background = element_blank(),
                       plot.background = element_rect(fill="white"),
                       panel.border = element_blank(),
                       panel.grid.major = element_blank(),
                       axis.text.x = element_blank(),
                       axis.text.y = element_blank(),
                       axis.ticks = element_blank(),
                       axis.title.x = element_blank(),
                       axis.title.y = element_blank(),
                       legend.background=element_blank(),
                       legend.position="none",
                       legend.text=element_text(size=8),
                       legend.title=element_blank(),
                       plot.title = element_text(size=12, face="bold", hjust = 0.5)))

###
# Measure tree species composition (beta-diversity)
###

# Parks

# Perform ordination:
parcs.dist<-vegdist(parcs, "bray")
ord <- cmdscale(parcs.dist)

# Keep dissimilarities for later:
pk<-as.matrix(parcs.dist)
pk1<-pk[lower.tri(pk, diag=F)]

# Subset values for "envfit":
bcn2 <- bcn[bcn$codi_barri %in% rownames(parcs), ]
bcn2 <- subset(bcn2, select=c(log_RDLpc_2019, densitat, vida_esp, edat))
fit <- envfit(ord, bcn2, perm = 999)
fit

# Transform ordination into table:
ord.tb<-as.data.frame(ord)

# Prepare legend of colors:

# Get minimum and maximum values:
x_min <- min(ord.tb$V1)-0.01
x_max <- max(ord.tb$V1)+0.01
y_min <- min(ord.tb$V2)-0.01
y_max <- max(ord.tb$V2)+0.01

# Create raster:
x <- raster(xmn=x_min, xmx=x_max, ymn=y_min, ymx=y_max, res=0.01)
values(x) <- 1:ncell(x)

# Create the 3 palettes:
pal1 <- colorRampPalette(c("black", "cyan"))(ncol(x))
pal2 <- colorRampPalette(c("black", "magenta"))(nrow(x))
pal3 <- colorRampPalette(c("cyan", "yellow"))(nrow(x))

# Create palettes with intermediate colors:
for (i in 2:nrow(x)) {
  pal <- colorRampPalette(c(pal2[i], pal3[i]))(ncol(x))
  pal1<- c(pal1, pal)
}
#plot(x, col=pal1) # Plots raster
cols<-pal1

# Extract colors:
ord.tb$comb <- raster::extract(x, ord.tb[c(1,2)])
a <- sort(ord.tb$comb)
ord.tb$colors<-cols[c(a)]

# Set vectors:
vectors<-as.data.frame(scores(fit, "vectors"))
as.factor(ord.tb$comb)

# Plot ordination:
pp<-ggplot() +
  geom_point(data=ord.tb, aes(x=V1, y=V2, color=as.factor(comb)), size=2) + # Show dots
  scale_color_manual(values = ord.tb$colors)+
  labs(x="PCoA 1", y="PCoA 2") +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        strip.text = element_text(size = 15),
        legend.position="none",
        panel.border = element_rect(fill=NA, colour = "black", size=1),
        legend.title = element_blank())

# Prepare data to plot map:
nc<-nt
nc$BARRI<-revalue(nc$BARRI, c("01" = "1", "02" = "2", "03"="3", "04"="4", "05" = "5",
                                            "06" = "6", "07" = "7", "08"= "8", "09"="9"))
ord.tb$BARRI<-rownames(ord.tb)
nc<-merge(nc, ord.tb, by="BARRI", all.x=T)
nc<-nc[order(as.numeric(nc$BARRI)),]

# Plot map:
ppm<-ggplot(nc) + 
  geom_sf(aes(fill = as.factor(comb)), size=0.25)+
  scale_fill_manual(values=c(ord.tb$colors), na.value = "white")+
  theme_opts

# Merge plots:
ppc<-ggarrange(pp, ppm, labels=c("a)", "b)"), ncol = 2, nrow = 1)


# Streets

# Perform ordination:
parcs.dist<-vegdist(viari, "bray")
ord <- cmdscale(parcs.dist)

# Keep dissimilarities for later:
pk<-as.matrix(parcs.dist)
pk2<-pk[lower.tri(pk, diag=F)]

# Subset values for "envfit":
bcn2 <- bcn[bcn$codi_barri %in% rownames(viari), ]
bcn2 <- subset(bcn2, select=c(log_RDLpc_2019, densitat, vida_esp, edat))
fit <- envfit(ord, bcn2, perm = 999)
fit

# Transform ordination into table:
ord.tb<-as.data.frame(ord)

# Prepare legend of colors:

# Get minimum and maximum values:
x_min <- min(ord.tb$V1)-0.01
x_max <- max(ord.tb$V1)+0.01
y_min <- min(ord.tb$V2)-0.01
y_max <- max(ord.tb$V2)+0.01

# Create raster:
x <- raster(xmn=x_min, xmx=x_max, ymn=y_min, ymx=y_max, res=0.01)
values(x) <- 1:ncell(x)

# Create the 3 palettes:
pal1 <- colorRampPalette(c("black", "cyan"))(ncol(x))
pal2 <- colorRampPalette(c("black", "magenta"))(nrow(x))
pal3 <- colorRampPalette(c("cyan", "yellow"))(nrow(x))

# Create palettes with intermediate colors:
for (i in 2:nrow(x)) {
  pal <- colorRampPalette(c(pal2[i], pal3[i]))(ncol(x))
  pal1<- c(pal1, pal)
}
cols<-pal1

# Extract colors:
ord.tb$comb <- raster::extract(x, ord.tb[c(1,2)])
ord.tb<-ord.tb[order(ord.tb$comb),]
ord.tb$colors<-cols[c(ord.tb$comb)]

# Set vectors:
vectors<-as.data.frame(scores(fit, "vectors"))
vectors<-vectors[-c(2,4),]

# Plot ordination:
rr<-ggplot() +
  geom_point(data=ord.tb, aes(x=V1, y=V2, color=as.factor(comb)), size=2) + # Show dots
  scale_color_manual(values = ord.tb$colors)+
  geom_segment(data=vectors, aes(x=0, y=0, xend =Dim1, yend = Dim2), color="red", arrow = arrow(length = unit(0.03, "npc"))) +
  labs(x="PCoA 1", y="PCoA 2") +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        strip.text = element_text(size = 15),
        legend.position="none",
        panel.border = element_rect(fill=NA, colour = "black", size=1),
        legend.title = element_blank())

# Prepare data to plot map:
nc<-nt
nc$BARRI<-revalue(nc$BARRI, c("01" = "1", "02" = "2", "03"="3", "04"="4", "05" = "5",
                              "06" = "6", "07" = "7", "08"= "8", "09"="9"))
ord.tb$BARRI<-rownames(ord.tb)
nc<-merge(nc, ord.tb, by="BARRI", all.x=T)
nc<-nc[order(as.numeric(nc$BARRI)),]

# Plot map:
rrm<-ggplot(nc) + 
  geom_sf(aes(fill = as.factor(comb)), size=0.25)+
  scale_fill_manual(values=c(ord.tb$colors), na.value = "white")+
  theme_opts

# Merge plots:
prr<-ggarrange(rr, rrm, labels=c("c)", "d)"), ncol = 2, nrow = 1)


# Zonal

# Perform ordination:
parcs.dist<-vegdist(zona, "bray")
ord <- cmdscale(parcs.dist)

# Keep dissimilarities for later:
pk<-as.matrix(parcs.dist)
pk3<-pk[lower.tri(pk, diag=F)]

# Subset values for "envfit":
bcn2 <- bcn[bcn$codi_barri %in% rownames(zona), ]
bcn2 <- subset(bcn2, select=c(log_RDLpc_2019, densitat, vida_esp, edat))
fit <- envfit(ord, bcn2, perm = 999)
fit

# Transform ordination into table:
ord.tb<-as.data.frame(ord)

# Prepare legend of colors:

# Get minimum and maximum values:
x_min <- min(ord.tb$V1)-0.01
x_max <- max(ord.tb$V1)+0.01
y_min <- min(ord.tb$V2)-0.01
y_max <- max(ord.tb$V2)+0.01

# Create raster:
x <- raster(xmn=x_min, xmx=x_max, ymn=y_min, ymx=y_max, res=0.01)
values(x) <- 1:ncell(x)

# Create the 3 palettes:
pal1 <- colorRampPalette(c("black", "cyan"))(ncol(x))
pal2 <- colorRampPalette(c("black", "magenta"))(nrow(x))
pal3 <- colorRampPalette(c("cyan", "yellow"))(nrow(x))

# Create palettes with intermediate colors:
for (i in 2:nrow(x)) {
  pal <- colorRampPalette(c(pal2[i], pal3[i]))(ncol(x))
  pal1<- c(pal1, pal)
}
cols<-pal1

# Extract colors:
ord.tb$comb <- raster::extract(x, ord.tb[c(1,2)])
ord.tb<-ord.tb[order(ord.tb$comb),]
ord.tb$colors<-cols[c(ord.tb$comb)]
ord.tb$comb[40]<-2350

# Set vectors:
vectors<-as.data.frame(scores(fit, "vectors"))
vectors<-vectors[-c(1,3,4),]

# Plot ordination:
zz<-ggplot() +
  geom_point(data=ord.tb, aes(x=V1, y=V2, color=as.factor(comb)), size=2) + # Show dots
  scale_color_manual(values = ord.tb$colors)+
  geom_segment(data=vectors, aes(x=0, y=0, xend =Dim1, yend = Dim2), color="red", arrow = arrow(length = unit(0.03, "npc"))) +
  labs(x="PCoA 1", y="PCoA 2") +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        strip.text = element_text(size = 15),
        legend.position="none",
        panel.border = element_rect(fill=NA, colour = "black", size=1),
        legend.title = element_blank())

# Prepare data to plot map:
nc<-nt
nc$BARRI<-revalue(nc$BARRI, c("01" = "1", "02" = "2", "03"="3", "04"="4", "05" = "5",
                              "06" = "6", "07" = "7", "08"= "8", "09"="9"))
ord.tb$BARRI<-rownames(ord.tb)
nc<-merge(nc, ord.tb, by="BARRI", all.x=T)
nc<-nc[order(as.numeric(nc$BARRI)),]

# Plot map:
zzm<-ggplot(nc) + 
  geom_sf(aes(fill = as.factor(comb)), size=0.25)+
  scale_fill_manual(values=c(ord.tb$colors), na.value = "white")+
  theme_opts

# Merge plots:
pzz<-ggarrange(zz, zzm, labels=c("e)", "f)"), ncol = 2, nrow = 1)


# Total

# Perform ordination:
parcs.dist<-vegdist(total, "bray")
ord <- cmdscale(parcs.dist)

# Keep dissimilarities for later:
pk<-as.matrix(parcs.dist)
pk4<-pk[lower.tri(pk, diag=F)]

# Subset values for "envfit":
bcn2 <- bcn[bcn$codi_barri %in% rownames(total), ]
bcn2 <- subset(bcn2, select=c(log_RDLpc_2019, densitat, vida_esp, edat))
fit <- envfit(ord, bcn2, perm = 999)
fit

# Transform ordination into table:
ord.tb<-as.data.frame(ord)

# Prepare legend of colors:

# Get minimum and maximum values:
x_min <- min(ord.tb$V1)-0.01
x_max <- max(ord.tb$V1)+0.01
y_min <- min(ord.tb$V2)-0.01
y_max <- max(ord.tb$V2)+0.01

# Create raster:
x <- raster(xmn=x_min, xmx=x_max, ymn=y_min, ymx=y_max, res=0.01)
values(x) <- 1:ncell(x)

# Create the 3 palettes:
pal1 <- colorRampPalette(c("black", "cyan"))(ncol(x))
pal2 <- colorRampPalette(c("black", "magenta"))(nrow(x))
pal3 <- colorRampPalette(c("cyan", "yellow"))(nrow(x))

# Create palettes with intermediate colors:
for (i in 2:nrow(x)) {
  pal <- colorRampPalette(c(pal2[i], pal3[i]))(ncol(x))
  pal1<- c(pal1, pal)
}
cols<-pal1

# Extract colors:
ord.tb$comb <- raster::extract(x, ord.tb[c(1,2)])
ord.tb<-ord.tb[order(ord.tb$comb),]
ord.tb$colors<-cols[c(ord.tb$comb)]
ord.tb$comb[7]<-1283
ord.tb$comb[52]<-3104

# Set vectors:
vectors<-as.data.frame(scores(fit, "vectors"))
vectors<-vectors[-c(1,3,4),]

# Plot ordination:
tt<-ggplot() +
  geom_point(data=ord.tb, aes(x=V1, y=V2, color=as.factor(comb)), size=2) + # Show dots
  scale_color_manual(values = ord.tb$colors)+
  geom_segment(data=vectors, aes(x=0, y=0, xend =Dim1, yend = Dim2), color="red", arrow = arrow(length = unit(0.03, "npc"))) +
  labs(x="PCoA 1", y="PCoA 2") +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        strip.text = element_text(size = 15),
        legend.position="none",
        panel.border = element_rect(fill=NA, colour = "black", size=1),
        legend.title = element_blank())

# Prepare data to plot map:
nc<-nt
nc$BARRI<-revalue(nc$BARRI, c("01" = "1", "02" = "2", "03"="3", "04"="4", "05" = "5",
                              "06" = "6", "07" = "7", "08"= "8", "09"="9"))
ord.tb$BARRI<-rownames(ord.tb)
nc<-merge(nc, ord.tb, by="BARRI", all.x=T)
nc<-nc[order(as.numeric(nc$BARRI)),]

# Plot map:
ttm<-ggplot(nc) + 
  geom_sf(aes(fill = as.factor(BARRI)), size=0.25)+
  scale_fill_manual(values=nc$colors)+
  theme_opts

# Merge plots:
ptt<-ggarrange(tt, ttm, labels=c("g)", "h)"), ncol = 2, nrow = 1)

# Merge plots of the different land-use types:
p<-ggarrange(annotate_figure(ppc, top = text_grob("Parks", face = "bold", size = 14)),  
             annotate_figure(prr, top = text_grob("Streets", face = "bold", size = 14)),  
             annotate_figure(pzz, top = text_grob("Zonal", face = "bold", size = 14)),
             annotate_figure(ptt, top = text_grob("Total", face = "bold", size = 14)), 
             ncol = 1, nrow = 4)


###
# Create color legend
###

# Create raster:
x <- raster(xmn=1, xmx=10, ymn=1, ymx=10, res=.1)
values(x) <- 1:ncell(x)

# Set number of cells:
num<-sqrt(ncell(x))

# Create first the 3 palettes:
pal1 <- colorRampPalette(c("black", "cyan"))(num)
pal2 <- colorRampPalette(c("black", "magenta"))(num)
pal3 <- colorRampPalette(c("cyan", "yellow"))(num)

# Create palettes with intermediate colors:
for (i in 2:num) {
  pal <- colorRampPalette(c(pal2[i], pal3[i]))(num)
  pal1<-c(pal1, pal)
}

# Plot legend
par(mar=c(0, 0, 0, 0), xpd = NA); plot(x, col=pal1, legend=FALSE, axes = 0, frame.plot=0, box = FALSE, useRaster=0)
p1<- ~plot(x, col=pal1, legend=FALSE, axes = 0, frame.plot=0, box = FALSE, useRaster=0)
p2<- as_grob(p1)

# Save final result and add vector names manually:
png("results/Fig5.png",
    res=600, height=10,width=5,units="in"); 

p + annotation_custom(p2, xmin  = 0.33, xmax = 0.82, ymin = 0.693, ymax = .92)

grid.text("HDIpc (log)", x = unit(0.173, "npc"), y = unit(0.555, "npc"), 
          gp=gpar(fontsize=8, fontface="bold", col="red"))
grid.text("Life\nexpectancy", x = unit(0.36, "npc"), y = unit(0.56, "npc"), 
          gp=gpar(fontsize=8, fontface="bold", col="red"))

grid.text("Population\ndensity", x = unit(0.181, "npc"), y = unit(0.412, "npc"), 
          gp=gpar(fontsize=8, fontface="bold", col="red"))

grid.text("Population\ndensity", x = unit(0.39, "npc"), y = unit(0.06, "npc"), 
          gp=gpar(fontsize=8, fontface="bold", col="red"))

grid.text("Color\nlegend", x = unit(0.54, "npc"), y = unit(0.86, "npc"), 
          gp=gpar(fontsize=9, fontface="bold"))

dev.off()


###
# Violin plots for tree composition
###

# Get mean beta-diversity:
ct<-data.frame(value=c(pk1, pk2, pk3, pk4),
               cat=c(rep("Parks", length(pk1)), rep("Streets", length(pk2)), rep("Zonal", length(pk3)), rep("Total", length(pk4))))

# Set factor levels:
ct$cat<- factor(ct$cat, levels = c("Parks", "Streets", "Zonal", "Total"))

# Kruskal-Wallis test and posthoc comparisons:
kruskal.test(value ~ cat, data = ct)
dunnTest(value~cat,data=ct)

# Plot:
dp2 <- ggplot(ct, aes(x=cat, y=value, fill=cat)) + 
  geom_violin(trim=FALSE)+
  geom_boxplot(width=0.1, fill="white")+
  labs(x="", y = "Bray Curtis distances")+
  geom_vline(xintercept = 3.5, linetype = "longdash") +
  theme_classic()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        strip.text = element_text(size = 15),
        legend.position="none",
        legend.title = element_blank())

# Save result and manually add letters:
png("results/Fig4.png",
    res=600, height=3,width=3.5,units="in"); 
dp2
grid.text("a", x = unit(0.26, "npc"), y = unit(0.97, "npc"), 
          gp=gpar(fontsize=10, fontface="bold"))
grid.text("c", x = unit(0.46, "npc"), y = unit(0.97, "npc"), 
          gp=gpar(fontsize=10, fontface="bold"))
grid.text("b", x = unit(0.66, "npc"), y = unit(0.97, "npc"), 
          gp=gpar(fontsize=10, fontface="bold"))
dev.off()

# Clean-up environment:
rm(list = ls())

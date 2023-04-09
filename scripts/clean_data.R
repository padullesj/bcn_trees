
####
# Title: Script to clean-up the data and calculate tree diverstiy and abundance in the neighorhoods for
#        the paper "Socioeconomics predict tree richness, abundance, and composition in Barcelona, Spain"
# Date: 15/02/2022
####

# Load libraries:
library(data.table)
library(plyr)
library(openxlsx)
library(stringr)
library(spdep)
library(spatialreg)
library(ggplot2)
library(ggpubr)
library(grid)
library(vegan)
library(reshape2)

###
# Clean up data:
###

# Load data. It can be downloaded from https://opendata-ajuntament.barcelona.cat/data/en/dataset?tags=Arbres:
parcs<-unique(fread("raw_data/2021_1T_OD_Arbrat_Parcs_BCN.csv", select=c(1,10,19),  data.table = T))
viari<-unique(fread("raw_data/2021_1T_OD_Arbrat_Viari_BCN.csv", select=c(1,10,19),  data.table = T))
zona<-unique(fread("raw_data/2021_1T_OD_Arbrat_Zona_BCN.csv", select=c(1,10,19),  data.table = T))

# Add land-use identifier to the tables:
parcs$type<-"parcs"
viari$type<-"viari"
zona$type<-"zona"
total<-unique(rbind(parcs, viari, zona)) # Merge tables

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
total$codi<-NULL # Remove individual species identifier

# Get community data for each land-use type:
cm<-subset(total, type=="parcs") # For parks
com_par<-dcast(cm, codi_barri ~ nom_cientific) # Convert table into wide format
rownames(com_par)<-com_par$codi_barri # Assign rownames
com_par$codi_barri<-NULL # Delete column with the codes of the neighborhoods

cm<-subset(total, type=="viari") # For streets
com_via<-dcast(cm, codi_barri ~ nom_cientific)
com_via<-com_via[-c(74),] 
rownames(com_via)<-com_via$codi_barri
com_via$codi_barri<-NULL

cm<-subset(total, type=="zona") # For zonal areas
com_zon<-dcast(cm, codi_barri ~ nom_cientific)
com_zon<-com_zon[-c(74),]
rownames(com_zon)<-com_zon$codi_barri
com_zon$codi_barri<-NULL

com<-dcast(total, codi_barri ~ nom_cientific) # Across land-use types
com<-com[-c(74),]
rownames(com)<-com$codi_barri
com$codi_barri<-NULL

###
# Get diversity indices
###

# Create list with community data:
coms<-list(com_par, com_via, com_zon, com)
cod<-c("parcs", "viari", "zona", "total")

# Loop for each land-use type:
out<-list()
for (i in 1:length(coms)) {
  
  # Calculate tree richness, abundance, and the Shannon index:
  sr<-specnumber(coms[[i]]) # Species richness
  h<-diversity(coms[[i]], index = "shannon") # Shannon index
  nt<-rowSums(coms[[i]]) # Tree abundance

  # Save into data frame:
  out[[i]]<-data.frame(codi_barri=names(sr), sr=sr, h=h, nt=nt)
  
  # Rename columns:
  colnames(out[[i]])<-c("codi_barri", paste("SR_",cod[i],sep=""), 
                        paste("shan_",cod[i],sep=""), paste("NT_",cod[i],sep=""))
}

###
# Combine with demographic and socioeconomic variables:
###

# Load data for demographic and socioeconomic variables:
bcn<-read.csv2("raw_data/bcn_renta2.csv", header=T, sep=";", dec=",")

# Merge:
bcn<-merge(bcn, out[[1]], by="codi_barri", all.x=T)
bcn<-merge(bcn, out[[2]], by="codi_barri", all.x=T)
bcn<-merge(bcn, out[[3]], by="codi_barri", all.x=T)
bcn<-merge(bcn, out[[4]], by="codi_barri", all.x=T)

# Replace cases with no trees as zeros:
bcn$NT_parcs[is.na(bcn$NT_parcs)] <- 0
bcn$SR_parcs[is.na(bcn$SR_parcs)] <- 0

# Save output:
write.table(bcn, "processed_data/bcn_div.csv")

# Clean up environment:
rm(list = ls())

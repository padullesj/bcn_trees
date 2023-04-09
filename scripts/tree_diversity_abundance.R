
####
# Title: Script to analyze data on tree species richness and abundance for the paper
#        "Socioeconomics predict tree richness, abundance, and composition in Barcelona, Spain"
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
library(FSA)

# Load data:
bcn<-read.table("processed_data/bcn_div.csv")

# Count trees and proportions:
sum(bcn$NT_total)
sum(bcn$NT_parcs)/sum(bcn$NT_total)
sum(bcn$NT_viari)/sum(bcn$NT_total)
sum(bcn$NT_zona)/sum(bcn$NT_total)

# Test for correlations between predictors:
cor(bcn[,c(5,7,8,9)], method = c("spearman"))

# Get distances between neighborhoods
geod<-as.matrix(dist(bcn[,c(4,5)]))
geod.inv <- 1/geod
diag(geod.inv) <- 0

###
# Models for tree species richness:
###

# Run GLMs:
smp <- glm(SR_parcs ~ scale(log(RDLpc_2019)) + scale(densitat) + scale(edat) +  scale(vida_esp)  + scale(log(area)), data=bcn, family = quasipoisson(link = "log"))
smv <- glm(SR_viari ~ scale(log(RDLpc_2019)) + scale(densitat) + scale(edat) +  scale(vida_esp)  + scale(log(area)), data=bcn, family = quasipoisson(link = "log"))
smz <- glm(SR_zona ~ scale(log(RDLpc_2019)) + scale(densitat) + scale(edat) +  scale(vida_esp)  + scale(log(area)), data=bcn, family = quasipoisson(link = "log"))
smt <- glm(SR_total ~ scale(log(RDLpc_2019)) + scale(densitat) + scale(edat) +  scale(vida_esp)  + scale(log(area)), data=bcn, family = quasipoisson(link = "log"))

# Get model summaries:
summary(smp)
summary(smv)
summary(smz)
summary(smt)

# Example to check diagnostic plots:
smp.diag <- boot::glm.diag(smp)
boot::glm.diag.plots(smp, smp.diag)

# Tests for spatial autocorrelation:
ape::Moran.I(residuals(smp), geod.inv)
ape::Moran.I(residuals(smv), geod.inv)
ape::Moran.I(residuals(smz), geod.inv)
ape::Moran.I(residuals(smt), geod.inv)

# Get McFadden's R-squared for model
with(summary(smp), 1 - deviance/null.deviance)
with(summary(smv), 1 - deviance/null.deviance)
with(summary(smz), 1 - deviance/null.deviance)
with(summary(smt), 1 - deviance/null.deviance)

###
# Models for tree abundance:
###

# Run GLMs:
tmp <- glm(NT_parcs ~ scale(log(RDLpc_2019)) + scale(densitat) + scale(edat) +  scale(vida_esp)  + scale(log(area)), data=bcn, family = quasipoisson(link = "log"))
tmv <- glm(NT_viari ~ scale(log(RDLpc_2019)) + scale(densitat) + scale(edat) +  scale(vida_esp)  + scale(log(area)), data=bcn, family = quasipoisson(link = "log"))
tmz <- glm(NT_zona ~ scale(log(RDLpc_2019)) + scale(densitat) + scale(edat) +  scale(vida_esp)  + scale(log(area)), data=bcn, family = quasipoisson(link = "log"))
tmt <- glm(NT_total ~ scale(log(RDLpc_2019)) + scale(densitat) + scale(edat) +  scale(vida_esp)  + scale(log(area)), data=bcn, family = quasipoisson(link = "log"))

# Get model summaries:
summary(tmp)
summary(tmv)
summary(tmz)
summary(tmt)

# Test for spatial autocorrelation:
ape::Moran.I(residuals(tmp), geod.inv)
ape::Moran.I(residuals(tmv), geod.inv)
ape::Moran.I(residuals(tmz), geod.inv)
ape::Moran.I(residuals(tmt), geod.inv)

# Get McFadden's R-squared for model
with(summary(tmp), 1 - deviance/null.deviance)
with(summary(tmv), 1 - deviance/null.deviance)
with(summary(tmz), 1 - deviance/null.deviance)
with(summary(tmt), 1 - deviance/null.deviance)

###
# Plot results from the models:
###

# Create dual indices of response variables:
index1<-c("Parks", "Streets", "Zonal", "Total",
         "Parks", "Streets", "Zonal", "Total")

index2<-c("Species richness", "Species richness", "Species richness", "Species richness",
          "Tree abundance", "Tree abundance", "Tree abundance", "Tree abundance")

# Create list with model output results:
models<-list(smp, smv, smz, smt, tmp, tmv, tmz, tmt)

# Loop for each GLM:
out<-NULL
for (i in 1:length(index1)) {
  
  # Retrieve object:
  m<-models[[i]]

  # Prepare summary table:
  rat<-as.data.frame(summary(m)$coefficients) # Extract coefficients
  names(rat)[2]<-"SD" # Get standard deviation
  rat<-rat[-c(1), ] # Remove Intercept
  rat$var<-c("log(HDIpc)", "Population density", "Mean Age", "Life expectancy", "log(Area)") # Assign predictor names
  rat$var <- factor(rat$var, levels = c("log(Area)", "Life expectancy", "Mean Age", "Population density", "log(HDIpc)")) # Change order of factors
  
  # Set confidence interval:
  rat$low<-rat$Estimate - (rat$SD* 1.959)
  rat$up<-rat$Estimate + (rat$SD* 1.959)
  
  # Define groups of variables based on the sign of the effect:
  rat$col <- ifelse((rat$low < 0) & (rat$up < 0),"Negative", ifelse((rat$low > 0) & (rat$up > 0),"Positive", ifelse("Not significant")))
  rat$col[is.na(rat$col)] <- "Not significant"
  
  # Add variable name:
  rat$index1<-index1[i]
  rat$index2<-index2[i]
  
  # Save final table:
  out<-rbind(out, rat)
}

# Set order of the differet factors:
out$index1 <- factor(out$index1, levels = c("Parks", "Streets", "Zonal", "Total")) #changer order of factors
out$index2 <- factor(out$index2, levels = c("Species richness", "Tree abundance")) #changer order of factors
out$col <- factor(out$col, levels = c("Not significant", "Negative", "Positive")) #changer order of factors

# Plot species richness:
sr<-subset(out, index2=="Species richness")
p_sr<-ggplot(sr, aes(x=var, y=Estimate, ymin=low, ymax=up, shape=col, fill=col, col=col)) +  
  ggh4x::facet_nested(. ~ index2 + index1) +
  geom_errorbar(size=0.5, width = 0.2) + 
  geom_point(size=3.5) +
  geom_vline(xintercept=1.5, linetype="dashed") +
  scale_shape_manual(values=c(21, 21, 21))+
  scale_color_manual(values=c("black", "darkorange", "blue")) +
  scale_fill_manual(values=c("white", "darkorange", "blue")) +
  geom_hline(yintercept=0, lty=2, size=0.25) +
  coord_flip() +
  xlab("") + ylab("Estimates (± 95% confident interval)") + 
  theme_bw() +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        legend.title = element_blank(),
        legend.position="top",
        axis.title.x = element_text(size = 10),
        axis.text.x = element_text(angle=45, vjust=.5),
        #panel.spacing.x = unit(4, "mm"),
        plot.title = element_text(color="black", hjust = 0.5, size=8, face="bold"),
        plot.margin = unit(c(0.1,0.5,0.1,0.1), "cm"))

# Plot tree abundance:
ta<-subset(out, index2=="Tree abundance")
p_ta<-ggplot(ta, aes(x=var, y=Estimate, ymin=low, ymax=up, shape=col, fill=col, col=col)) +  
  ggh4x::facet_nested(. ~ index2 + index1) +
  geom_errorbar(size=0.5, width = 0.2) + 
  geom_point(size=3.5) +
  geom_vline(xintercept=1.5, linetype="dashed") +
  scale_shape_manual(values=c(21, 21, 21))+
  scale_color_manual(values=c("black", "darkorange", "blue")) +
  scale_fill_manual(values=c("white", "darkorange", "blue")) +
  geom_hline(yintercept=0, lty=2, size=0.25) +
  coord_flip() +
  xlab("") + ylab("Estimates (± 95% confident interval)") + 
  theme_bw() +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        legend.title = element_blank(),
        legend.position="top",
        axis.title.x = element_text(size = 10),
        axis.text.x = element_text(angle=45, vjust=.5),
        #panel.spacing.x = unit(4, "mm"),
        plot.title = element_text(color="black", hjust = 0.5, size=8, face="bold"),
        plot.margin = unit(c(0.1,0.5,0.1,0.1), "cm"))

# Put plots together:
p<-ggarrange(p_sr, p_ta,
             #labels=c("a)", "b)"),
             common.legend = T,
             ncol = 1, nrow = 2)

# Save result and add pseudo-R squared manually:
png("results/Fig3.png",
    res=600, height=6,width=5.5,units="in"); 
p

grid.text("a)", x = unit(0.05, "npc"), y = unit(0.91, "npc"), 
          gp=gpar(fontsize=11, fontface="bold"))
grid.text("b)", x = unit(0.05, "npc"), y = unit(0.445, "npc"), 
          gp=gpar(fontsize=11, fontface="bold"))


grid.text(expression(R^2*" = 0.20"), x = unit(0.365, "npc"), y = unit(0.838, "npc"), 
          gp=gpar(fontsize=7.5))
grid.text(expression(R^2*" = 0.53"), x = unit(0.553, "npc"), y = unit(0.838, "npc"), 
          gp=gpar(fontsize=7.5))
grid.text(expression(R^2*" = 0.35"), x = unit(0.739, "npc"), y = unit(0.838, "npc"), 
          gp=gpar(fontsize=7.5))
grid.text(expression(R^2*" = 0.36"), x = unit(0.925, "npc"), y = unit(0.838, "npc"), 
          gp=gpar(fontsize=7.5))

grid.text(expression(R^2*" = 0.25"), x = unit(0.365, "npc"), y = unit(0.371, "npc"), 
          gp=gpar(fontsize=7.5))
grid.text(expression(R^2*" = 0.53"), x = unit(0.553, "npc"), y = unit(0.371, "npc"), 
          gp=gpar(fontsize=7.5))
grid.text(expression(R^2*" = 0.43"), x = unit(0.739, "npc"), y = unit(0.371, "npc"), 
          gp=gpar(fontsize=7.5))
grid.text(expression(R^2*" = 0.52"), x = unit(0.925, "npc"), y = unit(0.371, "npc"), 
          gp=gpar(fontsize=7.5))
dev.off()



###
#Violin plots for tree species richness and abundance:
###

# For tree species richness:

# Create table with variables
ct<-data.frame(value=c(bcn$SR_total, bcn$SR_parcs, bcn$SR_viari, bcn$SR_zona),
                  cat=c(rep("Total", nrow(bcn)), rep("Parks", nrow(bcn)), 
                        rep("Streets", nrow(bcn)), rep("Zonal", nrow(bcn))))

# Assign order of factors:
ct$cat<- factor(ct$cat, levels = c("Parks", "Streets", "Zonal", "Total"))

# Perform Kruskal-Wallis test with posthoc comparisons:
kruskal.test(value ~ cat, data = ct)
dunnTest(value~cat,data=ct)

# Plot:
dp1 <- ggplot(ct, aes(x=cat, y=value, fill=cat)) + 
  geom_violin(trim=FALSE)+
  geom_boxplot(width=0.1, fill="white")+
  labs(x="", y = "Species richness")+
  geom_vline(xintercept = 3.5, linetype = "longdash") +
  theme_classic()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        strip.text = element_text(size = 15),
        legend.position="none",
        legend.title = element_blank())

# For tree abundance:

# Create table with variables
ct<-data.frame(value=c(bcn$NT_total, bcn$NT_parcs, bcn$NT_viari, bcn$NT_zona),
               cat=c(rep("Total", nrow(bcn)), rep("Parks", nrow(bcn)), 
                     rep("Streets", nrow(bcn)), rep("Zonal", nrow(bcn))))

# Assign order of factors:
ct$cat<- factor(ct$cat, levels = c("Parks", "Streets", "Zonal", "Total"))

# Perform Kruskal-Wallis test with posthoc comparisons:
kruskal.test(value ~ cat, data = ct)
dunnTest(value~cat,data=ct)

# Plot:
dp2 <- ggplot(ct, aes(x=cat, y=value, fill=cat)) + 
  geom_violin(trim=FALSE)+
  geom_boxplot(width=0.1, fill="white")+
  labs(x="", y = "Tree abundance")+
  geom_vline(xintercept = 3.5, linetype = "longdash") +
  theme_classic()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        strip.text = element_text(size = 15),
        legend.position="none",
        legend.title = element_blank())

# Put plots together:
p<-ggarrange(dp1, dp2,
             labels=c("a)", "b)"),
             ncol = 2, nrow = 1)

# Save result and add letters manually:
png("results/Fig2.png",
    res=600, height=3,width=7,units="in"); 
p
grid.text("a", x = unit(0.135, "npc"), y = unit(0.963, "npc"), 
          gp=gpar(fontsize=10, fontface="bold"))
grid.text("b", x = unit(0.236, "npc"), y = unit(0.963, "npc"), 
          gp=gpar(fontsize=10, fontface="bold"))
grid.text("b", x = unit(0.333, "npc"), y = unit(0.963, "npc"), 
          gp=gpar(fontsize=10, fontface="bold"))

grid.text("a", x = unit(0.654, "npc"), y = unit(0.963, "npc"), 
          gp=gpar(fontsize=10, fontface="bold"))
grid.text("b", x = unit(0.748, "npc"), y = unit(0.963, "npc"), 
          gp=gpar(fontsize=10, fontface="bold"))
grid.text("a", x = unit(0.839, "npc"), y = unit(0.963, "npc"), 
          gp=gpar(fontsize=10, fontface="bold"))
dev.off()


# Clean up environment:
rm(list = ls())

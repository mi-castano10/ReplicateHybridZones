```{r}
## Load the package
library(hzar)
library(tidyverse)
library(raster)
library(rgeos)
library(drc)
library(Hmisc)
library(dr4pl)
library(ggrepel)

### CHANGE THIS ACCORDING TO THE TRANSECT/SEX YOU WANT TO PLOT
rampho<-rampho_MT2
### List here the traits that you are interested in plotting
traits<-c("Q", "FixedSNP1", "FixedSNP2")
traits_morpho<-c("rump.color.hue")

## Compare the 15 cline models from molecular data to the null model graphically for males
for (trait in traits){
  hzar.plot.cline(rampho[[trait]]$analysis$oDG);
}
for (trait in traits_morpho){
  hzar.plot.cline(rampho[[trait]]$analysis$oDG);
}

# Create an empty dataframe to store the molecular results
results_df <- data.frame()

for (trait in traits) {
  model_name <- rampho[[trait]]$analysis$model.name
  if (model_name != "nullModel") {
    param_variation <- hzar.getLLCutParam(rampho[[trait]]$analysis$model.selected,
                                          names(rampho[[trait]]$analysis$model.selected$data.param))
    ml_cline <- hzar.get.ML.cline(rampho[[trait]]$analysis$model.selected)
  
    # Extract center and width from MLcline and pMin and pMax
    center <- ml_cline$param.free$center
    width <- ml_cline$param.free$width
    pMin <- ml_cline$param.all$pMin
    pMax <- ml_cline$param.all$pMax

    # Store the results in the dataframe
    results_df <- bind_rows(results_df, data.frame(Trait = trait,
                                                   Model= rampho[[trait]]$analysis$model.name,
                                                   Parameter_Variation = param_variation,
                                                   Center = center,
                                                   Width = width,
                                                   pMin = pMin,
                                                   pMax = pMax))
  } else {
    cat("Skipping trait", trait, "as the selected model is nullModel\n")
  }
}

# Print the results dataframe
results_df2<-dplyr::select(results_df, Trait, Model, Center, Parameter_Variation.center2LLLow, Parameter_Variation.center2LLHigh, Width, Parameter_Variation.width2LLLow, Parameter_Variation.width2LLHigh, pMin, pMax)
colnames(results_df2)<-c("Trait","Model","Center","CenterLow", "CenterHigh","Width","WidthLow","WidthHigh", "pMin", "pMax")
print(results_df2)

# Create an empty dataframe to store the morphological results
results_df_morpho <- data.frame()

## Look at the variation in parameters for the selected model
for (trait in traits_morpho) {
  param_variation_morpho <- hzar.getLLCutParam(rampho[[trait]]$analysis$model.selected,
                            names(rampho[[trait]]$analysis$model.selected$data.param))
  ml_cline_morpho <- hzar.get.ML.cline(rampho[[trait]]$analysis$model.selected)
  
  # Extract center and width from MLcline
  center <- ml_cline_morpho$param.free$center
  width <- ml_cline_morpho$param.free$width
  
  # Store the results in the dataframe
  results_df_morpho <- bind_rows(results_df_morpho, data.frame(Trait = trait,
                                                               Model = rampho[[trait]]$analysis$model.name,
                                                 Parameter_Variation = param_variation_morpho,
                                                 Center = center,
                                                 Width = width))
}

results_df_morpho2<-dplyr::select(results_df_morpho, Trait, Model,Center, Parameter_Variation.center2LLLow, Parameter_Variation.center2LLHigh, Width, Parameter_Variation.width2LLLow, Parameter_Variation.width2LLHigh)

# Print the results dataframe
colnames(results_df_morpho2)<-c("Trait","Model","Center","CenterLow", "CenterHigh","Width","WidthLow","WidthHigh")
print(results_df_morpho2)

## APPEND THE RESULTS OF BOTH TABLES
ResultsTableMT2<-bind_rows(results_df_morpho2,results_df2) 
ResultsTableMT2$Trait <- factor(ResultsTableMT2$Trait, levels = unique(ResultsTableMT2$Trait))
ResultsTable<-ResultsTableMT2

pdf('./T2/Figures/FINAL_MT2_Clines.pdf', width=7, height = 5)
hzar.plot.cline(rampho$Q$analysis$model.selected,col="black",pch=17, lwd=2, ylim=c(0,1.1), xlim=c(-100,100))
for (trait in T1traits) {
  if (!(trait %in% c("Q"))){ #plot all the traits that are not the hybrid index in the same color
  hzar.plot.cline(rampho[[trait]]$analysis$model.selected, add=T, col="darkgrey",pch=17, lwd=2)
  }
}
hzar.plot.cline(rampho$rump.color.hue$analysis$model.selected,add=T,col="darkorange2",pch=17, lwd=2)
hzar.plot.cline(rampho$Q$analysis$model.selected,col="black",pch=17, lwd=2, add=T)

#Add a legend
#legend("topleft", legend=c("Q","Rump Color (Hue)","Fixed SNPs","0.5 Isocline", "Elevation", "Tree line"), col=c("black","darkorange2","darkgrey","red","blue","forestgreen"), lty=c(1,1,1,1,3,5), lwd =c(2,2,2,2,2,1.5), cex = 0.9, bg="white")

abline(v=ResultsTable[ResultsTable$Trait == "Q", "Center"], col="black", lty=2, lwd=1.5)
abline(v=ResultsTable[ResultsTable$Trait == "rump.color.hue", "Center"], col="darkorange2", lty=2, lwd=1.5)
abline(v=0, col="red", lty=1, lwd=3)
abline(v=-11.88, col="blue", lty=3, lwd=3)
abline(v=3.13, col="forestgreen", lty=5, lwd=2.5)

# Add arrows above specific x-values
#triangle_x_values <- ResultsTable$Center[[5:nrow(ResultsTable)]] # Specify x-values where you want to add arrows
triangle_x_values <- ResultsTable$Center
triangle_y <- 1.12  # Specify y-coordinate of the arrows
points(triangle_x_values, rep(triangle_y, length(triangle_x_values)), pch = 25, bg= c("darkorange","black","darkgrey","darkgrey"), col = "black", cex=1.3)

dev.off()
```

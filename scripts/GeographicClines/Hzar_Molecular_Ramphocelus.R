### Modified code from original script of the R package hzar: hybrid zone analysis using an R software package 

library(hzar)
library(doMC)
library(foreach)
library(tidyverse)
#library(raster)
#library(rgeos)
#library(drc)
#library(Hmisc)
#library(dr4pl)
#library(ggrepel)
#library(sf)
#library(geosphere)

## If you have doMC, use for each in parallel mode
## to speed up computation.
if(require(doMC)){
  registerDoMC(cores=8)
} else {
  ## Use foreach in sequential mode
  registerDoSEQ();
}

# Example to run clines for Males of Transect 2
# READ IN THE DATA (this should be a .txt in the format shown in the Hzar tutorial as "manakinMolecular").
# THIS FILE SHOULD CONTAIN THE ALLELE FREQUENCIES FOR ALL THE SNPs YOU ARE GOING TO INCLUDE IN THE ANALYSIS! 
ramphocelusMolecular<-read_delim("/scratch/juy3_lab/GBS_Ramphocelus/ClineAnalysis/T2/ramphocelusMolecular_MT2.txt", 
                                 delim = "\t", escape_double = FALSE, trim_ws = TRUE)

# Molecular Analysis
# We can do as many traits/ specific alleles at one locus as we want. 
# So I did this for loop to create a list to stay organized.

rampho<-list()

# First pick the traits you want to build your clines for
# In this case the Q population value and
# I'm using the Major allele of every SNP
# ua.nSamples <- "Q.nSamples"

# Get all column names in the ramphocelusMolecular dataframe
all_colnames <- colnames(ramphocelusMolecular)

# To do the analysis on all the traits present in the table "ramphocelusMolecular" exclude the first two columns and columns ending with .nSamples

traits <- all_colnames[!grepl(".nSamples$", all_colnames) & !(all_colnames %in% all_colnames[1:4])]

# Alternatively you can use a list of SNPs (in a .txt file) as the "traits" variable or tell it specifically which traits to use then comment out the line above and uncomment either of the following:
#traits <- readLines("/scratch/juy3_lab/GBS_Ramphocelus/ClineAnalysis/T2/FstMT2.txt")
#traits <- c("Q", "HI") 
# The names of the traits should be exactly the same as they are in the original ramphocelusMolecular table

# Blank out space in memory to hold molecular analysis
if(length(apropos("^rampho$",ignore.case=FALSE)) == 0 ||
   !is.list(rampho) ) rampho <- list()

# By executing this code, you will have the desired structure with five empty 
# lists for each trait in the rampho list.

for (trait in traits) {
  rampho[[trait]]<-list(
    obs=list(),
    models=list(),
    fitRs=list(),
    runs=list(),
    analysis=list()
  )
}

# Each trait in the list rampho has:
# Space to hold the observed data
# Space to hold the models to fit
# Space to hold the compiled fit requests
# Space to hold the output data chains
# Space to hold the analysed data

# Iterate over each trait
for (trait in traits) {
  # Remove the last letter from the trait name
  trait_cleaned <- sub("\\..*$", "", trait)
  
  # Perform the operation for the current trait
  rampho[[trait]][["obs"]] <- hzar.doMolecularData1DPops(ramphocelusMolecular$distance,
                                                         ramphocelusMolecular[[trait]],
                                                         ramphocelusMolecular[[paste0(trait_cleaned, ".nSamples")]])}

# Next we want to fit different cline models for molecular markers that influence the shape of the cline tails.
# if the scaling is fixed the mean is fixed at 0 or 1 and variance fixed at 0 for either end of the cline. 

# we make a list of all the models we want to test for each trait
# so that each model has a different combination of parameters as follows:
# Model 1 (G) has fixed scaling and no tails.
# Model 2 (G) has fixed scaling and mirrored tails. 
# Model 3 (G) has fixed scaling and both tails. 
# Model 4 (G+M) has free scaling and no tails. 
# Model 5 (G+M) has free scaling and mirrored tails. 
# Model 6 (G+M) has free scaling and both tails.
# Mu represents the mean allele frequencies and phenotypic trait values for the left and right cline tails. 
# Delta and tau are exponential decay curve parameters for left and right tails.

paramListPerTrait<-list(
  list(scaling = "none", tails = "none", id = "model1"),
  list(scaling = "none", tails = "right", id = "model2"),
  list(scaling = "none", tails = "left", id = "model3"),
  list(scaling = "none", tails = "mirror", id = "model4"),
  list(scaling = "none", tails = "both", id = "model5"),
  list(scaling = "fixed", tails = "none", id = "model6"),
  list(scaling = "fixed", tails = "right", id = "model7"),
  list(scaling = "fixed", tails = "left", id = "model8"),
  list(scaling = "fixed", tails = "mirror", id = "model9"),
  list(scaling = "fixed", tails = "both", id = "model10"),
  list(scaling = "free", tails = "none", id = "model11"),
  list(scaling = "free", tails = "right", id = "model12"),
  list(scaling = "free", tails = "left", id = "model13"),
  list(scaling = "free", tails = "mirror", id = "model14"),
  list(scaling = "free", tails = "both", id = "model15")
)

# A list to keep track of traits to be removed
traits_to_remove <- c()

for (trait in traits) {
  # Check for and handle missing values in the dist and obsFreq columns
  if (any(is.na(rampho[[trait]]$obs$dist)) || any(is.na(rampho[[trait]]$obs$obsFreq))) {
    cat("Warning: Missing values found in trait", trait, "- cleaning data.\n")
    rampho[[trait]]$obs <- na.omit(rampho[[trait]]$obs)
  }
  
  all_models_successful <- TRUE
  
  for (params in paramListPerTrait) {
    model_id <- paste0(params$id) # Creates a unique model ID for each model within each trait

    # Create a metadata object using the hzar.makeCline1DFreq function
    metadata <- tryCatch({
      hzar.makeCline1DFreq(rampho[[trait]]$obs, params$scaling, params$tails)
    }, error = function(e) {
      cat("Error in hzar.makeCline1DFreq for trait", trait, model_id, ":", e$message, "\n")
      NULL
    })
    
    # Store the metadata object in the respective trait if it was created successfully
    if (!is.null(metadata)) {
      rampho[[trait]]$models[[model_id]] <- metadata
    } else {
      cat("Failed to create metadata for model ID:", model_id, "\n")
      all_models_successful <- FALSE
      break
    }
  }
  
  # If any model creation failed, mark the trait for removal
  if (!all_models_successful) {
    traits_to_remove <- c(traits_to_remove, trait)
  }
}

# Remove the traits marked for removal from the rampho list
for (trait in traits_to_remove) {
  cat("Removing trait:", trait, "due to unsuccessful model creation.\n")
  rampho[[trait]] <- NULL
}

# Filter out NULL elements from rampho
rampho <- rampho[!sapply(rampho, is.null)]


# Modify all models to focus on the region where the data was collected
# Making them all go from -100 to 100 so that we have homogeneous plots

for (trait in traits) {
  rampho[[trait]]$models <- sapply(rampho[[trait]]$models,
                                   hzar.model.addBoxReq,
                                   -100 , 100,
                                   simplify=FALSE)
  
}

# Compile each of the models of each trait to prepare for fitting

for (trait in traits) {
  rampho[[trait]]$fitRs$init <- sapply(rampho[[trait]]$models,
                                       hzar.first.fitRequest.old.ML,
                                       obsData=rampho[[trait]]$obs,
                                       verbose=FALSE,
                                       simplify=FALSE)
}

# You can update the settings for the fitter of each model if desired.

## A typical chain length.  This value is the default setting in the package.
chainLength=1e6;                

## Make each model run off a separate seed
mainSeed<-list(
  A=c(596,528,124,978,544,99),
  B=c(528,124,978,544,99,596),
  C=c(124,978,544,99,596,528),
  D=c(546,578,194,678,533,89),
  E=c(569,303,818,538,334,607),
  G=c(786,48,570,774,302,450),
  H=c(291,635,447,381,33,234),
  I=c(826,870,697,18,918,260),
  J=c(761,982,412,905,501,206),
  K=c(776,204,628,907,891,950),
  L=c(49,176,882,187,81,273),
  M=c(639,594,263,249,408,436),
  N=c(988,749,652,707,204,505),
  O=c(897,885,924,674,991,978),
  P=c(17,680,428,555,239,765),
  Q=c(799,245,674,785,625,442))

# update the SEED, the chain_length and the burning and make each model run on a different Seed. 
# this will modify 5 seed values for each letter that will be assigned to each 
# model in each trait

for (trait in traits) {
  for (i in 1:length(paramListPerTrait)) {
    model_id <- paramListPerTrait[[i]]$id
    seed <- mainSeed[[i]]
    # Update seed each model with the list of seed values for each one (in the mainSeed list)
    rampho[[trait]]$fitRs$init[[model_id]]$mcmcParam$seed[[1]] <- seed
    # Update chain length
    rampho[[trait]]$fitRs$init[[model_id]]$mcmcParam$chainLength <- chainLength
    # Update burnin
    rampho[[trait]]$fitRs$init[[model_id]]$mcmcParam$burnin <- chainLength %/% 10
    # Update burnin
    rampho[[trait]]$fitRs$init[[model_id]]$mcmcParam$thin <- 100
  }
}

## Create an empty "init" sublist inside the runs list for every trait
for (trait in traits) {
  rampho[[trait]]$runs$init <- list()
}

#Iterate through the model list and populate it with results of the initial chain
for (j in seq_along(traits)) {
  trait <- traits[j]
  for (i in seq_along(rampho[[trait]]$models)) {
    model <- rampho[[trait]]$models[i]
    rampho[[j]]$runs$init[[i]] <- hzar.doFit(rampho[[j]]$fitRs$init[[i]])
  }
}

#Change the names of the init items back to the model names
for (trait in traits) {
  names(rampho[[trait]]$runs$init) <- c("model1",
                                        "model2",
                                        "model3",
                                        "model4",
                                        "model5",
                                        "model6",
                                        "model7",
                                        "model8",
                                        "model9",
                                        "model10",
                                        "model11",
                                        "model12",
                                        "model13",
                                        "model14",
                                        "model15")
}

## Compile a new set of fit requests using the initial chains 
for (trait in traits) {
  rampho[[trait]]$fitRs$chains <-
    lapply(rampho[[trait]]$runs$init, 
           hzar.next.fitRequest)
}

## Replicate each fit request 3 times, keeping the original
## seeds while switching to a new seed channel.
for (trait in traits) {
  rampho[[trait]]$fitRs$chains <-
    hzar.multiFitRequest(rampho[[trait]]$fitRs$chains,
                         each=3,
                         baseSeed=NULL)
}

## Just to be thorough, randomize the initial value for each fit

#CENTER FOR ALL MODELS
# Set the range for the randomization
min_value <- -100
max_value <- 100

# Iterate through each trait in the rampho list
for (trait in traits) {
  # Iterate through each model inside the chains list within the trait
  for (i in seq_along(rampho[[trait]]$fitRs$chains)) {
    # Generate a random value between min_value and max_value
    random_value <- runif(1, min_value, max_value)
    # Assign the random value to the "center" parameter of the current model
    rampho[[trait]]$fitRs$chains[[i]]$modelParam$init["center"] <- random_value
  }
}

#WIDTH FOR ALL MODELS
# Set the range for the randomization
min_value <- 0
max_value <- 50

# Iterate through each trait in the rampho list
for (trait in traits) {
  # Iterate through each model inside the chains list within the trait
  for (i in seq_along(rampho[[trait]]$fitRs$chains)) {
    # Generate a random value between min_value and max_value
    random_value <- runif(1, min_value, max_value)
    # Assign the random value to the "center" parameter of the current model
    rampho[[trait]]$fitRs$chains[[i]]$modelParam$init["width"] <- random_value
  }
}

#pMin (Mu) FOR ALL MODELS THAT INCLUDE IT 
# Set the range for the randomization
min_value <- 0
max_value <- 1

# Iterate through each trait in the rampho list
for (trait in traits) {
  # Iterate through each model inside the chains list within the trait
  for (i in seq_along(rampho[[trait]]$fitRs$chains)) {
    # Check if the "pMin" parameter exists in the current model
    if ("pMin" %in% names(rampho[[trait]]$fitRs$chains[[i]]$modelParam$init)) {
      # Generate a random value between min_value and max_value
      random_value <- runif(1, min_value, max_value)
      # Assign the random value to the "center" parameter of the current model
      rampho[[trait]]$fitRs$chains[[i]]$modelParam$init["pMin"] <- random_value
    }
  }
}

#pMax (Mu) FOR ALL MODELS THAT INCLUDE IT 
# Set the range for the randomization
min_value <- 0
max_value <- 1

# Iterate through each trait in the rampho list
for (trait in traits) {
  # Iterate through each model inside the chains list within the trait
  for (i in seq_along(rampho[[trait]]$fitRs$chains)) {
    # Check if the "pMax" parameter exists in the current model
    if ("pMax" %in% names(rampho[[trait]]$fitRs$chains[[i]]$modelParam$init)) {
      # Generate a random value between min_value and max_value
      random_value <- runif(1, min_value, max_value)
      # Assign the random value to the "center" parameter of the current model
      rampho[[trait]]$fitRs$chains[[i]]$modelParam$init["pMax"] <- random_value
    }
  }
}

#deltaL FOR ALL MODELS THAT INCLUDE IT 
# Set the range for the randomization
min_value <- 0
max_value <- 1

# Iterate through each trait in the rampho list
for (trait in traits) {
  # Iterate through each model inside the chains list within the trait
  for (i in seq_along(rampho[[trait]]$fitRs$chains)) {
    # Check if the "deltaL" parameter exists in the current model
    if ("deltaL" %in% names(rampho[[trait]]$fitRs$chains[[i]]$modelParam$init)) {
      # Generate a random value between min_value and max_value
      random_value <- runif(1, min_value, max_value)
      # Assign the random value to the "center" parameter of the current model
      rampho[[trait]]$fitRs$chains[[i]]$modelParam$init["deltaL"] <- random_value
    }
  }
}

#tauL FOR ALL MODELS THAT INCLUDE IT 
# Set the range for the randomization
min_value <- 0
max_value <- 1

# Iterate through each trait in the rampho list
for (trait in traits) {
  # Iterate through each model inside the chains list within the trait
  for (i in seq_along(rampho[[trait]]$fitRs$chains)) {
    # Check if the "tauL" parameter exists in the current model
    if ("tauL" %in% names(rampho[[trait]]$fitRs$chains[[i]]$modelParam$init)) {
      # Generate a random value between min_value and max_value
      random_value <- runif(1, min_value, max_value)
      # Assign the random value to the "center" parameter of the current model
      rampho[[trait]]$fitRs$chains[[i]]$modelParam$init["tauL"] <- random_value
    }
  }
}

#deltaR FOR ALL MODELS THAT INCLUDE IT 
# Set the range for the randomization
min_value <- 0
max_value <- 1

# Iterate through each trait in the rampho list
for (trait in traits) {
  # Iterate through each model inside the chains list within the trait
  for (i in seq_along(rampho[[trait]]$fitRs$chains)) {
    # Check if the "deltaR" parameter exists in the current model
    if ("deltaR" %in% names(rampho[[trait]]$fitRs$chains[[i]]$modelParam$init)) {
      # Generate a random value between min_value and max_value
      random_value <- runif(1, min_value, max_value)
      # Assign the random value to the "center" parameter of the current model
      rampho[[trait]]$fitRs$chains[[i]]$modelParam$init["deltaR"] <- random_value
    }
  }
}

#tauR FOR ALL MODELS THAT INCLUDE IT 
# Set the range for the randomization
min_value <- 0
max_value <- 1

# Iterate through each trait in the rampho list
for (trait in traits) {
  # Iterate through each model inside the chains list within the trait
  for (i in seq_along(rampho[[trait]]$fitRs$chains)) {
    # Check if the "tauR" parameter exists in the current model
    if ("tauR" %in% names(rampho[[trait]]$fitRs$chains[[i]]$modelParam$init)) {
      # Generate a random value between min_value and max_value
      random_value <- runif(1, min_value, max_value)
      # Assign the random value to the "center" parameter of the current model
      rampho[[trait]]$fitRs$chains[[i]]$modelParam$init["tauR"] <- random_value
    }
  }
}

## Go ahead and run a chain of 3 runs for every fit request
for (trait in traits) {
  rampho[[trait]]$runs$chains <-  
    hzar.doChain.multi(rampho[[trait]]$fitRs$chains,
                       doPar=TRUE,
                       inOrder=FALSE,
                       count=3)                     
}

# Start aggregation of data for analysis

# Create a model data group for the null model (expected allele
# frequency independent of distance along cline) to include in
# analysis.

for (trait in traits) {
  rampho[[trait]]$analysis$initDGs <- list(
    nullModel =  hzar.dataGroup.null(rampho[[trait]]$obs))
}

# Create a model data group (hzar.dataGroup object) for each model from the initial runs.

for (trait in traits){
  null_model_dg<-rampho[[trait]]$analysis$initDGs[[1]]
  model_dgs<-list(null_model_dg)
  for (j in seq_along(rampho[[trait]]$models)){
    model<-rampho[[trait]]$models[j]
    model_dg <- hzar.dataGroup.add(rampho[[trait]]$runs$init[[j]])
    model_dgs<-c(model_dgs, list(model_dg))
  }
  rampho[[trait]]$analysis$initDGs <- model_dgs 
}

#name the models again
for (trait in traits) {
  names(rampho[[trait]]$analysis$initDGs) <- c("nullModel",
                                        "model1",
                                        "model2",
                                        "model3",
                                        "model4",
                                        "model5",
                                        "model6",
                                        "model7",
                                        "model8",
                                        "model9",
                                        "model10",
                                        "model11",
                                        "model12",
                                        "model13",
                                        "model14",
                                        "model15")
}

# Create a hzar.obsDataGroup object from the four hzar.dataGroup just created, 
# copying the naming scheme (nullModel, modelI, modelII, modelIII).

for (trait in traits){
  rampho[[trait]]$analysis$oDG <-
    hzar.make.obsDataGroup(rampho[[trait]]$analysis$initDGs)
}

for (trait in traits){
  rampho[[trait]]$analysis$oDG <-
    hzar.copyModelLabels(rampho[[trait]]$analysis$initDGs,
                         rampho[[trait]]$analysis$oDG)
  
}


# Convert all runs to hzar.dataGroup objects, adding them to
# the hzar.obsDataGroup object.
for (trait in traits){
  rampho[[trait]]$analysis$oDG <-
    hzar.make.obsDataGroup(lapply(rampho[[trait]]$runs$chains,
                                  hzar.dataGroup.add),
                           rampho[[trait]]$analysis$oDG);
}

# Do model selection based on the AICc scores
for (trait in traits){
  print(rampho[[trait]]$analysis$AICcTable <- 
          hzar.AICc.hzar.obsDataGroup(rampho[[trait]]$analysis$oDG));
}

# Print out the model with the minimum AICc score
for (trait in traits){
  print(rampho[[trait]]$analysis$model.name <-
          rownames(rampho[[trait]]$analysis$AICcTable
          )[[ which.min(rampho[[trait]]$analysis$AICcTable$AICc )]])
}

# Extract the hzar.dataGroup object for the selected model
for (trait in traits){
  rampho[[trait]]$analysis$model.selected <-
    rampho[[trait]]$analysis$oDG$data.groups[[rampho[[trait]]$analysis$model.name]]
}

#save your (potentially very heavy) Rdata object to plot later.
rampho_MT2<-rampho
save(rampho_MT2, file="/scratch/juy3_lab/GBS_Ramphocelus/ClineAnalysis/T2/rampho_MT2.RData")

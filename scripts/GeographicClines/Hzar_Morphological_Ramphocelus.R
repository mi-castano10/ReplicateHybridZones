### Modified code from original script of the R package hzar: hybrid zone analysis using an R software package 

library(hzar)
library(doMC)
library(foreach)
library(tidyverse)


# If you have doMC, use for each in parallel mode
# to speed up computation.
if(require(doMC)){
  registerDoMC(cores=8)
} else {
  ## Use foreach in sequential
  mode
  registerDoSEQ();
}

# READ IN THE DATA
# Note that this is the individual data, with labels identifying the locality of each sample. This should be the same format as "manakinMorphological" in the original Hzar tutorial.
Morphological<-read.csv("/scratch/juy3_lab/GBS_Ramphocelus/ClineAnalysis/T2/ramphocelusMorphological_MT2.txt", h=T, row.names = 1)
Morphological$locality<-as.factor(Morphological$locality)

# Load locality data matching each localility to a site ID and a transect distance. This should be the same format as "manakinLocality" in the original Hzar tutorial.
Locations<-read.table("/scratch/juy3_lab/GBS_Ramphocelus/ClineAnalysis/T2/T2_Locality.txt", h=T)
Locations$locality<-as.factor(Locations$locality) 

# Morphological Analysis
# We can do as many traits as we want. 
# So I did this for loop to create a list to stay organized.

rampho_morpho<-list()

# First pick the traits you want to build your clines for
# In this case color and body size. 
# This should have the exact same names as the traits in manakinMorphological
traits_morpho<-c("rump.color.hue","rump.color.chroma","body.size")

# Blank out space in memory to hold molecular analysis
if(length(apropos("^rampho_morpho$",ignore.case=FALSE)) == 0 ||
   !is.list(rampho_morpho) ) rampho_morpho <- list()

# By executing this code, you will have the desired structure with five empty 
# lists for each trait in the rampho_morpho list.

for (trait in traits_morpho) {
  rampho_morpho[[trait]]<-list(
    obs=list(),
    models=list(),
    fitRs=list(),
    runs=list(),
    analysis=list()
  )
}
## Each trait in the list rampho has:
## Space to hold the observed data
## Space to hold the models to fit
## Space to hold the compiled fit requests
## Space to hold the output data chains
## Space to hold the analysed data

# Iterate over each transect (change transect name accordingly)

for (trait in traits_morpho){
  rampho_morpho[[trait]]$obs <-
    hzar.doNormalData1DRaw(hzar.mapSiteDist(Locations$locality,
                                            Locations$distance),
                           Morphological$locality,
                           Morphological[[trait]])
}

# We want to fit 3 different cline models for morphological (quantitative traits).
# All models estimate trait mean and variance on the left and right and additional variance in the centre as well as centre and width. 
# we make a list of all the models we want to test for each trait
# so that each model has a different combination of parameters as follows:

# Model 4 (G+M) has free scaling and no tails. 
# Model 5 (G+M) has free scaling and mirrored tails. 
# Model 6 (G+M) has free scaling and both tails.
paramListPerTrait_morpho<-list(
  list(scaling = "free", tails = "none", id = "model1"),
  list(scaling = "free", tails = "right", id = "model2"),
  list(scaling = "free", tails = "left", id = "model3"),
  list(scaling = "free", tails = "mirror", id = "model4"),
  list(scaling = "free", tails = "both", id = "model5")
)

for (trait in traits_morpho) {
  for (params in paramListPerTrait_morpho) {
    model_id <- paste0(params$id) #creates a unique model ID for each model within each trait
    
    metadata <- hzar.makeCline1DNormal(rampho_morpho[[trait]]$obs, params$tails) 
    
    #store the metadata object in the respective trait 
    rampho_morpho[[trait]]$models[[model_id]] <- metadata
  }
}

# Modify all models to focus on the region where the observed
# data were collected.

for (trait in traits_morpho) {
  rampho_morpho[[trait]]$models <- sapply(rampho_morpho[[trait]]$models,
                                          hzar.model.addBoxReq,
                                          -100 , 100,
                                          simplify=FALSE)
  
}

# Due to the large number of free variables, it is prudent to
# reduce the tune setting of all models from 1.5 to 1.2

for (trait in traits_morpho){
  for (j in seq_along(rampho_morpho[[trait]]$models)){
    hzar.meta.tune(rampho_morpho[[trait]]$models[[j]])<-1.2
  }
}

# Compile each of the models to prepare for fitting
# Note that we are using hzar.first.fitRequest.gC for fitting
# guassian (aka "normal") clines.

for (trait in traits_morpho) {
  rampho_morpho[[trait]]$fitRs$init <- sapply(rampho_morpho[[trait]]$models,
                                              hzar.first.fitRequest.gC,
                                              obsData=rampho_morpho[[trait]]$obs,
                                              verbose=FALSE,
                                              simplify=FALSE)
}

# You can update the settings for the fitter if desired.
# If you do do it for each model! 

## A typical chain length.  This value is the default setting in the package.
chainLength=1e6;                       

## Make each model run off a separate seed
mainSeed<-list(
  A=c(596,528,124,978,544,99),
  B=c(528,124,978,544,99,596),
  C=c(124,978,544,99,596,528),
  D=c(546,578,194,678,533,89),
  E=c(569,303,818,538,334,607))

# update the SEED and make each model run on a different Seed. 
# this will modify 5 seed values for each letter that will be assigned to each model in each trait

for (trait in traits_morpho) {
  for (i in 1:length(paramListPerTrait_morpho)) {
    model_id<- paramListPerTrait_morpho[[i]]$id
    seed<-mainSeed[[i]] 
    # Update seed each model with the list of seed values for each one (in the mainSeed list)
    rampho_morpho[[trait]]$fitRs$init[[model_id]]$mcmcParam$seed[[1]]<-seed
    # Update chain length
    rampho_morpho[[trait]]$fitRs$init[[model_id]]$mcmcParam$chainLength <- chainLength
    # Update burnin
    rampho_morpho[[trait]]$fitRs$init[[model_id]]$mcmcParam$burnin <- chainLength %/% 10
    # Update burnin
    rampho_morpho[[trait]]$fitRs$init[[model_id]]$mcmcParam$thin <- 100
  }
}

# Create an empty "init" sublist inside the runs list for every trait
for (trait in traits_morpho) {
  rampho_morpho[[trait]]$runs$init <- list()
}

# I terate through the model list and populate it with results from the run of the initial chain
for (j in seq_along(traits_morpho)) {
  trait <- traits_morpho[j]
  for (i in seq_along(rampho_morpho[[trait]]$models)) {
    model <- rampho_morpho[[trait]]$models[i]
    rampho_morpho[[j]]$runs$init[[i]] <- hzar.doFit(rampho_morpho[[j]]$fitRs$init[[i]])
  }
}

# Change the names of the init items back to the model names
for (trait in traits_morpho) {
  names(rampho_morpho[[trait]]$runs$init) <- c("model1",
                                               "model2",
                                               "model3",
                                               "model4",
                                               "model5")
}

# Compile a new set of fit requests using the initial chains 
for (trait in traits_morpho) {
  rampho_morpho[[trait]]$fitRs$chains <-
    lapply(rampho_morpho[[trait]]$runs$init, 
           hzar.next.fitRequest)
}

# Replicate each fit request 3 times, keeping the original
# seeds while switching to a new seed channel.
for (trait in traits_morpho) {
  rampho_morpho[[trait]]$fitRs$chains <-
    hzar.multiFitRequest(rampho_morpho[[trait]]$fitRs$chains,
                         each=3,
                         baseSeed=NULL)
}

# Just to be thorough, randomize the initial value for each fit

#CENTER FOR ALL MODELS
# Set the range for the randomization
min_value <- -100
max_value <- 100

# Iterate through each trait in the rampho list
for (trait in traits_morpho) {
  # Iterate through each model inside the chains list within the trait
  for (i in seq_along(rampho_morpho[[trait]]$fitRs$chains)) {
    # Generate a random value between min_value and max_value
    random_value <- runif(1, min_value, max_value)
    # Assign the random value to the "center" parameter of the current model
    rampho_morpho[[trait]]$fitRs$chains[[i]]$modelParam$init["center"] <- random_value
  }
}

#WIDTH FOR ALL MODELS
# Set the range for the randomization
min_value <- 0
max_value <- 50

# Iterate through each trait in the rampho list
for (trait in traits_morpho) {
  # Iterate through each model inside the chains list within the trait
  for (i in seq_along(rampho_morpho[[trait]]$fitRs$chains)) {
    # Generate a random value between min_value and max_value
    random_value <- runif(1, min_value, max_value)
    # Assign the random value to the "center" parameter of the current model
    rampho_morpho[[trait]]$fitRs$chains[[i]]$modelParam$init["width"] <- random_value
  }
}

#varH FOR ALL MODELS THAT INCLUDE IT 
# Set the range for the randomization
min_value <- 0
max_value <- 5

# Iterate through each trait in the rampho list
for (trait in traits_morpho) {
  # Iterate through each model inside the chains list within the trait
  for (i in seq_along(rampho_morpho[[trait]]$fitRs$chains)) {
    # Check if the "varH" parameter exists in the current model
    if ("varH" %in% names(rampho_morpho[[trait]]$fitRs$chains[[i]]$modelParam$init)) {
      # Generate a random value between min_value and max_value
      random_value <- runif(1, min_value, max_value)
      # Assign the random value to the "center" parameter of the current model
      rampho_morpho[[trait]]$fitRs$chains[[i]]$modelParam$init["varH"] <- random_value
    }
  }
}

#muL FOR ALL MODELS THAT INCLUDE IT 
# Set the range for the randomization
min_value <- 0
max_value <- 1

# Iterate through each trait in the rampho list
for (trait in traits_morpho) {
  # Iterate through each model inside the chains list within the trait
  for (i in seq_along(rampho_morpho[[trait]]$fitRs$chains)) {
    # Check if the "muL" parameter exists in the current model
    if ("muL" %in% names(rampho_morpho[[trait]]$fitRs$chains[[i]]$modelParam$init)) {
      # Generate a random value between min_value and max_value
      random_value <- runif(1, min_value, max_value)
      # Assign the random value to the "center" parameter of the current model
      rampho_morpho[[trait]]$fitRs$chains[[i]]$modelParam$init["muL"] <- random_value
    }
  }
}

#muR FOR ALL MODELS THAT INCLUDE IT 
# Set the range for the randomization
min_value <- 0
max_value <- 1

# Iterate through each trait in the rampho list
for (trait in traits_morpho) {
  # Iterate through each model inside the chains list within the trait
  for (i in seq_along(rampho_morpho[[trait]]$fitRs$chains)) {
    # Check if the "muR" parameter exists in the current model
    if ("muR" %in% names(rampho_morpho[[trait]]$fitRs$chains[[i]]$modelParam$init)) {
      # Generate a random value between min_value and max_value
      random_value <- runif(1, min_value, max_value)
      # Assign the random value to the "center" parameter of the current model
      rampho_morpho[[trait]]$fitRs$chains[[i]]$modelParam$init["muR"] <- random_value
    }
  }
}

#varL FOR ALL MODELS THAT INCLUDE IT 
# Set the range for the randomization
min_value <- 0
max_value <- 5

# Iterate through each trait in the rampho list
for (trait in traits_morpho) {
  # Iterate through each model inside the chains list within the trait
  for (i in seq_along(rampho_morpho[[trait]]$fitRs$chains)) {
    # Check if the "tauL" parameter exists in the current model
    if ("varL" %in% names(rampho_morpho[[trait]]$fitRs$chains[[i]]$modelParam$init)) {
      # Generate a random value between min_value and max_value
      random_value <- runif(1, min_value, max_value)
      # Assign the random value to the "center" parameter of the current model
      rampho_morpho[[trait]]$fitRs$chains[[i]]$modelParam$init["varL"] <- random_value
    }
  }
}

#varR FOR ALL MODELS THAT INCLUDE IT 
# Set the range for the randomization
min_value <- 0
max_value <- 5

# Iterate through each trait in the rampho list
for (trait in traits_morpho) {
  # Iterate through each model inside the chains list within the trait
  for (i in seq_along(rampho_morpho[[trait]]$fitRs$chains)) {
    # Check if the "varR" parameter exists in the current model
    if ("varR" %in% names(rampho_morpho[[trait]]$fitRs$chains[[i]]$modelParam$init)) {
      # Generate a random value between min_value and max_value
      random_value <- runif(1, min_value, max_value)
      # Assign the random value to the "center" parameter of the current model
      rampho_morpho[[trait]]$fitRs$chains[[i]]$modelParam$init["varR"] <- random_value
    }
  }
}

#deltaL FOR ALL MODELS THAT INCLUDE IT 
# Set the range for the randomization
min_value <- 0
max_value <- 1

# Iterate through each trait in the rampho_morpho list
for (trait in traits_morpho) {
  # Iterate through each model inside the chains list within the trait
  for (i in seq_along(rampho_morpho[[trait]]$fitRs$chains)) {
    # Check if the "deltaL" parameter exists in the current model
    if ("deltaL" %in% names(rampho_morpho[[trait]]$fitRs$chains[[i]]$modelParam$init)) {
      # Generate a random value between min_value and max_value
      random_value <- runif(1, min_value, max_value)
      # Assign the random value to the "center" parameter of the current model
      rampho_morpho[[trait]]$fitRs$chains[[i]]$modelParam$init["deltaL"] <- random_value
    }
  }
}

#tauL FOR ALL MODELS THAT INCLUDE IT 
# Set the range for the randomization
min_value <- 0
max_value <- 1

# Iterate through each trait in the rampho_morpho list
for (trait in traits_morpho) {
  # Iterate through each model inside the chains list within the trait
  for (i in seq_along(rampho_morpho[[trait]]$fitRs$chains)) {
    # Check if the "tauL" parameter exists in the current model
    if ("tauL" %in% names(rampho_morpho[[trait]]$fitRs$chains[[i]]$modelParam$init)) {
      # Generate a random value between min_value and max_value
      random_value <- runif(1, min_value, max_value)
      # Assign the random value to the "center" parameter of the current model
      rampho_morpho[[trait]]$fitRs$chains[[i]]$modelParam$init["tauL"] <- random_value
    }
  }
}

#deltaR FOR ALL MODELS THAT INCLUDE IT 
# Set the range for the randomization
min_value <- 0
max_value <- 1

# Iterate through each trait in the rampho_morpho list
for (trait in traits_morpho) {
  # Iterate through each model inside the chains list within the trait
  for (i in seq_along(rampho_morpho[[trait]]$fitRs$chains)) {
    # Check if the "deltaR" parameter exists in the current model
    if ("deltaR" %in% names(rampho_morpho[[trait]]$fitRs$chains[[i]]$modelParam$init)) {
      # Generate a random value between min_value and max_value
      random_value <- runif(1, min_value, max_value)
      # Assign the random value to the "center" parameter of the current model
      rampho_morpho[[trait]]$fitRs$chains[[i]]$modelParam$init["deltaR"] <- random_value
    }
  }
}

#tauR FOR ALL MODELS THAT INCLUDE IT 
# Set the range for the randomization
min_value <- 0
max_value <- 1

# Iterate through each trait in the rampho_morpho list
for (trait in traits_morpho) {
  # Iterate through each model inside the chains list within the trait
  for (i in seq_along(rampho_morpho[[trait]]$fitRs$chains)) {
    # Check if the "tauR" parameter exists in the current model
    if ("tauR" %in% names(rampho_morpho[[trait]]$fitRs$chains[[i]]$modelParam$init)) {
      # Generate a random value between min_value and max_value
      random_value <- runif(1, min_value, max_value)
      # Assign the random value to the "center" parameter of the current model
      rampho_morpho[[trait]]$fitRs$chains[[i]]$modelParam$init["tauR"] <- random_value
    }
  }
}

#deltaM FOR ALL MODELS THAT INCLUDE IT 
# Set the range for the randomization
min_value <- 0
max_value <- 1

# Iterate through each trait in the rampho_morpho list
for (trait in traits_morpho) {
  # Iterate through each model inside the chains list within the trait
  for (i in seq_along(rampho_morpho[[trait]]$fitRs$chains)) {
    # Check if the "deltaM" parameter exists in the current model
    if ("deltaM" %in% names(rampho_morpho[[trait]]$fitRs$chains[[i]]$modelParam$init)) {
      # Generate a random value between min_value and max_value
      random_value <- runif(1, min_value, max_value)
      # Assign the random value to the "center" parameter of the current model
      rampho_morpho[[trait]]$fitRs$chains[[i]]$modelParam$init["deltaM"] <- random_value
    }
  }
}

#tauM FOR ALL MODELS THAT INCLUDE IT 
# Set the range for the randomization
min_value <- 0
max_value <- 1

# Iterate through each trait in the rampho_morpho list
for (trait in traits_morpho) {
  # Iterate through each model inside the chains list within the trait
  for (i in seq_along(rampho_morpho[[trait]]$fitRs$chains)) {
    # Check if the "tauM" parameter exists in the current model
    if ("tauM" %in% names(rampho_morpho[[trait]]$fitRs$chains[[i]]$modelParam$init)) {
      # Generate a random value between min_value and max_value
      random_value <- runif(1, min_value, max_value)
      # Assign the random value to the "center" parameter of the current model
      rampho_morpho[[trait]]$fitRs$chains[[i]]$modelParam$init["tauM"] <- random_value
    }
  }
}

# Go ahead and run a chain of 3 runs for every fit request
for (trait in traits_morpho) {
  rampho_morpho[[trait]]$runs$chains <-  
    hzar.doChain.multi(rampho_morpho[[trait]]$fitRs$chains,
                       doPar=TRUE,
                       inOrder=FALSE,
                       count=3)                     
}

# Start aggregation of data for analysis

# Clear out a spot to collect the data for analysis (note that
# there is currently no "null model" to compare against).

for (trait in traits_morpho) {
  rampho_morpho[[trait]]$analysis$initDGs <- list(
    #nullModel =  hzar.dataGroup.null(rampho_morpho[[trait]]$obs)
  )
}

# Create a model data group (hzar.dataGroup object) for each model from the initial runs.

for (trait in traits_morpho){
  #null_model_dg<-rampho_morpho[[trait]]$analysis$initDGs[[1]]
  #model_dgs<-list(null_model_dg)
  for (j in seq_along(rampho_morpho[[trait]]$models)){
    rampho_morpho[[trait]]$analysis$initDGs[[j]] <-
      hzar.dataGroup.add(rampho_morpho[[trait]]$runs$init[[j]])
   #model_dg <- hzar.dataGroup.add(rampho_morpho[[trait]]$runs$init[[j]])
   #model_dgs<-c(model_dgs, list(model_dg))
  }
  #rampho_morpho[[trait]]$analysis$initDGs <- model_dgs 
}

# name the models again
for (trait in traits_morpho) {
  names(rampho_morpho[[trait]]$analysis$initDGs) <- c(
    "model1",
    "model2",
    "model3",
    "model4",
    "model5")
}

# Create a hzar.obsDataGroup object from the five hzar.dataGroup
# just created, copying the naming scheme (model1, model 2, model3...).

for (trait in traits_morpho){
  rampho_morpho[[trait]]$analysis$oDG <-
    hzar.make.obsDataGroup(rampho_morpho[[trait]]$analysis$initDGs)
}


for (trait in traits_morpho){
  rampho_morpho[[trait]]$analysis$oDG <-
    hzar.copyModelLabels(rampho_morpho[[trait]]$analysis$initDGs,
                         rampho_morpho[[trait]]$analysis$oDG)
  
}

# Convert all runs to hzar.dataGroup objects, adding them to
# the hzar.obsDataGroup object.
for (trait in traits_morpho){
  rampho_morpho[[trait]]$analysis$oDG <-
    hzar.make.obsDataGroup(lapply(rampho_morpho[[trait]]$runs$chains,
                                  hzar.dataGroup.add),
                           rampho_morpho[[trait]]$analysis$oDG);
}

# Do model selection based on the AICc scores
for (trait in traits_morpho){
  print(rampho_morpho[[trait]]$analysis$AICcTable <- 
          hzar.AICc.hzar.obsDataGroup(rampho_morpho[[trait]]$analysis$oDG));
}

# Print out the model with the minimum AICc score
for (trait in traits_morpho){
  print(rampho_morpho[[trait]]$analysis$model.name <-
          rownames(rampho_morpho[[trait]]$analysis$AICcTable
          )[[ which.min(rampho_morpho[[trait]]$analysis$AICcTable$AICc )]])
}

# Extract the hzar.dataGroup object for the selected model
for (trait in traits_morpho){
  rampho_morpho[[trait]]$analysis$model.selected <-
    rampho_morpho[[trait]]$analysis$oDG$data.groups[[rampho_morpho[[trait]]$analysis$model.name]]
}

# Save the results in an Rdata object to plot later
rampho_morpho_MT2<-rampho_morpho
save(rampho_morpho_MT2, file="/scratch/juy3_lab/GBS_Ramphocelus/ClineAnalysis/T2/rampho_morpho_MT2.RData")

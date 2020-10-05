###########################################
#                                         #
#       Red-list diversity analysis       #
#     following Boccard et al 2018        #
#                                         #
#                  by                     #
#           Edouard Lavergne              #
#                 2020                    #
#                                         #
###########################################

 
# Diversity analyses in this script follow "Numerical Ecology with R"
# by Daniel Borcard, Francois Gillet and Pierre Legendre

# All analyses presented in the manuscript can be found in this script
# Four data tables (.csv files) are also included with this script in order
# to re-run the analyses. Do not forget to change the path of the folder where
# you saved the .csv data files using the function setwd() for the script to work
# Data files: "2020-06 - 22rivers_HT_final.csv"
#             "2020-06 - 22rivers_LT_final.csv"
#             "2019 - Species red list.csv"
#             "2020-06 - 22rivers_Mouth_Env_data.csv"
# R version 3.6.1 (2019-07-05)

# Load the required packages
library(BiodiversityR)
library(ClassDiscovery)
library(Hmisc)
library(plotrix)
library(pvclust)
library(vegan)

# Clear all objects in R
rm(list=ls())


#########################
# Load Custom Functions #
#########################

# The "red.list" function transforms the P/A community matrix into a red-list matrix
# where lines represent sites and columns represent red-list categories.
# X is a P/A community matrix,  where lines represent sites and columns represent species
# Y is a vector indicating the red-list category (DD, LC, NT, LP, VU, EN and CR in our case) of each species
red.list <- function(X, Y) {
  tmp <- matrix(data = 0, nrow = dim(X)[1], ncol = length(levels(Y)), byrow = FALSE)
  rownames(tmp) <- rownames(X)
  colnames(tmp) <- levels(Y)
  for (i in 1:length(levels(Y))){
    print(names(Y)[which(Y == levels(Y)[i])])
    tmp.1 <- names(Y)[which(Y == levels(Y)[i])]
    tmp.2 <- 0
    for (j in 1:length(tmp.1)){
      if (length(which(colnames(X)==tmp.1[j]))==0){
        tmp.2 <- c(tmp.2,"NA")
      } else {
        tmp.2 <- c(tmp.2, which(colnames(X)==tmp.1[j]))
      }
    }
    print(tmp.2)
    tmp.2 <- as.integer(tmp.2)
    tmp.2 <- tmp.2[-1]
    if (any(is.na(tmp.2))){
      tmp.2 <- tmp.2[-which(is.na(tmp.2))]
    }
    if (length(tmp.2)>1){
      tmp[,i] <- rowSums(X[,tmp.2])
    } else{
      tmp[,i] <- X[,tmp.2]
    }
  }
  X <- tmp
  X
}

#############################
# Setting working directory #
#############################

# Enter the path of the folder where you saved the .csv data files
setwd("D:/Edouard/Dropbox/[R] programs/Stat/eDNA/")

################
# Loading data #
################

# Load High Tide (HT) community data (here eDNA count data but could already be P/A data)
# Lines represent sites (rivers) and columns represent species
spe.ht <- read.table("2020-06 - 22rivers_HT_final.csv",
                     header=TRUE, sep=",", na.strings="NA",
                     dec=".", strip.white=TRUE)
row.names(spe.ht) <- spe.ht[,1] # Load row names from the first column
spe.ht <- spe.ht[,-1] # delete first column with names

# Load LT community data
spe.lt <- read.table("2020-06 - 22rivers_LT_final.csv",
                     header=TRUE, sep=",", na.strings="NA",
                     dec=".", strip.white=TRUE)
row.names(spe.lt) <- spe.lt[,1] # Load row names from the first column
spe.lt <- spe.lt[,-1] # delete first column with names


# Load species red-list status (here from the Japanese Ministry of Environment JMOE)
# Lines represent species and columns represent red-list classifications for different organizations
# The first column is for JMOE categories (DD, LC, NT, LP, VU, EN and CR)
# The second column is for IUCN categories (DD, LC, NT, VU, EN and CR)
# Here we will use only the JMOE list
spe.rl <- read.table("2019 - Species red list.csv",
                     header=TRUE, sep=",", na.strings="NA",
                     dec=".", strip.white=TRUE)
spe.rl[,1] <- gsub(" ", ".", spe.rl[,1], ignore.case = FALSE, perl = FALSE,
                   fixed = FALSE, useBytes = FALSE)
row.names(spe.rl) <- spe.rl[,1] # Load row names from the first column
spe.rl <- spe.rl[,-1] # delete first column with names

# Load environmental data
# Lines represent sites (rivers) and columns represent environmental, land use or population data
env <- read.table("2020-06 - 22rivers_Mouth_Env_data.csv", header=TRUE,
                  sep=",", na.strings="NA", dec=".",
                  strip.white=TRUE)

# TO load row names make sure the sites are arranged in the same order as spe.ht and spe.lt
row.names(env) <- env[,1] # Load row names from the first column
env <- env[,-1] # delete first column with names
env <- env[,-1] # delete the column containing ID as I did not used it here


# Check for rivers that have not been sampled for both
# HT and LT data sets
# Rivers for HT
(tmp.1 <- which(rowMeans(spe.ht, na.rm = FALSE, dims = 1)==0))
# Rivers for LT
(tmp.2 <- which(rowMeans(spe.lt, na.rm = FALSE, dims = 1)==0))

#################################################
# SST and salinity comparison between HT and LT #
#################################################

# SST
#####

# Test of normality an homocedasticity
shapiro.test(env$SST_HT)
shapiro.test(env$SST_LT)
bartlett.test(c(env$SST_HT, env$SST_LT) ~ c(rep("H", length(env$SST_HT)), rep("L", length(env$SST_HT))))

# Two sided t test is performed
t.test(env$SST_HT, env$SST_LT)

# Average
paste("Average HT SST:",
      paste(round(mean(env$SST_HT),1), round(sd(env$SST_HT),1), sep = " \u00B1 "),
      paste("(Mean", "SD)", sep = " \u00B1 "), sep = " ")
paste("Average LT SsT:",
      paste(round(mean(env$SST_LT),1), round(sd(env$SST_LT),1), sep = " \u00B1 "),
      paste("(Mean", "SD)", sep = " \u00B1 "), sep = " ")
paste("Combined data set, Average SST:",
      paste(round(mean(env$SST_Av),1), round(sd(env$SST_Av),1), sep = " \u00B1 "),
      paste("(Mean", "SD)", sep = " \u00B1 "), sep = " ")


# Salinity
###########

# Test of normality and homocedasticity
shapiro.test(env$Salinity_HT)
shapiro.test(env$Salinity_LT)
bartlett.test(c(env$Salinity_HT, env$Salinity_LT) ~ c(rep("H", length(env$Salinity_HT)), rep("L", length(env$Salinity_HT))))

# Two sided test Wilcoxon test is performed
wilcox.test(env$Salinity_HT, env$Salinity_LT, correct=TRUE) 

# Average
paste("Average HT salinity:",
      round(mean(env$Salinity_HT),1),
      sep = " ")
paste("LT salinity:",
      round(mean(env$Salinity_LT),1),
      sep = " ")
paste("Combined data set, Average salinity:",
      round(mean(env$Salinity_Av),1),
      sep = " ")

# CONCLUSION: We decide to pool HT and LT data


#################################################
# Selecting environmental variables of interest #
#################################################

# SST and Salinity were both measured in situ at HT and LT, as no difference was found between
# HT and LT for both parameters, average SST and average salinity beween HT and LT were used
env.hl <- env[, c("Latitude", "Longitude", "SST_Av", "Salinity_Av", "SS", "Length_km",
                  "DENSITY2015", "PF2014", "AGRI2014", "FOREST2014",
                  "AB_LAND2014", "U_AREA2014", "RL2014",
                  "Ramsar", "GOLF2014", "W_AREA", "POP2015",
                  "pH", "CAR", "RAR", "CRAR", "DOB", "DOS",
                  "TN", "discharge")]

# Selectiing environmental variables (including High tide SST and High tide Salinity)
env.ht <- env[, c("Latitude", "Longitude", "SST_HT", "Salinity_HT", "SS", "Length_km",
                  "DENSITY2015", "PF2014", "AGRI2014", "FOREST2014",
                  "AB_LAND2014", "U_AREA2014", "RL2014",
                  "Ramsar", "GOLF2014", "W_AREA", "POP2015",
                  "pH", "CAR", "RAR", "CRAR", "DOB", "DOS",
                  "TN", "discharge")]

# Selectiing environmental variables (including Low tide SST and Low tide Salinity)
env.lt <- env[, c("Latitude", "Longitude", "SST_LT", "Salinity_LT", "SS", "Length_km",
                  "DENSITY2015", "PF2014", "AGRI2014", "FOREST2014",
                  "AB_LAND2014", "U_AREA2014", "RL2014",
                  "Ramsar", "GOLF2014", "W_AREA", "POP2015",
                  "pH", "CAR", "RAR", "CRAR", "DOB", "DOS",
                  "TN", "discharge")]



############################################
# Descriptive statistics of community data #
############################################

# Combining HT and LT community data set 
#######################################

# Fusion of HT and LT data
spe <- spe.ht + spe.lt

# Our species list is part of a larger study and includes the names of species were detected
# in other sites (i.e. coastal areas) than the ones relevent for the present study.
# We first remove all species entry that were not detected in the sites of interest for this study.
# Every lines (species) with 0 detection at all sites will be discarded.
# If a matrix presents no line with only 0 entries, then nothing will happen. 

# Removing from the list species not present in the sites of interest
(tmp <- which(colMeans(spe, na.rm = FALSE, dims = 1)==0))
if (length(tmp)>0){
  spe <- spe[,-tmp]
}


# eDNA results present reading counts for each species per sites
# Lets transform the data to presence/absence data (P/A = 1/0)
spe.pa <- decostand(spe, "pa")
paste("Total number of species:", ncol(spe.pa), sep = " ")
paste("Number of species present at all rivers:", length(which(colSums(spe.pa) == nrow(spe.pa))), sep = " ")
if (length(which(colSums(spe.pa) == nrow(spe.pa)))>0){
  paste("Species names:")
  which(colSums(spe.pa) == nrow(spe.pa))
}
paste("Number of species present in at least half of the rivers:", length(which(colSums(spe.pa) > nrow(spe.pa)/2)), sep = " ")
if (length(which(colSums(spe.pa) > nrow(spe.pa)/2))>0){
  paste("Species names:")
  which(colSums(spe.pa) > nrow(spe.pa)/2)
}
paste("Number of species present in only a single river:", length(which(colSums(spe.pa) == 1)), sep = " ")
if (length(which(colSums(spe.pa) == 1))>0){
  paste("Species names:")
  which(colSums(spe.pa) == 1)
}

# Species richness per river
############################

(s.rich <- rowSums(spe.pa))

# Test of normality
shapiro.test(s.rich)

# Range of species richness
range(s.rich)

# Average of species richness per site
paste("Average species richness per site:",
      paste(round(mean(s.rich), 1), round(sd(s.rich),1), sep = " \u00B1 "),
      paste("(Mean", "SD)", sep = " \u00B1 "), sep = " ")


# Red-list Species richness per site (NT, LP, VU, EN, CR)
#########################################################

(s.rich.rl <- rowSums(spe.pa[,rl.count]))

# Test of normality
shapiro.test(s.rich.rl)

# Red list Range of species richness per site
range(s.rich.rl)

# Average of Red list Species richness per site
paste("Average Red-listed species richness per site:",
      round(mean(s.rich.rl), 1),
      sep = " ")

# How many red-listed species
spe.only.rl <- spe.rl[which(spe.rl$JMOE != "DD" & spe.rl$JMOE != "LC"),]
rl.count <- 0
for (i in 1:nrow(spe.only.rl)){
  tmp <- which(colnames(spe.pa) == rownames(spe.only.rl)[i])
  rl.count <- c(rl.count, tmp)
}
rl.count <- rl.count[-1]
paste("Number of Red-listed species:", length(rl.count), sep = " ")

# How many red-listed species appear at one site only
rl.single.count <- 0
for (i in 1:length(which(colSums(spe.pa[,rl.count]) == 1))){
  tmp <- which(rownames(spe.only.rl)==names(which(colSums(spe.pa[, rl.count])==1))[i])
  rl.single.count <- c(rl.single.count, tmp)
}
rl.single.count <- rl.single.count[-1]
paste("Number of Red-listed species at one site only:", length(rl.single.count), sep = " ")

# Red-listed species occurence
range(colSums(spe.pa[,rl.count]))
mean(colSums(spe.pa[,rl.count]))
colSums(spe.pa[,rl.count])[which(colSums(spe.pa[,rl.count])>10)] # occurence > 10

# HT community data 
###################

# Removing from the list species not present in the sites of interest
# and species present only at LT
(tmp <- which(colMeans(spe.ht, na.rm = FALSE, dims = 1)==0))
if (length(tmp)>0){
  spe.ht <- spe.ht[,-tmp]
}

# eDNA results present reading counts for each species per sites
# Lets transform the data to presence/absence data (P/A = 1/0)
spe.pa.ht <- decostand(spe.ht, "pa")

paste("Total number of species at HT:", ncol(spe.pa.ht), sep = " ")
paste("Number of species present in all rivers at HT:",
      length(which(colSums(spe.pa.ht) == nrow(spe.pa.ht))), sep = " ")
if (length(which(colSums(spe.pa.ht) == nrow(spe.pa.ht)))>0){
  paste("Species names:")
  which(colSums(spe.pa.ht) == nrow(spe.pa.ht))
}
paste("Number of species present in at least half of the rivers at HT:",
      length(which(colSums(spe.pa.ht) > nrow(spe.pa.ht)/2)), sep = " ")
if (length(which(colSums(spe.pa.ht) > nrow(spe.pa.ht)/2))>0){
  paste("Species names:")
  which(colSums(spe.pa.ht) > nrow(spe.pa.ht)/2)
}
paste("Number of species present in only a single river at HT:",
      length(which(colSums(spe.pa.ht) == 1)), sep = " ")
if (length(which(colSums(spe.pa.ht) == 1))>0){
  paste("Species names:")
  which(colSums(spe.pa.ht) == 1)
}


# Species richness per site at HT
##################################

(s.rich.ht <- rowSums(spe.pa.ht))

# Test of normality
shapiro.test(s.rich.ht)

# Range of species richness per site at HT
range(s.rich.ht)

# Average
paste("Average species richness per sites at HT:",
      paste(round(mean(s.rich.ht), 1), round(sd(s.rich.ht),1), sep = " \u00B1 "),
      paste("(Mean", "SD)", sep = " \u00B1 "), sep = " ")

# How many red-listed species at HT
rl.count.ht <- 0
for (i in 1:nrow(spe.only.rl)){
  tmp <- which(colnames(spe.pa.ht) == rownames(spe.only.rl)[i])
  rl.count.ht <- c(rl.count.ht, tmp)
}
rl.count.ht <- rl.count.ht[-1]
paste("Number of Red-listed species at HT:", length(rl.count.ht), sep = " ")

# How many red-listed species appear in only one site at HT
rl.single.count.ht <- 0
for (i in 1:length(which(colSums(spe.pa.ht[,rl.count.ht]) == 1))){
  tmp <- which(rownames(spe.only.rl)==names(which(colSums(spe.pa.ht[, rl.count.ht])==1))[i])
  rl.single.count.ht <- c(rl.single.count.ht, tmp)
}
rl.single.count.ht <- rl.single.count.ht[-1]
paste("Number of Red-listed species in only one site at HT:", length(rl.single.count.ht), sep = " ")


# Red listed species occurences at HT
range(colSums(spe.pa.ht[,rl.count.ht]))
mean(colSums(spe.pa.ht[,rl.count.ht]))

# Red-list Species richness per site at HT
(s.rich.rl.ht <- rowSums(spe.pa.ht[,rl.count.ht]))

# Test of normality
shapiro.test(s.rich.rl.ht)

# Red list Range of species richness per site at HT
range(s.rich.rl.ht)

# Average of Red list Species richness per site at HT
paste("Average Red-listed species richness per site at HT:",
      paste(round(mean(s.rich.rl.ht), 1), round(sd(s.rich.rl.ht),1), sep = " \u00B1 "),
      paste("(Mean", "SD)", sep = " \u00B1 "), sep = " ")


# LT community data 
###################

# Removing from the list species not present in the sites of interest
# and species present only at HT
(tmp <- which(colMeans(spe.lt, na.rm = FALSE, dims = 1)==0))
if (length(tmp)>0){
  spe.lt <- spe.lt[,-tmp]
}

# eDNA results present reading counts for each species per sites
# Lets transform the data to presence/absence data (P/A = 1/0)
spe.pa.lt <- decostand(spe.lt, "pa")

paste("Total number of species at LT:", ncol(spe.pa.lt), sep = " ")
paste("Number of species present in all rivers at LT:",
      length(which(colSums(spe.pa.lt) == nrow(spe.pa.lt))), sep = " ")
if (length(which(colSums(spe.pa.lt) == nrow(spe.pa.lt)))>0){
  paste("Species names:")
  which(colSums(spe.pa.lt) == nrow(spe.pa.lt))
}
paste("Number of species present in at least half of the rivers at LT:",
      length(which(colSums(spe.pa.lt) > nrow(spe.pa.lt)/2)), sep = " ")
if (length(which(colSums(spe.pa.lt) > nrow(spe.pa.lt)/2))>0){
  paste("Species names:")
  which(colSums(spe.pa.lt) > nrow(spe.pa.lt)/2)
}
paste("Number of species present in only a single river at LT:",
      length(which(colSums(spe.pa.lt) == 1)), sep = " ")
if (length(which(colSums(spe.pa.lt) == 1))>0){
  paste("Species names:")
  which(colSums(spe.pa.lt) == 1)
}


# Species richness per site at LT
##################################

(s.rich.lt <- rowSums(spe.pa.lt))

# Test of normality
shapiro.test(s.rich.lt)

# Range of species richness per site at LT
range(s.rich.lt)

# Average
paste("Average species richness per sites at LT:",
      paste(round(mean(s.rich.lt), 1), round(sd(s.rich.lt),1), sep = " \u00B1 "),
      paste("(Mean", "SD)", sep = " \u00B1 "), sep = " ")

# How many red-listed species at LT
rl.count.lt <- 0
for (i in 1:nrow(spe.only.rl)){
  tmp <- which(colnames(spe.pa.lt) == rownames(spe.only.rl)[i])
  rl.count.lt <- c(rl.count.lt, tmp)
}
rl.count.lt <- rl.count.lt[-1]
paste("Number of Red-listed species at LT:",
      length(rl.count.lt)
      , sep = " ")

# How many red listed species appear in only one site at LT
rl.single.count.lt <- 0
for (i in 1:length(which(colSums(spe.pa.lt[,rl.count.lt]) == 1))){
  tmp <- which(rownames(spe.only.rl)==names(which(colSums(spe.pa.lt[, rl.count.lt])==1))[i])
  rl.single.count.lt <- c(rl.single.count.lt, tmp)
}
rl.single.count.lt <- rl.single.count.lt[-1]
paste("Number of Red-listed species in only one site at LT:",
      length(rl.single.count.lt), sep = " ")


# Red listed species occurences at LT
range(colSums(spe.pa.lt[,rl.count.lt]))
mean(colSums(spe.pa.lt[,rl.count.lt]))

# Red-list Species richness per site at LT
(s.rich.rl.lt <- rowSums(spe.pa.lt[,rl.count.lt]))

# Test of normality
shapiro.test(s.rich.rl.lt)

# Red list Range of species richness per site at LT
range(s.rich.rl.lt)

# Average of Red list Species richness per site at LT
paste("Average Red-listed species richness per site at LT:",
      paste(round(mean(s.rich.rl.lt), 1), round(sd(s.rich.rl.lt),1), sep = " \u00B1 "),
      paste("(Mean", "SD)", sep = " \u00B1 "), sep = " ")


#################################################
# Species richness comparison between HT and LT #
#################################################

# Species richness
##################

# Test of normality an homocedasticity
shapiro.test(s.rich.ht)
shapiro.test(s.rich.lt)
bartlett.test(c(s.rich.ht, s.rich.lt) ~ c(rep("H", length(s.rich.ht)), rep("L", length(s.rich.lt))))

# Two sided t test is performed
t.test(s.rich.ht, s.rich.lt)

# Red-listed species richness
#############################

# Test of normality an homocedasticity
shapiro.test(s.rich.rl.ht)
shapiro.test(s.rich.rl.lt)
bartlett.test(c(s.rich.rl.ht, s.rich.rl.lt) ~ c(rep("H", length(s.rich.rl.ht)), rep("L", length(s.rich.rl.lt))))

# Two sided t test is performed
t.test(s.rich.rl.ht, s.rich.rl.lt)

# CONCLUSION: No significant differences in species richness and red-listed species richness
# between high and low tides

#################################
# Producing the red List matrix #
#################################

spe.rl2 <- spe.rl[,1] # Keep only the list of interest (here the Japanese one JMOE
                      # not the IUCN one) 
names(spe.rl2) <- rownames(spe.rl)
spe.pa.rl <- red.list(spe.pa, spe.rl2)[,c(4,6,5,7,3,1)] # c.f. Custom Functions
                                                        # The order of the column (4,6,5,7,3,1)
                                                        # was rearranged from alphabetic order
                                                        # to threat level order as follow LC, NT, LP,
                                                        # VU, EN and CR. DD are not included in the analysis
                                                        # It will produce NAs warning messages
spe.pa.ht.rl <- red.list(spe.pa.ht, spe.rl2)[,c(4,6,5,7,3,1)]
spe.pa.lt.rl <- red.list(spe.pa.lt, spe.rl2)[,c(4,6,5,7,3,1)]

#################################################
# red list species comparison between HT and LT # NT, LP, VU, EN and CR
#################################################

# Test of normality an homocedasticity
shapiro.test(rowSums(spe.pa.rl[,c(2:6)]))
shapiro.test(rowSums(spe.pa.ht.rl[,c(2:6)]))
shapiro.test(rowSums(spe.pa.lt.rl[,c(2:6)]))
bartlett.test(c(rowSums(spe.pa.ht.rl[,c(2:6)]),
                rowSums(spe.pa.lt.rl[,c(2:6)])) ~ 
                c(rep("H", length(rowSums(spe.pa.ht.rl))),
                  rep("L", length(rowSums(spe.pa.lt.rl)))))

# Two sided t test is performed
t.test(rowSums(spe.pa.ht.rl[,c(2:6)]),
            rowSums(spe.pa.lt.rl[,c(2:6)]))

# Average red-listed species number per site (NT, LP, VU, EN and CR)
paste("Average Red-Listed species per river at HT:",
      paste(round(mean(rowSums(spe.pa.ht.rl[,c(2:6)])), 1),
            round(sd(rowSums(spe.pa.ht.rl[,c(2:6)])),1),
            sep = " \u00B1 "),
      paste("(Mean", "SD)", sep = " \u00B1 "), sep = " ")
paste("Average Red-Listed species per river at LT:",
      paste(round(mean(rowSums(spe.pa.lt.rl[,c(2:6)])), 1),
            round(sd(rowSums(spe.pa.lt.rl[,c(2:6)])),1),
            sep = " \u00B1 "),
      paste("(Mean", "SD)", sep = " \u00B1 "), sep = " ")
paste("Average Red-List species per river at both HT and LT:",
      round(mean(rowSums(spe.pa.rl[,c(2:6)])), 1),
      sep = " ")


###################################################
# Threatened species comparison between HT and LT # LP, VU, EN and CR
###################################################

# Test of normality an homocedasticity
shapiro.test(rowSums(spe.pa.rl[,c(3:6)]))
shapiro.test(rowSums(spe.pa.ht.rl[,c(3:6)]))
shapiro.test(rowSums(spe.pa.lt.rl[,c(3:6)]))
bartlett.test(c(rowSums(spe.pa.ht.rl[,c(3:6)]),
                rowSums(spe.pa.lt.rl[,c(3:6)])) ~ 
                c(rep("H", length(rowSums(spe.pa.ht.rl))),
                  rep("L", length(rowSums(spe.pa.lt.rl)))))
wilcox.test(rowSums(spe.pa.ht.rl[,c(3:6)]),
            rowSums(spe.pa.lt.rl[,c(3:6)]))

# Averagethreatened species number per site (LP, VU, EN and CR)
paste("Average Threatened species per river at HT:",
      round(mean(rowSums(spe.pa.ht.rl[,c(3:6)])), 1),
      sep = " ")
paste("Average Threatened species per river at LT:",
      paste(round(mean(rowSums(spe.pa.lt.rl[,c(3:6)])), 1),
            round(sd(rowSums(spe.pa.lt.rl[,c(3:6)])),1),
            sep = " \u00B1 "),
      paste("(Mean", "SD)", sep = " \u00B1 "), sep = " ")

paste("Average Threatened species per river at both HT and LT:",
      paste(round(mean(rowSums(spe.pa.rl[,c(3:6)])), 1),
            round(sd(rowSums(spe.pa.rl[,c(3:6)])),1),
            sep = " \u00B1 "),
      paste("(Mean", "SD)", sep = " \u00B1 "), sep = " ")

# Number species per red-list category
threat.cat <- levels(spe.rl$JMOE)[c(2,4,6,5,7,3,1)]
for (i in 1:length(threat.cat)){
  print(threat.cat[i])
  print(tmp <- length(which(spe.rl[colnames(spe.pa),1] == threat.cat[i])))
}

# Number species per red-list category at HT
for (i in 1:length(threat.cat)){
  print(threat.cat[i])
  print(tmp <- length(which(spe.rl[colnames(spe.pa.ht),1] == threat.cat[i])))
}

# Number species per red-list category at LT
for (i in 1:length(threat.cat)){
  print(threat.cat[i])
  print(tmp <- length(which(spe.rl[colnames(spe.pa.lt),1] == threat.cat[i])))
}


# Number of Red-Listed and threatened species
RL <- length(which(spe.rl[colnames(spe.pa),1] == "NT" | spe.rl[colnames(spe.pa),1] == "LP" | spe.rl[colnames(spe.pa),1] == "VU" | spe.rl[colnames(spe.pa),1] == "EN" | spe.rl[colnames(spe.pa),1] == "CR"))
threatened <- length(which(spe.rl[colnames(spe.pa),1] == "LP" | spe.rl[colnames(spe.pa),1] == "VU" | spe.rl[colnames(spe.pa),1] == "EN" | spe.rl[colnames(spe.pa),1] == "CR"))
RL.ht <- length(which(spe.rl[colnames(spe.pa.ht),1] == "NT" | spe.rl[colnames(spe.pa.ht),1] == "LP" | spe.rl[colnames(spe.pa.ht),1] == "VU" | spe.rl[colnames(spe.pa.ht),1] == "EN" | spe.rl[colnames(spe.pa.ht),1] == "CR"))
threatened.ht <- length(which(spe.rl[colnames(spe.pa.ht),1] == "LP" | spe.rl[colnames(spe.pa.ht),1] == "VU" | spe.rl[colnames(spe.pa.ht),1] == "EN" | spe.rl[colnames(spe.pa.ht),1] == "CR"))
RL.lt <- length(which(spe.rl[colnames(spe.pa.lt),1] == "NT" | spe.rl[colnames(spe.pa.lt),1] == "LP" | spe.rl[colnames(spe.pa.lt),1] == "VU" | spe.rl[colnames(spe.pa.lt),1] == "EN" | spe.rl[colnames(spe.pa.lt),1] == "CR"))
threatened.lt <- length(which(spe.rl[colnames(spe.pa.lt),1] == "LP" | spe.rl[colnames(spe.pa.lt),1] == "VU" | spe.rl[colnames(spe.pa.lt),1] == "EN" | spe.rl[colnames(spe.pa.lt),1] == "CR"))

paste("Number of red-listed and threatened species in the combined data set:", RL, "and", threatened, "species", sep = " ")
paste("Number of red-listed and threatened species at HT:", RL.ht, "and", threatened.ht, "species", sep = " ")
paste("Number of red-listed and threatened species at LT:", RL.lt, "and", threatened.lt, "species", sep = " ")

# Simplifing matrix names
X <- spe.pa.rl
Y <- env.hl

# Compute distance matrix among sites (Hellinger transformation)
RL.hel <- decostand(X, "hellinger")
RL.ch <- vegdist(RL.hel, "euc")

# Attach site names to object of class 'dist'
attr(RL.ch, "labels") <- rownames(X)

#Compute Ward's minimum variance clustering
RL.ch.ward <- hclust(RL.ch, method = "ward.D2")

# cophenetic coefficient
RL.ch.ward.coph <- cophenetic(RL.ch.ward)
cor(RL.ch, RL.ch.ward.coph, method = "spearman")

# Gower (1983) distance
(gow.dist.ward <- sum((RL.ch - RL.ch.ward.coph) ^ 2))

# Compute p-values for all clusters (edges) of a dendrogram
spech.pv <-  pvclust(t(RL.hel),
                     method.hclust = "ward.D2",
                     method.dist = "euc",
                     parallel = TRUE,
                     nboot=10000)

# Plot dendrogram with p-values
par(mfrow=c(1,1))
plot(spech.pv)

# Highlight clusters with high AU p-values
pvrect(spech.pv, alpha = 0.95, pv = "au")

# hclust class dendogram
spech.den <-hclust(RL.ch, method = "ward.D2")

# Iidentify river group
gr <- cutree(spech.den, k=2)
gr <- -1*(gr-3) # This was to invert the cluster numbers (i.e. 1 and 2)
                # to fit the cluster name used in teh manuscript 
gr2 <- gr # Groups backup


# Build the clustered dendogram
###############################

# Plot the final dendrogram with group colors according to 
# human-use intensity (need RDA results)
col1 <- c("firebrick2","dodgerblue2")
par(mfrow=c(1,1))
plotColoredClusters(RL.ch.ward, labs = rownames(X),
                    cols = col1[gr], cex = 0.8, main = "",
                    line = 0)

# Add columns with species number for Threaten species (VU, EN, CR) and Red-listes species
X2 <- cbind(X, rowSums(X[,3:6]), rowSums(X[,2:6]))
colnames(X2)[c(7,8)] <- c("TH", "RL")

# How many red-listed species (NT, LP, VU, EN, CR)
rl.count.t <- 0
for (i in 1:nrow(spe.only.rl)){
  tmp <- which(colnames(spe.pa) == rownames(spe.only.rl)[i])
  rl.count.t <- c(rl.count.t, tmp)
}
rl.count.t <- rl.count.t[-1]
paste("Number of red-listed species:", length(rl.count.t), sep = " ")

# How many red-listed species (NT, LP, VU, EN, CR) for cluster 1
rl.count.gr1 <- 0
for (i in 1:nrow(spe.only.rl)){
  tmp <- which(names(which(colSums(spe.pa[gr==1,])!=0)) == rownames(spe.only.rl)[i])
  rl.count.gr1 <- c(rl.count.gr1, tmp)
}
rl.count.gr1 <- rl.count.gr1[-1]
paste("Number of red-listed species in Cluster 1:", length(rl.count.gr1), sep = " ")

# How many red listed species (NT, LP, VU, EN, CR) for cluster 2
rl.count.gr2 <- 0
for (i in 1:nrow(spe.only.rl)){
  tmp <- which(names(which(colSums(spe.pa[gr==2,])!=0)) == rownames(spe.only.rl)[i])
  rl.count.gr2 <- c(rl.count.gr2, tmp)
}
rl.count.gr2 <- rl.count.gr2[-1]
paste("Number of red-listed species in Cluster 2:", length(rl.count.gr2), sep = " ")


# Red-listed species occurrences
range(colSums(spe.pa[, rl.count.t]))
mean(colSums(spe.pa[, rl.count.t]))

# HT
range(colSums(spe.pa.ht[, rl.count.ht]))
mean(colSums(spe.pa.ht[, rl.count.ht]))

# LT
range(colSums(spe.pa.lt[, rl.count.lt]))
mean(colSums(spe.pa.lt[, rl.count.lt]))

# All data Cluster 1
range(colSums(spe.pa[gr==1, rl.count.gr1]))
mean(colSums(spe.pa[gr==1, rl.count.gr1]))

# All data Cluster 2
range(colSums(spe.pa[gr==2, rl.count.gr2]))
mean(colSums(spe.pa[gr==2, rl.count.gr2]))


# Number of species per red-list category for each cluster group
spe.pa.clust1 <- spe.pa[gr==1, which(colSums(spe.pa[gr==1,])!=0)]
spe.pa.clust2 <- spe.pa[gr==2, which(colSums(spe.pa[gr==2,])!=0)]

# Cluster 1
for (i in 1:length(threat.cat)){
  print(threat.cat[i])
  print(tmp <- length(which(spe.rl[colnames(spe.pa.clust1),1] == threat.cat[i])))
}

# Cluster 2
for (i in 1:length(threat.cat)){
  print(threat.cat[i])
  print(tmp <- length(which(spe.rl[colnames(spe.pa.clust2),1] == threat.cat[i])))
}


# Environmental data
####################

# Calculating total area from land use data
# The Ramsar variable represents protected wetlands and this variable is omitted
# in the total area calculation because protected wetlands areas are part of other
# land use areas
tot.area <- rowSums(Y[,c("PF2014", "AGRI2014",
                              "FOREST2014", "AB_LAND2014",
                              "U_AREA2014", "RL2014", "GOLF2014")])

# Calculating the proportion of watershed area cover for each land use
Y[,c("PF2014", "AGRI2014",
          "FOREST2014", "AB_LAND2014",
          "U_AREA2014", "RL2014",
          "Ramsar", "GOLF2014")] <- Y[,c("PF2014", "AGRI2014",
                                          "FOREST2014", "AB_LAND2014",
                                          "U_AREA2014", "RL2014",
                                          "Ramsar", "GOLF2014")]*100/tot.area

# Env data correlations
library(corrgram)
corrgram(Y, order=NULL, lower.panel=panel.shade,
         upper.panel=NULL)


# Full data normality test
##########################

# raw data
for (i in 1:ncol(env)) {
  print("------------------------------------")
  print(colnames(env)[i])
  print(shapiro.test(env[,i]))
  print(mean(env[,i]))
  print(sd(env[,i]))
}

# Percentage data
for (i in 1:ncol(Y)) {
  print("------------------------------------")
  print(colnames(Y)[i])
  print(shapiro.test(Y[,i]))
  print(mean(Y[,i]))
  print(sd(Y[,i]))
}


# Cluster group comparisons
###########################

# Test assumptions for env data
# Normality of residuals and Homogeneity of variances
var.test <- rep(1, ncol(Y))
names(var.test) <- colnames(Y)
for (i in 1:ncol(Y)) {
  shap.tmp <- shapiro.test(resid(aov(Y[,i] ~ as.factor(gr))))
  bart.tmp <- bartlett.test(Y[,i] ~ as.factor(gr))
  # Select parametric or non parametric test
  if (shap.tmp$p.value <= 0.05 | bart.tmp$p.value <= 0.05){
    var.test[i] <- 2
  }
  print("-----------------------------------")
  print(colnames(Y)[i])
  print(shap.tmp)
  print(bart.tmp)
  print("-----------------------------------")
}

# H0 (two sided test) no difference in average env data between group 1 and 2 
# H0 (one sided test) average env data in group 1 is larger (or smaller) that the average in group 2 (one sided)
count.pval <- 0
var.mem <- 0
for (i in 1:ncol(Y)) {
  if (var.test[i] == 1){
    test.tmp <- t.test(Y[,i] ~ as.factor(gr), alternative=c("two.sided"))
    test.text <- "two sided test"
    test.av.gr1 <- round(apply(Y[gr==1,], 2, mean),3)[i]
    test.sd.gr1 <- round(apply(Y[gr==1,], 2, sd),3)[i]
    test.av.gr2 <- round(apply(Y[gr==2,], 2, mean),3)[i]
    test.sd.gr2 <- round(apply(Y[gr==2,], 2, sd),3)[i]
    if (test.tmp$p.value > 0.05 & test.tmp$p.value < 0.1) { # test for one sided test
      test.tmp2 <- t.test(Y[,i] ~ as.factor(gr), alternative=c("greater"))
      test.tmp3 <- t.test(Y[,i] ~ as.factor(gr), alternative=c("less"))
      if (test.tmp2$p.value <= 0.05){
        test.tmp <- test.tmp2
        test.text <- "one sided test"
        count.pval <- count.pval + 1
        var.mem <- cbind(var.mem, i)
        test.av.gr1 <- round(apply(Y[gr==1,], 2, mean),3)[i]
        test.sd.gr1 <- round(apply(Y[gr==1,], 2, sd),3)[i]
        test.av.gr2 <- round(apply(Y[gr==2,], 2, mean),3)[i]
        test.sd.gr2 <- round(apply(Y[gr==2,], 2, sd),3)[i]
      } else if (test.tmp3$p.value <= 0.05){
        test.tmp <- test.tmp3
        test.text <- "one sided test"
        count.pval <- count.pval + 1
        var.mem <- cbind(var.mem, i)
        test.av.gr1 <- round(apply(Y[gr==1,], 2, mean),3)[i]
        test.sd.gr1 <- round(apply(Y[gr==1,], 2, sd),3)[i]
        test.av.gr2 <- round(apply(Y[gr==2,], 2, mean),3)[i]
        test.sd.gr2 <- round(apply(Y[gr==2,], 2, sd),3)[i]
      }
    } else if (test.tmp$p.value <= 0.05) {
      count.pval <- count.pval + 1
      var.mem <- cbind(var.mem, i)
    }
  } else {
    test.tmp <- wilcox.test(Y[,i] ~ as.factor(gr), alternative=c("two.sided"))
    test.text <- "two sided test"
    test.av.gr1 <- round(apply(Y[gr==1,], 2, mean),3)[i]
    test.sd.gr1 <- round(apply(Y[gr==1,], 2, sd),3)[i]
    test.av.gr2 <- round(apply(Y[gr==2,], 2, mean),3)[i]
    test.sd.gr2 <- round(apply(Y[gr==2,], 2, sd),3)[i]
    if (test.tmp$p.value > 0.05 & test.tmp$p.value < 0.1) { # test for one sided test
      test.tmp2 <- wilcox.test(Y[,i] ~ as.factor(gr), alternative=c("greater"))
      test.tmp3 <- wilcox.test(Y[,i] ~ as.factor(gr), alternative=c("less"))
      if (test.tmp2$p.value <= 0.05){
        test.tmp <- test.tmp2
        test.text <- "one sided test"
        count.pval <- count.pval + 1
        var.mem <- cbind(var.mem, i)
        test.av.gr1 <- round(apply(Y[gr==1,], 2, mean),3)[i]
        test.sd.gr1 <- round(apply(Y[gr==1,], 2, sd),3)[i]
        test.av.gr2 <- round(apply(Y[gr==2,], 2, mean),3)[i]
        test.sd.gr2 <- round(apply(Y[gr==2,], 2, sd),3)[i]
      } else if (test.tmp3$p.value <= 0.05){
        test.tmp <- test.tmp3
        test.text <- "one sided test"
        count.pval <- count.pval + 1
        var.mem <- cbind(var.mem, i)
        test.av.gr1 <- round(apply(Y[gr==1,], 2, mean),3)[i]
        test.sd.gr1 <- round(apply(Y[gr==1,], 2, sd),3)[i]
        test.av.gr2 <- round(apply(Y[gr==2,], 2, mean),3)[i]
        test.sd.gr2 <- round(apply(Y[gr==2,], 2, sd),3)[i]
      }
    } else if (test.tmp$p.value <= 0.05) {
      count.pval <- count.pval + 1
      var.mem <- cbind(var.mem, i)
    }
  }
  print("-----------------------------------")
  print(colnames(Y)[i])
  print(test.text)
  print(test.tmp)
  print(test.av.gr1)
  print(test.sd.gr1)
  print(test.av.gr2)
  print(test.sd.gr2)
  print("-----------------------------------")
}
print(count.pval)
var.mem <- var.mem[-1]
print(var.mem)

# Plot boxplots
nc <- ceiling(length(var.mem)/2)
par(mfrow=c(2,nc))
for (i in 1:length(var.mem)){
  boxplot(Y[,var.mem[i]] ~ as.factor(gr),
          col = col1,
          main= colnames(Y)[var.mem[i]],
          ylab="Cover (%)")
}


# Group comparison for red-list categories
##########################################

# Normality of residuals and Homogeneity of variances
colnames(X2) <- c(colnames(X), "threatened", "red listed")
cat.test <- rep(1, ncol(X2))
for (i in 1:length(cat.test)) {
    shap.tmp <- shapiro.test(resid(aov(X2[,i] ~ as.factor(gr))))
    bart.tmp <- bartlett.test(X2[,i] ~ as.factor(gr))
    # Select parametric or non parametric test
    if (shap.tmp$p.value <= 0.05 | bart.tmp$p.value <= 0.05) {
      cat.test[i] <- 2
    }
  print("-----------------------------------")
  print(names(cat.test)[i])
  print(shap.tmp)
  print(bart.tmp)
  print("-----------------------------------")
}

# H0 (two sided test) no difference in average env data between group 1 and 2 
# H0 (one sided test) average env data in group 1 is larger (or smaller) that the average in group 2 
count.pval <- 0
cat.mem <- 0
for (i in 1:length(cat.test)) {
  if (cat.test[i] == 1){
    test.tmp <- t.test(X2[,i] ~ as.factor(gr), alternative=c("two.sided"))
    test.text <- "two sided test"
    if (test.tmp$p.value > 0.05 & test.tmp$p.value < 0.2) { # test for one sided test
      test.tmp2 <- t.test(X2[,i] ~ as.factor(gr), alternative=c("greater"))
      test.tmp3 <- t.test(X2[,i] ~ as.factor(gr), alternative=c("less"))
      if (test.tmp2$p.value <= 0.05){
        test.tmp <- test.tmp2
        test.text <- "one sided test"
        count.pval <- count.pval + 1
        cat.mem <- cbind(cat.mem, i)
      } else if (test.tmp3$p.value <= 0.05){
        test.tmp <- test.tmp3
        test.text <- "one sided test"
        count.pval <- count.pval + 1
        cat.mem <- cbind(cat.mem, i)
      }
    } else if (test.tmp$p.value <= 0.05) {
      count.pval <- count.pval + 1
      cat.mem <- cbind(cat.mem, i)
    }
  } else {
    test.tmp <- wilcox.test(X2[,i] ~ as.factor(gr), alternative=c("two.sided"))
    test.text <- "two sided test"
    if (test.tmp$p.value > 0.05 & test.tmp$p.value < 0.2) { # test for one sided test
      test.tmp2 <- wilcox.test(X2[,i] ~ as.factor(gr), alternative=c("greater"))
      test.tmp3 <- wilcox.test(X2[,i] ~ as.factor(gr), alternative=c("less"))
      if (test.tmp2$p.value <= 0.05){
        test.tmp <- test.tmp2
        test.text <- "one sided test"
        count.pval <- count.pval + 1
        cat.mem <- cbind(cat.mem, i)
      } else if (test.tmp3$p.value <= 0.05){
        test.tmp <- test.tmp3
        test.text <- "one sided test"
        count.pval <- count.pval + 1
        cat.mem <- cbind(cat.mem, i)
      }
    } else if (test.tmp$p.value <= 0.05) {
      count.pval <- count.pval + 1
      cat.mem <- cbind(cat.mem, i)
    }
  }
  print("-----------------------------------")
  print(colnames(X2)[i])
  print(test.text)
  print(test.tmp)
  
  print("Total")
  print(mean(X2[, i]))
  print(sd(X2[, i]))
  
  print("Cluster 1")
  print(mean(X2[gr==1, i]))
  print(sd(X2[gr==1, i]))
  print(range(X2[gr==1, i]))
  
  print("Cluster 2")
  print(mean(X2[gr==2, i]))
  print(sd(X2[gr==2, i]))
  print(range(X2[gr==2, i]))
  
  print("-----------------------------------")
}
print(count.pval)
cat.mem <- cat.mem[-1]
print(cat.mem)

# Check for differences in species richness between clusters
shapiro.test(resid(aov(s.rich ~ as.factor(gr))))
bartlett.test(s.rich ~ as.factor(gr))
t.test(s.rich ~ as.factor(gr), alternative=c("two.sided"))
mean(s.rich)
sd(s.rich)
mean(s.rich[gr==1])
sd(s.rich[gr==1])
mean(s.rich[gr==2])
sd(s.rich[gr==2])

print("Species richness per Cluster")
length(which(colSums(spe.pa) > 0))
length(which(colSums(spe.pa[gr==1,]) > 0))
length(which(colSums(spe.pa[gr==2,]) > 0))


# Plot boxplots
par(mfrow=c(2,4))
for (i in 1:length(cat.test)){
  boxplot(X2[,i] ~ as.factor(gr),
          col = col1,
          main= colnames(X2)[i],
          ylab="Number of species")
}


# Number of sites containing species of a particular red list category
site.cat.sen.th <- matrix(data = NA, nrow = 2, ncol = 8)
colnames(site.cat.sen.th) <- colnames(X2)
rownames(site.cat.sen.th) <- c("gr1", "gr2")
tmp1 <- X2[which(gr==1),]
tmp2 <- X2[which(gr==2),]
for(i in 1:dim(X2)[2]){
  site.cat.sen.th[1,i] <- dim(tmp1)[1] - length(which(tmp1[,i]==0))
  site.cat.sen.th[2,i] <- dim(tmp2)[1] - length(which(tmp2[,i]==0))
}
site.cat.sen.th

# Linear relations between Forest and Red-listed species number
fit <- lm((X2[,8]) ~ Y$FOREST2014)
summary(fit)
summary1 = paste0("Red listed species = ", round(fit$coefficients[2],5), " x Forest % + ",round(fit$coefficients[1],2))
summary2 = expression(paste("Adjusted  ", R^2, " =", sep = " "))

par(mfrow=c(1,1))
plot(Y$FOREST2014, X2[,8], pch=19,
     col= col1[gr], xlim = c(40,95), ylim=c(0,12),
     axes=FALSE, ann = FALSE)
text(65, 11, summary1, pos=2)
text(62, 10, summary2, pos=2)
text(65, 10, round(summary(fit)$adj.r.squared, 2), pos=2)
mtext("Forest (%)", side=1, line=1.4, cex=1)
mtext("Number of red listed species", side=2, line=2, cex=1)
mgp.axis(1, mgp = c(0, 0.4, 0), cex.axis=1)
mgp.axis(2, mgp = c(0, 0.6, 0), cex.axis=1, las=2)
abline(fit, col="black")
F.lim <- (summary(rowSums(X[,c(2:6)]))[4] - fit$coefficients[1]) / fit$coefficients[2]
F.lim
abline(v = F.lim, lty = 3, col = "black")
text(Y$FOREST2014, X2[,8], rownames(X2),cex=0.8, pos=1)
box()


# Linear relations between urban area and Red-listed species number
fit <- lm(X2[,8] ~ Y$U_AREA2014)
summary(fit)
summary1 = paste0("Red liested species = ", round(fit$coefficients[2],5), " x Urban % + ",round(fit$coefficients[1],2))
summary2 = expression(paste("Adjusted  ", R^2, " =", sep = " "))

par(mfrow=c(1,1))
plot(Y$U_AREA2014, X2[,8], pch=19,
     col= col1[gr], xlim = c(0,40), ylim=c(0,12),
     axes=FALSE, ann = FALSE)
text(40, 12, summary1, pos=2)
text(38, 11, summary2, pos=2)
text(40, 11, round(summary(fit)$adj.r.squared, 2), pos=2)
mtext("Urban area (%)", side=1, line=1.4, cex=1)
mtext("Number of red listed species", side=2, line=2, cex=1)
mgp.axis(1, mgp = c(0, 0.4, 0), cex.axis=1)
mgp.axis(2, mgp = c(0, 0.6, 0), cex.axis=1, las=2)
abline(fit, col="black")
F.lim <- (summary(rowSums(X[,c(2:6)]))[4] - fit$coefficients[1]) / fit$coefficients[2]
F.lim
abline(v = F.lim, lty = 3, col = "black")
text(Y$U_AREA2014, X2[,8], rownames(X2),cex=0.8, pos=1)
box()

# Linear relations between abandoned land and Red-listed species number
fit <- lm(X2[,8] ~ Y$AB_LAND2014)
summary(fit)
summary1 = paste0("Red liested species = ", round(fit$coefficients[2],5), " x Ab land + ",round(fit$coefficients[1],2))
summary2 = expression(paste("Adjusted  ", R^2, " =", sep = " "))

par(mfrow=c(1,1))
plot(Y$AB_LAND2014, X2[,8], pch=19,
     col= col1[gr], xlim = c(0,10), ylim=c(0,12),
     axes=FALSE, ann = FALSE)
text(10, 12, summary1, pos=2)
text(9.5, 11, summary2, pos=2)
text(10, 11, round(summary(fit)$adj.r.squared, 2), pos=2)
mtext("ab land (%)", side=1, line=1.4, cex=1)
mtext("Number of red listed species", side=2, line=2, cex=1)
mgp.axis(1, mgp = c(0, 0.4, 0), cex.axis=1)
mgp.axis(2, mgp = c(0, 0.6, 0), cex.axis=1, las=2)
abline(fit, col="black")
F.lim <- (summary(rowSums(X[,c(2:6)]))[4] - fit$coefficients[1]) / fit$coefficients[2]
F.lim
abline(v = F.lim, lty = 3, col = "black")
box()

# Linear relations between pop and Red-listed species number
tmp <- Y$POP2015/1000
fit <- lm(X2[,8] ~ tmp)
summary(fit)
summary1 = paste0("Red liested species = ", round(fit$coefficients[2],5), " x pop + ",round(fit$coefficients[1],2))
summary2 = expression(paste("Adjusted  ", R^2, " =", sep = " "))

par(mfrow=c(1,1))
plot((Y$POP2015/1000), X2[,8], pch=19,
     col= col1[gr], xlim = c(0,2500), ylim=c(0,12),
     axes=FALSE, ann = FALSE)
text(2500, 12, summary1, pos=2)
text(2400, 11, summary2, pos=2)
text(2500, 11, round(summary(fit)$adj.r.squared, 2), pos=2)
mtext("pop", side=1, line=1.4, cex=1)
mtext("Number of red listed species", side=2, line=2, cex=1)
mgp.axis(1, mgp = c(0, 0.4, 0), cex.axis=1)
mgp.axis(2, mgp = c(0, 0.6, 0), cex.axis=1, las=2)
abline(fit, col="black")
F.lim <- (summary(rowSums(X[,c(2:6)]))[4] - fit$coefficients[1]) / fit$coefficients[2]
F.lim
abline(v = F.lim, lty = 3, col = "black")
box()


# Linear relations between Latitude and Red-listed species number
fit <- lm(X2[,8] ~ Y$Latitude)
summary(fit)

# Linear relations between longitude and Red-listed species number
fit <- lm(X2[,8] ~ Y$Longitude)
summary(fit)


# PCA analysis
##############

RL.h.pca <- rda(RL.hel)

# Biplots PCA scaling 1
pca.sum <- summary(RL.h.pca, scaling = 1)
(RL.h.pca.env <- envfit(pca.sum ~ ., data=Y[,-c(19,20)], perm = 9999))
RL.h.pca.env2 <- RL.h.pca.env # back up
par(mfrow = c(1, 1), mar=c(4,4,1,1))
explained.var <- round(pca.sum$cont$importance[2,],3)*100
cum.explained.var <- round(pca.sum$cont$importance[3,],3)*100
text1 <- paste("PCA 1:", explained.var[1],
               "% of the total explained variance", sep = " ")
text2 <- paste("PCA 2:", explained.var[2],
               "% of the total explained variance", sep = " ")
plot(pca.sum$sites, pch=19, col=col1[gr],
     xlab = text1, ylab = text2,
     asp = 1, axes = TRUE, main = "PCA scaling 1",
     xlim=c(-0.6,0.6), ylim=c(-0.6,0.6))
text(pca.sum$sites, rownames(pca.sum$sites),
     col=col1[gr],cex = 0.7, pos=4)
arrows(pca.sum$species[,1]*0, pca.sum$species[,2]*0,
       pca.sum$species[,1], pca.sum$species[,2],
       col= "black", length = 0)
text(pca.sum$species, rownames(pca.sum$species), col="black",
     cex = 1, pos=4)

# Plot significant variables with a user-selected colour
p.symb.1 <- which(RL.h.pca.env$vectors$pvals <= 0.001) # ***
p.symb.2 <- which(RL.h.pca.env$vectors$pvals > 0.001 & RL.h.pca.env$vectors$pvals <= 0.01) # **
p.symb.3 <- which(RL.h.pca.env$vectors$pvals > 0.01 & RL.h.pca.env$vectors$pvals <= 0.05) # *
p.symb.4 <- which(RL.h.pca.env$vectors$pvals > 0.05 & RL.h.pca.env$vectors$pvals <= 0.1) # .
symb <- rep("", length(RL.h.pca.env$vectors$pvals))
symb[p.symb.1] <- "***"
symb[p.symb.2] <- "**"
symb[p.symb.3] <- "*"
symb[p.symb.4] <- "."
rownames(RL.h.pca.env$vectors$arrows) <- paste(rownames(RL.h.pca.env2$vectors$arrows), symb, sep=" ")
plot(RL.h.pca.env, p.max = 0.05, col = "forestgreen")
abline(h = 0, lty = 3)
abline(v = 0, lty = 3)

# Equilibrium circle: sqrt(d/p)
d <- 2 # dimension of the reduced space (usually, d = 2)
p <-  length(which(pca.sum$cont$importance[1,]>0)) # dimensionality of the multivariate PCA
                                                    # space, which is the number of eigenvalues > 0;
                                                    # usually equal to the number of descriptors.
eq.radius <- sqrt(d/p) # Radius of equilibrium circle

# Draw circle
draw.circle(0,0, eq.radius, nv=100, border="black", lty=3)


# Broken stick analysis
(PCA.sig <- PCAsignificance(RL.h.pca,axes=6))

# Kaiser Guttman criterion and broken stick model
bp <- barplot(PCA.sig[2,], xlab = "PCA axes", ylab = "Relative Eigenvalues")
abline(h=mean(PCA.sig[2,]), lty=2, lwd=2, col="blue")
b_stick <- PCA.sig[4,]
points(bp, b_stick, pch=19, col="red")
lines(bp, b_stick, type="l", lwd=2, lty=3, col="red")


# Biplots PCA scaling 2
pca.sum2 <- summary(RL.h.pca, scaling = 2)
(RL.h.pca.env <- envfit(pca.sum2 ~ ., data=Y[,-c(19,20)], perm = 9999))
par(mfrow = c(1, 1), mar=c(4,4,1,1))
explained.var <- round(pca.sum2$cont$importance[2,],3)*100
cum.explained.var <- round(pca.sum2$cont$importance[3,],3)*100
text1 <- paste("PCA 1:", explained.var[1],
               "% of the total explained variance", sep = " ")
text2 <- paste("PCA 2:", explained.var[2],
               "% of the total explained variance", sep = " ")
plot(pca.sum2$sites, pch=19, col=col1[gr],
     xlab = text1, ylab = text2,
     asp = 1, axes = TRUE, main = "PCA scaling 2",
     xlim=c(-0.3,0.3), ylim=c(-0.4,0.6))
text(pca.sum2$sites, rownames(pca.sum2$sites),
     col=col1[gr],cex = 0.7, pos=4)
arrows(pca.sum2$species[,1]*0, pca.sum2$species[,2]*0,
       pca.sum2$species[,1], pca.sum2$species[,2],
       col= "black", length = 0)
text(pca.sum2$species, rownames(pca.sum2$species), col="black",
     cex = 1, pos=4)

# Plot significant variables with a user-selected colour
p.symb.1 <- which(RL.h.pca.env$vectors$pvals <= 0.001) # ***
p.symb.2 <- which(RL.h.pca.env$vectors$pvals > 0.001 & RL.h.pca.env$vectors$pvals <= 0.01) # **
p.symb.3 <- which(RL.h.pca.env$vectors$pvals > 0.01 & RL.h.pca.env$vectors$pvals <= 0.05) # *
p.symb.4 <- which(RL.h.pca.env$vectors$pvals > 0.05 & RL.h.pca.env$vectors$pvals <= 0.1) # .

symb <- rep("", length(RL.h.pca.env$vectors$pvals))
symb[p.symb.1] <- "***"
symb[p.symb.2] <- "**"
symb[p.symb.3] <- "*"
symb[p.symb.4] <- "."

rownames(RL.h.pca.env$vectors$arrows) <- paste(rownames(RL.h.pca.env2$vectors$arrows), symb, sep=" ")

plot(RL.h.pca.env, p.max = 0.05, col = "forestgreen")
abline(h = 0, lty = 3)
abline(v = 0, lty = 3)


## Same analysis including Rate of change (rc)
###############################################

# Rate of change for the different variables , the rc.period indicate the length of the period in years
rc.env.names <- c("POP15", "POP10", "POP5", "DENSITY15", "DENSITY10", "DENSITY5",
  "PF17", "PF8", "PF5", "AGRI17", "AGRI8", "AGRI5", "FOREST17",
  "FOREST8", "FOREST5", "AB_LAND17", "AB_LAND8", "AB_LAND5",
  "U_AREA17", "U_AREA8", "U_AREA5", "RL17", "RL8", "RL5",
  "GOLF17", "GOLF8", "GOLF5")
rc.period <- c(15,10,5,15,10,5,17,8,5,17,8,5,17,8,5,17,8,5,17,8,5,17,8,5,17,8,5)

rc.env <- matrix(0, nrow = nrow(env.hl), ncol = length(rc.env.names))
rownames(rc.env) <- row.names(env.hl)
colnames(rc.env) <-rc.env.names

tmp.env <- env[,-c(1:20)]
for (i in c(1:9)){
  rc.env[,(i*3)-2] <- ((tmp.env[,i*4]-tmp.env[,(i*4)-3])/tmp.env[,(i*4)-3])/rc.period[(i*3)-2]
  rc.env[,(i*3)-1] <- ((tmp.env[,i*4]-tmp.env[,(i*4)-2])/tmp.env[,(i*4)-2])/rc.period[(i*3)-2]
  rc.env[,(i*3)] <- ((tmp.env[,i*4]-tmp.env[,(i*4)-1])/tmp.env[,(i*4)-1])/rc.period[(i*3)]
} 
rc.env <- data.frame(rc.env)
# Because of some missing data some rate of change calculations return "inf"
# that we need to change to "NaN".
rc.env[7,25] <- NaN
rc.env[11,26] <- NaN

# back up
rc.env2 <- rc.env
# combine the previous data set with this rate of change data set  
rc.env <- cbind(Y[,-c(19,20)], rc.env2)

# Biplots PCA scaling 1
pca.sum3 <- summary(RL.h.pca, scaling = 1)
(pca.rc.env <- envfit(pca.sum3 ~ ., data = rc.env, perm = 9999, na.rm = TRUE))
par(mfrow = c(1, 1), mar=c(4,4,1,1))
explained.var <- round(pca.sum3$cont$importance[2,],3)*100
cum.explained.var <- round(pca.sum3$cont$importance[3,],3)*100
text1 <- paste("PCA 1:", explained.var[1],
               "% of the total explained variance", sep = " ")
text2 <- paste("PCA 2:", explained.var[2],
               "% of the total explained variance", sep = " ")
plot(pca.sum3$sites, pch=19, col=col1[gr],
     xlab = text1, ylab = text2,
     asp = 1, axes = TRUE, main = "PCA scaling 1",
     xlim=c(-0.8,0.8), ylim=c(-0.4,0.6))
text(pca.sum3$sites, rownames(pca.sum3$sites),
     col=col1[gr],cex = 0.7, pos=4)
arrows(pca.sum3$species[,1]*0, pca.sum3$species[,2]*0,
       pca.sum3$species[,1], pca.sum3$species[,2],
       col= "black", length = 0)
text(pca.sum3$species, rownames(pca.sum3$species), col="black",
     cex = 1, pos=4)

# Plot significant variables with a user-selected colour
p.symb.1 <- which(pca.rc.env$vectors$pvals <= 0.001) # ***
p.symb.2 <- which(pca.rc.env$vectors$pvals > 0.001 & pca.rc.env$vectors$pvals <= 0.01) # **
p.symb.3 <- which(pca.rc.env$vectors$pvals > 0.01 & pca.rc.env$vectors$pvals <= 0.05) # *
p.symb.4 <- which(pca.rc.env$vectors$pvals > 0.05 & pca.rc.env$vectors$pvals <= 0.1) # .

symb <- rep("", length(pca.rc.env$vectors$pvals))
symb[p.symb.1] <- "***"
symb[p.symb.2] <- "**"
symb[p.symb.3] <- "*"
symb[p.symb.4] <- "."

pca.rc.env2 <- pca.rc.env
rownames(pca.rc.env$vectors$arrows) <- paste(rownames(pca.rc.env2$vectors$arrows), symb, sep=" ")

plot(pca.rc.env, p.max = 0.05, col = "forestgreen")
abline(h = 0, lty = 3)
abline(v = 0, lty = 3)

# Equilibrium circle: sqrt(d/p)
d <- 2 # dimension of the reduced space (usually, d = 2)
p <-  length(which(pca.sum3$cont$importance[1,]>0)) # dimensionality of the multivariate PCA space, which is the
# number of eigenvalues > 0; usually equal to the number of descriptors.
eq.radius <- sqrt(d/p) # Radius of equilibrium circle

# Draw circle
draw.circle(0,0, eq.radius, nv=100, border="black", lty=3)

# Broken stick analysis
PCA.sig <- PCAsignificance(RL.h.pca,axes=6)
# Kaiser Guttman criterion and broken stick model
bp <- barplot(PCA.sig[2,], xlab = "PCA axes", ylab = "Relative Eigenvalues")
abline(h=mean(PCA.sig[2,]), lty=2, lwd=2, col="blue")
b_stick <- PCA.sig[4,]
points(bp, b_stick, pch=19, col="red")
lines(bp, b_stick, type="l", lwd=2, lty=3, col="red")

# Biplots PCA scaling 2
pca.sum4 <- summary(RL.h.pca, scaling = 2)
(pca.rc.env <- envfit(pca.sum4 ~ ., data = rc.env, perm = 9999, na.rm = TRUE))
par(mfrow = c(1, 1), mar=c(4,4,1,1))
explained.var <- round(pca.sum4$cont$importance[2,],3)*100
cum.explained.var <- round(pca.sum4$cont$importance[3,],3)*100
text1 <- paste("PCA 1:", explained.var[1],
               "% of the total explained variance", sep = " ")
text2 <- paste("PCA 2:", explained.var[2],
               "% of the total explained variance", sep = " ")
plot(pca.sum4$sites, pch=19, col=col1[gr],
     xlab = text1, ylab = text2,
     asp = 1, axes = TRUE, main = "PCA scaling 2",
     xlim=c(-0.3,0.3), ylim=c(-0.4,0.6))
text(pca.sum4$sites, rownames(pca.sum4$sites),
     col=col1[gr],cex = 0.7, pos=4)
arrows(pca.sum4$species[,1]*0, pca.sum4$species[,2]*0,
       pca.sum4$species[,1], pca.sum4$species[,2],
       col= "black", length = 0)
text(pca.sum4$species, rownames(pca.sum4$species), col="black",
     cex = 1, pos=4)

# Plot significant variables with a user-selected colour
p.symb.1 <- which(pca.rc.env$vectors$pvals <= 0.001) # ***
p.symb.2 <- which(pca.rc.env$vectors$pvals > 0.001 & pca.rc.env$vectors$pvals <= 0.01) # **
p.symb.3 <- which(pca.rc.env$vectors$pvals > 0.01 & pca.rc.env$vectors$pvals <= 0.05) # *
p.symb.4 <- which(pca.rc.env$vectors$pvals > 0.05 & pca.rc.env$vectors$pvals <= 0.1) # .

symb <- rep("", length(pca.rc.env$vectors$pvals))
symb[p.symb.1] <- "***"
symb[p.symb.2] <- "**"
symb[p.symb.3] <- "*"
symb[p.symb.4] <- "."

pca.rc.env2 <- pca.rc.env
rownames(pca.rc.env$vectors$arrows) <- paste(rownames(pca.rc.env2$vectors$arrows), symb, sep=" ")

plot(pca.rc.env, p.max = 0.05, col = "forestgreen")
abline(h = 0, lty = 3)
abline(v = 0, lty = 3)


# Linear relations between forest cover rate of change over a 5 year period
# and Red-listed species number
fit <- lm(X2[,8] ~ rc.env$FOREST5)
summary(fit)

# Linear relations between abandoned land cover rate of change over a 5 year period
# and Red-listed species number
fit <- lm(X2[,8] ~ rc.env$AB_LAND5)
summary(fit)


# Accumulation curves
######################

par(mfrow = c(1, 2), mar=c(4,4,1,1))
Accum.1 <- accumresult(x = spe.pa, method="exact")
Accum.1
plot(Accum.1, xlim = c(0, 22), ylim=c(0,200),
     axes=FALSE, ann = FALSE)
text(11, 210, "Species Accumulation Curve", pos=1)
text(22, 20, "Cluster 1", col = col1[1], pos=2)
text(22, 10, "Cluster 2", col = col1[2], pos=2)
text(22, 0, "Total", pos=2)
mtext("Sites", side=1, line=1.4, cex=1)
mtext("Species richness", side=2, line=2, cex=1)
mgp.axis(1, mgp = c(0, 0.4, 0), cex.axis=1)
mgp.axis(2, mgp = c(0, 0.6, 0), cex.axis=1, las=2)

Accum.1 <- accumresult(x = spe.pa[gr==2,], method="exact")
Accum.1
plot(Accum.1, col = "blue", add=T)

Accum.1 <- accumresult(x = spe.pa[gr==1,], method="exact")
Accum.1
plot(Accum.1, col = "red", add=T)
box()

Accum.2 <- accumresult(x = spe.pa[, count.t], method="exact")
Accum.2
plot(Accum.2, xlim = c(0, 22), ylim=c(0,50),
     axes=FALSE, ann = FALSE)
text(11, 52, "Red-listed Species Accumulation Curve", pos=1)
text(22, 5, "Cluster 1", col = col1[1], pos=2)
text(22, 2.5, "Cluster 2", col = col1[2], pos=2)
text(22, 0, "Total", pos=2)
mtext("Sites", side=1, line=1.4, cex=1)
mtext("Number of red-listed species", side=2, line=2, cex=1)
mgp.axis(1, mgp = c(0, 0.4, 0), cex.axis=1)
mgp.axis(2, mgp = c(0, 0.6, 0), cex.axis=1, las=2)

Accum.2 <- accumresult(x = spe.pa[gr==2, count.gr2], method="exact")
Accum.2
plot(Accum.2, col = "blue", add=T)

Accum.2 <- accumresult(x = spe.pa[gr==1, count.gr1], method="exact")
Accum.2
plot(Accum.2, col = "red", add=T)
box()

#Linear relationship between SS and abandoned land cover
fit <- lm(rc.env$SS ~ rc.env$AB_LAND2014)
summary(fit)
summary1 = paste0("SS = ", round(fit$coefficients[2],5), " x AB_LAND + ",round(fit$coefficients[1],2))
summary2 = expression(paste("Adjusted  ", R^2, " =", sep = " "))

par(mfrow=c(1,1))
plot(rc.env$AB_LAND2014, rc.env$SS, pch=19,
     col= col1[gr], xlim = c(0, 5), ylim=c(0,50),
     axes=FALSE, ann = FALSE)
text(5, 50, summary1, pos=2)
text(4.7, 48, summary2, pos=2)
text(5, 48, round(summary(fit)$adj.r.squared, 2), pos=2)
mtext("SS", side=2, line=2, cex=1)
mtext("Ab_land", side=1, line=1.4, cex=1)
mgp.axis(1, mgp = c(0, 0.4, 0), cex.axis=1)
mgp.axis(2, mgp = c(0, 0.6, 0), cex.axis=1, las=2)
abline(fit, col="black")
F.lim <- (summary(rc.env$AB_LAND2014)[4] - fit$coefficients[1]) / fit$coefficients[2]
F.lim
abline(v = F.lim, lty = 3, col = "black")
box()
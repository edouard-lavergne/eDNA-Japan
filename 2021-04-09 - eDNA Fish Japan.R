###########################################
#                                         #
#       Red-list diversity analysis       #
#      following Boccard et al. 2018      #
#                                         #
#                  by                     #
#           Edouard Lavergne              #
#                 2021                    #
#                                         #
###########################################

 
# Diversity analyses in this script follow "Numerical Ecology with R"
# by Daniel Borcard, Francois Gillet and Pierre Legendre

# All analyses presented in the manuscript can be found in this script
# Four data tables (.csv files) are also included with this script in order
# to re-run the analyses. Do not forget to change the working directory where
# you saved the .csv data files using the function setwd() for the script to work
# Data files: "2020-06 - 22rivers_HT_final.csv"
#             "2020-06 - 22rivers_LT_final.csv"
#             "2019 - Species red list.csv"
#             "2020-06 - 22rivers_Mouth_Env_data.csv"
# R version 4.0.4 (2021-04-05)

# Load the required packages
library(BiodiversityR)
library(Hmisc)
library(plotrix)
library(vegan)

# Clear all objects in R
rm(list=ls())


#########################
# Load Custom Functions #
#########################

# The "red.list" function transforms the P/A community matrix into a red-list matrix
# where lines represent sites and columns represent red-list categories.
# X is a P/A community matrix,  where lines represent sites and columns represent species
# Y is a vector indicating the red-list category of each species (DD, LC, NT, LP, VU, EN
# and CR in our case)

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
setwd("XXXX")


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


# Check for rivers that have not been sampled for both HT and LT data sets
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
                  "GOLF2014", "W_AREA", "POP2015",
                  "pH", "CRAR", "DOB", "DOS",
                  "TN", "discharge")]

# Selecting environmental variables (including High tide SST and High tide Salinity)
env.ht <- env[, c("Latitude", "Longitude", "SST_HT", "Salinity_HT", "SS", "Length_km",
                  "DENSITY2015", "PF2014", "AGRI2014", "FOREST2014",
                  "AB_LAND2014", "U_AREA2014", "RL2014",
                  "GOLF2014", "W_AREA", "POP2015",
                  "pH", "CRAR", "DOB", "DOS",
                  "TN", "discharge")]

# Selecting environmental variables (including Low tide SST and Low tide Salinity)
env.lt <- env[, c("Latitude", "Longitude", "SST_LT", "Salinity_LT", "SS", "Length_km",
                  "DENSITY2015", "PF2014", "AGRI2014", "FOREST2014",
                  "AB_LAND2014", "U_AREA2014", "RL2014",
                  "GOLF2014", "W_AREA", "POP2015",
                  "pH", "CRAR", "DOB", "DOS",
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

# How many red-listed species
spe.only.rl <- spe.rl[which(spe.rl$JMOE != "DD" & spe.rl$JMOE != "LC"),]
rl.count <- 0
for (i in 1:nrow(spe.only.rl)){
  tmp <- which(colnames(spe.pa) == rownames(spe.only.rl)[i])
  rl.count <- c(rl.count, tmp)
}
rl.count <- rl.count[-1]
paste("Number of Red-listed species:", length(rl.count), sep = " ")


(s.rich.rl <- rowSums(spe.pa[,rl.count]))

# Test of normality
shapiro.test(s.rich.rl)

# Red list Range of species richness per site
range(s.rich.rl)

# Average of Red list Species richness per site
paste("Average Red-listed species richness per site:",
      round(mean(s.rich.rl), 1),
      sep = " ")


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

# Red-listed species occurrences
range(colSums(spe.pa[, rl.count.t]))
mean(colSums(spe.pa[, rl.count.t]))

# HT
range(colSums(spe.pa.ht[, rl.count.ht]))
mean(colSums(spe.pa.ht[, rl.count.ht]))

# LT
range(colSums(spe.pa.lt[, rl.count.lt]))
mean(colSums(spe.pa.lt[, rl.count.lt]))

# Environmental data
####################

# Calculating total area from land use data
tot.area <- rowSums(Y[,c("PF2014", "AGRI2014",
                              "FOREST2014", "AB_LAND2014",
                              "U_AREA2014", "RL2014", "GOLF2014")])

# Calculating the proportion of watershed area cover for each land use
Y[,c("PF2014", "AGRI2014",
          "FOREST2014", "AB_LAND2014",
          "U_AREA2014", "RL2014",
          "GOLF2014")] <- Y[,c("PF2014", "AGRI2014",
                                          "FOREST2014", "AB_LAND2014",
                                          "U_AREA2014", "RL2014",
                                          "GOLF2014")]*100/tot.area

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

# Linear relationships between Forest and species Richness
fit <- lm(s.rich ~ Y$FOREST2014)
summary(fit)

# Linear relationships between Agri and species Richness
fit <- lm(s.rich ~ Y$AGRI2014)
summary(fit)

# Linear relations between Forest and Red-listed species number
fit <- lm((X2[,8]) ~ Y$FOREST2014)
summary(fit)
summary1 = paste0("Red listed species = ", round(fit$coefficients[2],5), " x Forest % + ",round(fit$coefficients[1],2))
summary2 = expression(paste("Adjusted  ", R^2, " =", sep = " "))

par(mfrow=c(1,1))
plot(Y$FOREST2014, X2[,8], pch=19,
     col= "black", xlim = c(40,95), ylim=c(0,12),
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

# Linear relations between Agri and Red-listed species number
fit <- lm((X2[,8]) ~ Y$AGRI2014)
summary(fit)
summary1 = paste0("Red listed species = ", round(fit$coefficients[2],5), " x Agri % + ",round(fit$coefficients[1],2))
summary2 = expression(paste("Adjusted  ", R^2, " =", sep = " "))

par(mfrow=c(1,1))
plot(Y$AGRI2014, X2[,8], pch=19,
     col= "black", xlim = c(0,13), ylim=c(0,12),
     axes=FALSE, ann = FALSE)
text(12, 11, summary1, pos=2)
text(11, 10, summary2, pos=2)
text(12, 10, round(summary(fit)$adj.r.squared, 2), pos=2)
mtext("Agri (%)", side=1, line=1.4, cex=1)
mtext("Number of red listed species", side=2, line=2, cex=1)
mgp.axis(1, mgp = c(0, 0.4, 0), cex.axis=1)
mgp.axis(2, mgp = c(0, 0.6, 0), cex.axis=1, las=2)
abline(fit, col="black")
F.lim <- (summary(rowSums(X[,c(2:6)]))[4] - fit$coefficients[1]) / fit$coefficients[2]
F.lim
abline(v = F.lim, lty = 3, col = "black")
text(Y$AGRI2014, X2[,8], rownames(X2),cex=0.8, pos=1)
box()

# Linear relations between urban area and Red-listed species number
fit <- lm(X2[,8] ~ Y$U_AREA2014)
summary(fit)

# Linear relations between abandoned land and Red-listed species number
fit <- lm(X2[,8] ~ Y$AB_LAND2014)
summary(fit)

# Linear relations between pop and Red-listed species number
tmp <- Y$POP2015/1000
fit <- lm(X2[,8] ~ tmp)
summary(fit)

# Linear relations between Latitude and Red-listed species number
fit <- lm(X2[,8] ~ Y$Latitude)
summary(fit)

# Linear relations between longitude and Red-listed species number
fit <- lm(X2[,8] ~ Y$Longitude)
summary(fit)

# RDA analysis
##############

# Run a partial RDA with land use and population data and latitude as covariate
(RL.rda <- rda(RL.hel ~ Length_km + DENSITY2015 + PF2014 +
                   AGRI2014 + FOREST2014 + AB_LAND2014 +
                   U_AREA2014 + RL2014 + GOLF2014 + W_AREA +
                   POP2015 + CRAR + SST_Av + SS + Salinity_Av +
                   pH + DOB + DOS + TN + discharge +
                 Condition(Latitude),
                 data = Y))

# Results of the RDA
(RL.rda.sum <- summary(RL.rda)) # Scaling 2 (default)

# Global test of the RDA result
anova(RL.rda, permutations = how(nperm = 9999))

# Tests of all canonical axes
anova(RL.rda, by = "axis", permutations = how(nperm = 9999))

# Unadjusted R^2 retrieved from the rda object
(R2 <- RsquareAdj(RL.rda)$r.squared)

# Adjusted R^2 retrieved from the rda object
(R2adj <- RsquareAdj(RL.rda)$adj.r.squared)

# Apply Kaiser-Guttman criterion to residual axes
# Might indicate that additional information are located
# in the residual
RL.rda$CA$eig[RL.rda$CA$eig > mean(RL.rda$CA$eig)]

# Variance inflation factors (VIF) (should be <10)
vif.cca(RL.rda)

# Forward selection using vegan's ordistep()
# This function allows the use of factors.
mod0 <- rda(RL.hel ~ Condition(Latitude), data = Y[,-2])
step.forward <- ordistep(mod0,
                         scope = formula(RL.rda),
                         Pout = 0.05,
                         direction = "forward",
                         permutations = how(nperm = 9999))
RsquareAdj(step.forward)

# Final RDA
(RL.rda.f <- rda(RL.hel ~ AGRI2014 + SS + FOREST2014 + POP2015 + DOS + Condition(Latitude),
               data = Y))

# Results of the RDA
(RL.rda.f.sum <- summary(RL.rda.f)) # Scaling 2 (default)

# Global test of the RDA result
anova(RL.rda.f, permutations = how(nperm = 9999))

# Tests of all canonical axes
anova(RL.rda.f, by = "axis", permutations = how(nperm = 9999))

# Unadjusted R^2 retrieved from the rda object
(R2 <- RsquareAdj(RL.rda.f)$r.squared)

# Adjusted R^2 retrieved from the rda object
(R2adj <- RsquareAdj(RL.rda.f)$adj.r.squared)

# Apply Kaiser-Guttman criterion to residual axes
# Might indicate that additional information are located
# in the residual
RL.rda.f$CA$eig[RL.rda.f$CA$eig > mean(RL.rda.f$CA$eig)]


# Variance inflation factors (VIF) (should be <10)
vif.cca(RL.rda.f)

RL.h.pca <- RL.rda.f
# Biplots PCA scaling 1
pca.sum <- summary(RL.h.pca, scaling = 1)
(RL.h.pca.env <- envfit(pca.sum ~ AGRI2014 + SS + FOREST2014 +
                          POP2015 + DOS,
                        data = Y, perm = 9999, scaling = 1))
RL.h.pca.env2 <- RL.h.pca.env # back up
par(mfrow = c(1, 1), mar=c(4,4,1,1))
(explained.var <- round(pca.sum$concont$importance[2,]*R2adj,3)*100)
(cum.explained.var <- round(pca.sum$concont$importance[3,]*R2adj,3)*100)
text1 <- paste("RDA 1:", explained.var[1],
               "% of the total explained variance", sep = " ")
text2 <- paste("RDA 2:", explained.var[2],
               "% of the total explained variance", sep = " ")
plot(pca.sum$sites, pch=19, cex = X2[,8]/5, col="black",
     xlab = text1, ylab = text2,
     asp = 1, axes = TRUE, main = "Constrained RDA scaling 1",
     xlim=c(-0.6,0.6), ylim=c(-0.6,0.6))
text(pca.sum$sites, rownames(pca.sum$sites),
     col="black",cex = 0.7, pos=4)
arrows(pca.sum$species[,1]*0, pca.sum$species[,2]*0,
       pca.sum$species[,1], pca.sum$species[,2],
       col= "firebrick2", length = 0)
text(pca.sum$species, rownames(pca.sum$species), col="firebrick2",
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
plot(RL.h.pca.env, p.max = 0.5, col = "dodgerblue2")
abline(h = 0, lty = 3)
abline(v = 0, lty = 3)

# Equilibrium circle: sqrt(d/p)
d <- 2 # dimension of the reduced space (usually, d = 2)
p <-  length(which(pca.sum$cont$importance[1,]>0)) # dimensionality of the multivariate PCA
                                                    # space, which is the number of eigenvalues > 0;
                                                    # usually equal to the number of descriptors.
eq.radius <- sqrt(d/p) # Radius of equilibrium circle

# Draw circle
draw.circle(0,0, eq.radius, nv=10000, border="firebrick2", lty=3)

# Accumulation curves
######################

par(mfrow = c(1, 2))
Accum.1 <- accumresult(x = spe.pa, method="exact")
Accum.1
plot(Accum.1, xlim = c(0, 22), ylim=c(0,200),
     axes=FALSE, ann = FALSE)
text(11, 210, "Species Accumulation Curve", pos=1)
mtext("Sites", side=1, line=1.4, cex=1)
mtext("Species richness", side=2, line=2, cex=1)
mgp.axis(1, mgp = c(0, 0.4, 0), cex.axis=1)
mgp.axis(2, mgp = c(0, 0.6, 0), cex.axis=1, las=2)
box()

Accum.2 <- accumresult(x = spe.pa[, rl.count], method="exact")
Accum.2
plot(Accum.2, xlim = c(0, 22), ylim=c(0,200),
     axes=FALSE, ann = FALSE)
text(11, 210, "Red-listed Species Accumulation Curve", pos=1)
mtext("Sites", side=1, line=1.4, cex=1)
mtext("Number of red-listed species", side=2, line=2, cex=1)
mgp.axis(1, mgp = c(0, 0.4, 0), cex.axis=1)
mgp.axis(2, mgp = c(0, 0.6, 0), cex.axis=1, las=2)
box()


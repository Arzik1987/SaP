
# check whether the needed packages are installed. Install if something is missing

list.of.packages <- c("lubridate", "anytime")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)


# load required packages
# library(utils)
# library(magrittr)
library(anytime)
library(lubridate)
library(bpriv)



# load the file specifying the datasets selected for the experiments. 
# One can inspect this file to see which data was used.
# One can also modify this to add another experiments
dir <- getwd()
data.info <- read.csv(paste0(dir, '/data_preparation.csv'), stringsAsFactors = FALSE)


# The function which changes discretization rate.
# freq - required sampling rate in seconds
# the function assumes the input is in Watts measured each second without missing measurements,
# and outputs the measurements in kWh. It also assumes that the measurements start right after midnight
# (signifficant for a feature creation process)

resample <- function(lp, freq) {
  llp <- length(lp)
  
  # we need features associated with datasets to compute one version of mutual information - MI_v
  features <- c(rep(1, 6.5*3600), rep(2, 3*3600), rep(3, 6*3600), rep(4, 7*3600), rep(1, 1.5*3600))
  if(llp%%length(features) > 0) stop("missing values in load profile")
  features <- rep(features, llp%/%length(features))
  
  # if load profile is not exactly divisible to freq, something might be wrong.
  ind <- (llp - (llp %% freq))
  if(ind < llp) warning("the length of load profile is not exactly divisible to frequency")
  lp <- lp[1:ind]
  features <- features[1:ind]
  
  lp <- lp/3600000 # convert to kWh
  
  # now we aggregate the subsequent values of load to obtained the required sampling rate
  # we also mak the feature vector sparser
  ind <- ((1:length(lp)) - 1) %/% freq
  lp <- as.numeric(tapply(lp, ind, sum))
  features <- as.numeric(tapply(features, ind, max))
  return (list(lpo = lp, features = features, interval = freq/60))
}


# the function inputs missing values

input.miss <- function(d){
  # If a very first value is missing - set it to the first not-missing value
  first <- min(which(!is.na(d)))
  d[1:(first - 1)] <- d[first]
  # for each missing measurement define it equal to the previous measurement
  for(i in 2:length(d)){
    if(is.na(d[i])) d[i] <- d[i - 1]
  }
  d
}


#### importing ECO data ####

# Import a single ECO profile
# ECO.dir - the path to the directory with ECO data
# j - rowindex in the file "data_preparation.csv, specifying the profile

import.ECO.profile <- function(j, ECO.dir){
  
  flstart <- anydate(data.info$Start_date[j]) - data.info$Warm_up_days[j]
  flend <- anydate(data.info$End_date[j])
  files <- paste0(anydate(flstart:flend), ".csv")
  
  unzip(paste0(ECO.dir, "/", data.info$Pathpart[j], "_sm_csv.zip"), 
        files = paste0(data.info$Pathpart[j], "/", files), 
        exdir = ECO.dir, overwrite = TRUE)
  cat(paste("load profile number", j, ":    ", files, "\n"))
  
  lp = c()
  for (i in 1:length(files)) {
    lp = c(lp, read.csv(paste0(ECO.dir, "/", data.info$Pathpart[j], "/", files[i]), header = FALSE, sep = ',')[, 1])
  }
  return(resample(lp, data.info$Sample_rate_seconds[j]))
}


# Import data

ECO.dir <- paste0(dir, '/data/ECO')
ECO.lps <- list()
ECOinds <- (1:nrow(data.info))[grepl("ECO", data.info$Name)]
ECO.lps <- lapply(ECOinds, import.ECO.profile, ECO.dir = ECO.dir)
names(ECO.lps) <- data.info$Name[ECOinds]


#### REDD ####

# Import a single REDD profile
# REDD.dir - the path to the directory with ECO data
# j - rowindex in the file "data_preparation.csv, specifying the profile

import.REDD.profile <- function(j, REDD.dir){
  
  lp <- read.table(paste0(REDD.dir, '/low_freq/', data.info$Pathpart[j], '/channel_1.dat'), header = FALSE)

  d = lp[anydate(lp[, 1]) >= (anydate(data.info$Start_date[j]) - data.info$Warm_up_days[j]),]
  d = d[anydate(d[, 1]) <= anydate(data.info$End_date[j]), ] 
  
  low <- as.numeric(anytime(data.info$Start_date[j]) - data.info$Warm_up_days[j]*3600*24)
  high <- as.numeric(anytime(data.info$End_date[j])) + 3600*24 - 1
  full <- as.data.frame(low:high)
  
  colnames(full) <- "dates"
  colnames(d) <- c("dates", "vals")
  d <- merge(full, d, all = TRUE)
  d <- input.miss(d[, 2])
  
  resample(d, data.info$Sample_rate_seconds[j])
}


# Import data

REDD.dir <- paste0(dir, '/data/REDD')
# IMPORTANT: first extract files from '.bz2' archive
untar(paste0(REDD.dir, '/low_freq.tar'), exdir = REDD.dir, 
      files = c('low_freq/house_1/channel_1.dat', 'low_freq/house_2/channel_1.dat'))

REDD.lps <- list()
REDDinds <- (1:nrow(data.info))[grepl("REDD", data.info$Name)]
REDD.lps <- lapply(REDDinds, import.REDD.profile, REDD.dir = REDD.dir)
names(REDD.lps) <- data.info$Name[REDDinds]


#### Smart ####

# Import a single Smart profile
# Smart.dir - the path to the directory with ECO data
# j - rowindex in the file "data_preparation.csv, specifying the profile


import.Smart.profile <- function(j, Smart.dir) {
  
  flstart <- anydate(data.info$Start_date[j]) - data.info$Warm_up_days[j]
  flend <- anydate(data.info$Start_date[j])
  files <- anydate(flstart:flend)
  lp <- data.frame(matrix(ncol = 2, nrow = 0))

  for(i in 1:length(files)){
    tmp <- paste0(year(files[i]), "-", substr(month.name, 1, 3)[month(files[i])], "-", day(files[i]), ".csv")
    lp <- rbind(lp, read.csv(paste0(Smart.dir, '/homeB-power/', tmp), header = FALSE))
  }
  
  low <- as.numeric(anytime(flstart))
  high <- as.numeric(anytime(flend)) + 3600*24 - 1
  full <- as.data.frame(low:high)
  
  colnames(full) <- "dates"
  colnames(lp) <- c("dates", "vals")
  lp <- merge(full, lp, all = TRUE)
  lp <- input.miss(lp[, 2])
  
  resample(lp, data.info$Sample_rate_seconds[j])
}


# Import data

Smart.dir <- paste0(dir, '/data/Smart')
untar(paste0(Smart.dir, '/homeB-power.tar.gz'), exdir = Smart.dir)

Smart.lps <- list()
Smartinds <- (1:nrow(data.info))[grepl("Smart", data.info$Name)]
Smart.lps <- lapply(Smartinds, import.Smart.profile, Smart.dir = Smart.dir)
names(Smart.lps) <- data.info$Name[Smartinds]

#### CER ####

# Import a single CER profile
# CER.dir - the path to the directory with ECO data
# j - rowindex in the file "data_preparation.csv, specifying the profile

import.CER.profile <- function(j, CER.dir) {

  lp = read.csv(paste0(CER.dir, "/", data.info$Pathpart[j], ".txt"), header = FALSE)
  if(data.info$Time[j] == "Winter"){
    lp = lp[6721:(6721 + 671), 1]
  } else {
    lp = lp[1:672, 1]
  }
  llp <- length(lp)
  features <- c(rep(1, 6.5*2), rep(2, 3*2), rep(3, 6*2), rep(4, 7*2), rep(1, 1.5*2))
  features <- rep(features, llp%/%length(features))
  
  return(list(lpo = lp, features = features, interval = 30))
}


# Import data

CER.dir <- paste0(dir, '/data/CER')
CERinds <- (1:nrow(data.info))[grepl("CER", data.info$Name)]
CER.lps <- lapply(CERinds, import.CER.profile, CER.dir = CER.dir)
names(CER.lps) <- data.info$Name[CERinds]


#### Combine all load profiles and save ####

all.lps <- c(ECO.lps, REDD.lps, Smart.lps, CER.lps)
save(all.lps, file = paste0(dir, "/all_lps.RData"))

for(i in 1:length(all.lps)){
  plot(all.lps[[i]][[1]], type = "l", col = "red", main = names(all.lps)[[i]],
       xlab = paste0("Interval = ", all.lps[[i]]$interval), ylab = "kWh")
}


# check whether the needed packages are installed. Install if something is missing

list.of.packages <- c("utils", "magrittr", "anytime")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)


# load required packages
# library(utils)
# library(magrittr)
library(anytime)
library(bpriv)


# load the file specifying the datasets selected for the experiments. 
# One can inspect this file to see which data was used.
# One can also modify this to add another experiments
dir <- getwd()
data.info <- read.csv(paste0(dir, '/data_preparation.csv'))


# The function which changes discretization rate.
# freq - required sampling rate in seconds
# the function assumes the input is in Watts measured each second without missing measurements,
# and outputs the measurements in kWt

resample <- function(lp, freq) {
  llp <- length(lp)
  features <- c(rep(1, 6.5*3600), rep(2, 3*3600), rep(3, 6*3600), rep(4, 7*3600), rep(1, 1.5*3600))
  if(llp%%length(features) > 0) stop("missing values in load profile")
  features <- rep(features, llp%/%length(features))
  
  lp <- lp[1:(llp - (llp %% freq))]
  features <- features[1:(llp - (llp %% freq))]
  
  lp <- lp/3600000
  ind <- ((1:length(lp)) - 1) %/% freq
  lp <- as.numeric(tapply(lp, ind, sum))
  features <- as.numeric(tapply(features, ind, max))
  return (list(lpo = lp, features = features, interval = freq/60))
}


# the function inputs missing values

input.miss <- function(d){
  first <- min(which(!is.na(d)))
  d[1:(first - 1)] <- d[first]
  for(i in 1:length(d)){
    if(is.na(d[i])) d[i] <- d[i - 1]
  }
  d
}


    #### ECO ####


# Import a single ECO profile

import.ECO.profile <- function(basepath, filepath) {
  return (read.csv(paste0(basepath, filepath), header = FALSE, sep = ',')[, 1])
}


# Import multiple consecutive ECO profiles; used for low sample rate

import.ECO.profiles <- function(basepath, filepaths, freq) {
  lp = c()
  for (i in filepaths) {
    lp = c(lp, import.ECO.profile(basepath, i))
  }
  return(resample(lp, freq))
}


# Import all ECO profiles

import.all.ECO.profiles <- function(ECO.dir, data.info, prefix, range) {
  ECO.file.names = data.info$File.name
  lps = list()
  for (i in 1:length(range)) {
    file.names = strsplit(ECO.file.names[range[i]], split = " ")[[1]]
    unzip(paste0(ECO.dir, "/", prefix, "_sm_csv.zip"), files = paste0(prefix, "/", file.names), 
          exdir = ECO.dir, overwrite = TRUE)
    cat(paste0(range[i], "/24 ", file.names, "\n"))
    if (length(file.names) == 1) {
      lps[[i]] = import.ECO.profile(basepath = paste0(ECO.dir, "/", prefix, "/"), filepath = file.names) %>% 
        resample(.,data.info$`Sample.rate.(s)`[range[i]])
    } else {
      lps[[i]] = import.ECO.profiles(basepath = paste0(ECO.dir, "/", prefix, "/"), 
                                     filepath = file.names, freq = data.info$`Sample.rate.(s)`[range[i]])
    }
  }
  return(lps)
}


# House 1: 4 residents
ECO.dir <- paste0(dir, '/data/ECO')
ECO.lps <- import.all.ECO.profiles(ECO.dir = ECO.dir, data.info = data.info, prefix = "01", range = 1:12)

# House 2: 2 residents
ECO.lps <- append(ECO.lps, import.all.ECO.profiles(ECO.dir = ECO.dir, data.info = data.info, 
                                                   prefix = "02", range = 13:24))


    #### REDD ####


REDD.dir <- paste0(dir, '/data/REDD')

# !!!! did not work with bz2...
untar(paste0(REDD.dir, '/low_freq.tar'), exdir = REDD.dir, 
      files = c('low_freq/house_1/channel_1.dat', 'low_freq/house_2/channel_1.dat'))
REDD.data <- read.table(paste0(REDD.dir, '/low_freq/house_1/channel_1.dat'), header = FALSE)
REDD.data[, 3] <- anydate(REDD.data[, 1])

import.REDD.profile <- function(date, name) {
  d <- REDD.data[REDD.data[, 3] == anydate(date), 1:2]
  low <- as.numeric(anytime(date))
  high <- as.numeric(anytime(date)) + 3600*24 - 1
  full <- as.data.frame(low:high)
  
  colnames(full) <- "dates"
  colnames(d) <- c("dates", "vals")
  d <- merge(full, d, all = TRUE)
  d <- input.miss(d[, 2])
  
  resample(d, data.info$`Sample.rate.(s)`[data.info$Name == name])
}

import.REDD.profiles <- function(start.date, end.date, name) {
  d = REDD.data[anydate(REDD.data[, 3]) >= anydate(start.date),]
  d = d[anydate(d[, 3]) <= anydate(end.date), 1:2] 
  
  low <- as.numeric(anytime(start.date))
  high <- as.numeric(anytime(end.date)) + 3600*24 - 1
  full <- as.data.frame(low:high)
  
  colnames(full) <- "dates"
  colnames(d) <- c("dates", "vals")
  d <- merge(full, d, all = TRUE)
  d <- input.miss(d[, 2])
  
  resample(d, data.info$`Sample.rate.(s)`[data.info$Name == name])
}

REDD.lps <- mapply(import.REDD.profile,
                   c('2011-04-19', '2011-04-23', '2011-04-24', '2011-04-20'),
                   c('REDD1', 'REDD2', 'REDD3', 'REDD4'), SIMPLIFY = FALSE)

REDD.data <- read.table(paste(REDD.dir,'/low_freq/house_2/channel_1.dat', sep = ""), header=FALSE)
REDD.data[, 3] <- anydate(REDD.data[, 1])
REDD.lps[[5]] <- import.REDD.profiles('2011-04-19', '2011-04-25', 'REDD5')



### new

# d <- REDD.data[REDD.data[, 3] == anydate('2011-04-24'), 1:2]
# low <- as.numeric(anytime('2011-04-24'))
# high <- as.numeric(anytime('2011-04-24')) + 3600*24 - 1
# full <- as.data.frame(low:high)
# 
# colnames(full) <- "dates"
# colnames(d) <- c("dates", "vals")
# d <- merge(full, d, all = TRUE)

### end new

    #### Smart ####


Smart.dir <- paste0(dir, '/data/Smart')
untar(paste(Smart.dir, '/homeB-power.tar.gz', sep = ""), exdir = Smart.dir)

import.Smart.profile <- function(d, name) {
  
  low <- as.numeric(anytime(anydate(d[1, 1])))
  high <- as.numeric(anytime(anydate(d[1, 1]))) + 3600*24 - 1
  full <- as.data.frame(low:high)
  
  colnames(full) <- "dates"
  colnames(d) <- c("dates", "vals")
  d <- merge(full, d, all = TRUE)
  d <- input.miss(d[, 2])
  
  resample(d, data.info$`Sample.rate.(s)`[data.info$Name == name])
}

Smart.lps = list()
Smart.data <- read.csv(paste0(Smart.dir, '/homeB-power/2012-Apr-15.csv'), header = FALSE)
Smart.lps[[1]] <- import.Smart.profile(Smart.data, 'SmartB1')
Smart.data <- read.csv(paste0(Smart.dir, '/homeB-power/2012-Apr-16.csv'), header = FALSE)
Smart.lps[[2]] <- import.Smart.profile(Smart.data, 'SmartB2')


    #### CER ####


CER.dir <- paste0(dir, '/data/CER')
import.CER.profiles <- function(ind) {
  filename <- data.info$File.name[ind]
  lp = read.csv(paste0(CER.dir, "/", filename), header = FALSE)
  if (data.info$`Location/Time`[ind] %>% grepl("Winter", .)) {
    lp = lp[6721:(6721 + 671), 1]
  } else {
    lp = lp[1:672, 1]
  }
  llp <- length(lp)
  features <- c(rep(1, 6.5*2), rep(2, 3*2), rep(3, 6*2), rep(4, 7*2), rep(1, 1.5*2))
  features <- rep(features, llp%/%length(features))
  
  return(list(lpo = lp, features = features, interval = 30))
}

CER.lps <- lapply(32:37, import.CER.profiles)

    #### Combine and save ####

all.lps <- ECO.lps %>% append(., REDD.lps) %>% append(., Smart.lps) %>% append(., CER.lps)
names(all.lps) <- data.info$Name
save(all.lps, file = paste0(dir, "/all_lps.RData"))



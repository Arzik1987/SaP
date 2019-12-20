# list.of.packages <- c("batchtools")
# new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
# if(length(new.packages)) install.packages(new.packages)

library(batchtools)
library(bpriv)

setwd("C:\\Projects\\4 Storage and Privacy\\Experiments_gitHub")

#### if you did something wrong and want to redo everything from scratch

# reg <- loadRegistry(file.dir = paste0(getwd(), "/registry"), work.dir = getwd(), writeable = TRUE)
# removeRegistry(reg = reg)

#### for the first start of this code you need to make a registry, where the results will reside

reg = makeExperimentRegistry(file.dir = paste0(getwd(), "/registry"), packages = "bpriv", seed = 1)

#### for subsequent uses of the code, when the registry exists

# reg <- loadRegistry(file.dir = paste0(getwd(), "/registry"), work.dir = getwd(), writeable = TRUE)

### if you want to remove experimental results, but keep the registry

# reg <- loadRegistry(file.dir = paste0(getwd(), "/registry"), work.dir = getwd(), writeable = TRUE)
# clearRegistry(reg = getDefaultRegistry())

#### parallelize

reg$cluster.functions = makeClusterFunctionsSocket()
saveRegistry()

#### BBLHs ("problems" in bacthtools terminology)

load(paste0(getwd(), "/all_lps.RData"))

modifylp <- function(data, job, algo.name, max.be, max.bp, be, mode, ...){
  possible.args <- list(max.be = max.be, max.bp = max.bp, be = be, mode = mode)
  arg.names <- names(as.list(args(algo.name)))
  args.to.use <- possible.args[names(possible.args) %in% arg.names]

  #### check if better load profiles exist
  data$lpo[data$lpo < 0] <- 0
  
  if("voltage" %in% arg.names){
    d <- data[c("lpo", "interval", "voltage")]
  } else {
    d <- data[c("lpo", "interval")]
  }
  
    res <- do.call(algo.name, append(d, args.to.use))
  
  start <- data$warmup + 1
  end <- length(data$lpo)

  list(lpo = data$lpo[start:end], lpm = res[start:end], 
       features = data$features[start:end], interval = data$interval)
}

for(i in 1:length(all.lps)){
  addProblem(name = names(all.lps)[i], data = all.lps[[i]], fun = modifylp)
}

#### privacy measures ("algorithms" in batchtools terminolohy)

# cluster similarity

meas.cs.orig <- function(data, job, instance, ...) {
  res <- priv.cs(lpo = instance$lpo, lpm = instance$lpm, ...)
  res
}

meas.cs.diff <- function(data, job, instance, ...) {
  res <- priv.cs(lpo = diff(instance$lpo), lpm = diff(instance$lpm), ...)
  res
}

# conditional entropy

meas.ce.orig <- function(data, job, instance, ...) {
  res <- priv.ce(lpo = instance$lpo, lpm = instance$lpm, ...)
  -res
}

meas.ce.diff <- function(data, job, instance, ...) {
  res <- priv.ce(lpo = diff(instance$lpo), lpm = diff(instance$lpm), ...)
  -res
}

# entropy ratio

meas.er.hist.diff <- function(data, job, instance, ...) {
  res <- priv.er.hist(lpo = diff(instance$lpo), lpm = diff(instance$lpm), ...)
  res
}

# feature mass

meas.fm.diff <- function(data, job, instance, regime, ...) {
  thr <- ifelse(regime == "rfm", 0, instance$interval*50/60000) # threshold = 50 Watt for "fm" and "ed" regimes
  res <- priv.fm(lpo = diff(instance$lpo), lpm = diff(instance$lpm), thr = thr,
                 regime = regime, ...) 
  res
}

# K-divergence

meas.kd.orig <- function(data, job, instance, ...) {
  res <- priv.kd(lpo = instance$lpo, lpm = instance$lpm, ...)
  -res
}

# differenced is commented out since it is only applied as conditional 9and we do not program it)
#
# meas.kd.diff <- function(data, job, instance, ...) {
#   res <- priv.kd(lpo = diff(instance$lpo), lpm = diff(instance$lpm), ...)
#   -res
# }

# KL divergence

meas.kl.hist.orig <- function(data, job, instance, ...) {
  res <- priv.kl.hist(lpo = instance$lpo, lpm = instance$lpm, ...)
  -res
}

meas.kl.hist.diff <- function(data, job, instance, ...) {
  res <- priv.kl.hist(lpo = diff(instance$lpo), lpm = diff(instance$lpm), ...)
  -res
}

# Load variance

meas.lv.orig <- function(data, job, instance, ...) {
  res <- priv.lv(lpo = instance$lpo, lpm = instance$lpm, ...)
  res
}

# Mutual information

meas.mi.hist.orig <- function(data, job, instance, ...) {
  res <- priv.mi.hist(lpo = instance$lpo, lpm = instance$lpm, features = instance$features, ...)
  res
}

meas.mi.hist.diff <- function(data, job, instance, ...) {
  res <- priv.mi.hist(lpo = diff(instance$lpo), lpm = diff(instance$lpm), features = instance$features, ...)
  res
}

# Reconstruction

meas.rc.orig <- function(data, job, instance, ...) {
  res <- priv.rc(lpo = instance$lpo, lpm = instance$lpm, ...)
  res
}

# coefficient of determination

meas.dc.orig <- function(data, job, instance, ...) {
  res <- priv.dc(lpo = instance$lpo, lpm = instance$lpm, ...)
  res
}

meas.dc.diff <- function(data, job, instance, ...) {
  res <- priv.dc(lpo = diff(instance$lpo), lpm = diff(instance$lpm), ...)
  res
}

# Total variation distance

meas.tvd.orig <- function(data, job, instance, ...) {
  res <- priv.tvd(lpo = instance$lpo, lpm = instance$lpm, ...)
  -res
}


####

addAlgorithm(name = "meas.cs.orig", fun = meas.cs.orig)
addAlgorithm(name = "meas.cs.diff", fun = meas.cs.diff)

addAlgorithm(name = "meas.ce.orig", fun = meas.ce.orig)
addAlgorithm(name = "meas.ce.diff", fun = meas.ce.diff)

addAlgorithm(name = "meas.er.hist.diff", fun = meas.er.hist.diff)

addAlgorithm(name = "meas.fm.diff", fun = meas.fm.diff)

addAlgorithm(name = "meas.kd.orig", fun = meas.kd.orig)

addAlgorithm(name = "meas.kl.hist.orig", fun = meas.kl.hist.orig)
addAlgorithm(name = "meas.kl.hist.diff", fun = meas.kl.hist.diff)

addAlgorithm(name = "meas.lv.orig", fun = meas.lv.orig)

addAlgorithm(name = "meas.mi.hist.orig", fun = meas.mi.hist.orig)
addAlgorithm(name = "meas.mi.hist.diff", fun = meas.mi.hist.diff)

addAlgorithm(name = "meas.rc.orig", fun = meas.rc.orig)

addAlgorithm(name = "meas.dc.orig", fun = meas.dc.orig)
addAlgorithm(name = "meas.dc.diff", fun = meas.dc.diff)

addAlgorithm(name = "meas.tvd.orig", fun = meas.tvd.orig)



#### check what is there

reg$problems
reg$algorithms

#### parameters (names of BBLHs to use and storage characteristics - inputs of BBLHs)

ades = list(
  meas.cs.orig = data.table(),
  meas.cs.diff = data.table(),
  meas.ce.orig = data.table(),
  meas.ce.diff = data.table(),
  meas.er.hist.diff = data.table(regime = c("zeros", "nozeros")),
  meas.fm.diff = data.table(regime = c("fm", "rfm", "ed")),
  meas.kd.orig = data.table(),
  meas.kl.hist.orig = data.table(),
  meas.kl.hist.diff = data.table(),
  meas.lv.orig = data.table(),
  meas.mi.hist.orig = data.table(regime = c("iid", "ms", "mns", "bin")),
  meas.mi.hist.diff = data.table(regime = c("iid", "ms", "bin")),
  meas.rc.orig = data.table(regime = c("v", "w")),
  meas.dc.orig = data.table(regime = c("reg", "lpo")),
  meas.dc.diff = data.table(regime = c("reg")),
  meas.tvd.orig = data.table()
)


tmp1 <- rbind(
  CJ(ind = 1:4, algo.name = "alg.be", mode = c("first", "second")),
  CJ(ind = 1:4, algo.name = "alg.nill", mode = "na"),
  CJ(ind = 1:4, algo.name = "alg.stepping", mode = c("ls1", "ls2", "lc", "rc"))
)

tmp2 <- data.table(ind = 1:4,
                   max.be = c(1.2, 4.0, 2.4, 6.5),
                   max.bp = c(0.27, 2.0, 4.0, 5.0),
                   be = c(0.5, 0.5, 0.5, 0.5))

pars <- merge(tmp1, tmp2, by = "ind")[, -1]

pdes = list()
for(i in 1:length(all.lps)){
  pdes <- append(pdes, list(pars))
}
names(pdes) <- names(all.lps)


addExperiments(pdes, ades, repls = 1)


#### check what is where

summarizeExperiments()
unwrap(getJobPars())

#### run experiments

submitJobs()
waitForJobs()




# check whether the needed packages are installed. Install if something is missing

list.of.packages <- c("batchtools", "devtools", "purrr", "infotheo", "cluster", "entropy", "FNN", "wavelets")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

if(!("bpriv" %in% installed.packages()[,"Package"])){
  install.packages("bpriv_0.0.1.tar.gz", repos = NULL, type = "source")
}

library(batchtools)
library(data.table)
library(bpriv)


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



######################
#### BBLHs ("problems" in bacthtools terminology)
######################

load(paste0(getwd(), "/all_lps.RData")) # load pre-processed load profiles

modifylp <- function(data, job, algo.name, max.be, max.bp, be, mode, ...){
  possible.args <- list(max.be = max.be, max.bp = max.bp, be = be, mode = mode)
  arg.names <- names(as.list(args(algo.name)))
  args.to.use <- possible.args[names(possible.args) %in% arg.names]

  #### if there are some errors...
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


######################
#### privacy measures ("algorithms" in batchtools terminolohy)
######################


# cluster similarity

cs.orig <- function(data, job, instance, ...) {
  res <- priv.cs(lpo = instance$lpo, lpm = instance$lpm, ...)
  res
}

cs.diff <- function(data, job, instance, ...) {
  res <- priv.cs(lpo = diff(instance$lpo), lpm = diff(instance$lpm), ...)
  res
}

# conditional entropy

ce.orig <- function(data, job, instance, ...) {
  res <- priv.ce(lpo = instance$lpo, lpm = instance$lpm, ...)
  -res
}

ce.diff <- function(data, job, instance, ...) {
  res <- priv.ce(lpo = diff(instance$lpo), lpm = diff(instance$lpm), ...)
  -res
}

# entropy ratio

er.diff <- function(data, job, instance, ...) {
  res <- priv.er.hist(lpo = diff(instance$lpo), lpm = diff(instance$lpm), ...)
  res
}

# feature mass

fm.diff <- function(data, job, instance, regime, ...) {
  thr <- ifelse(regime == "rfm", 0, instance$interval*50/60000) # threshold = 50 Watt for "fm" and "ed" regimes
  res <- priv.fm(lpo = diff(instance$lpo), lpm = diff(instance$lpm), thr = thr,
                 regime = regime, ...) 
  res
}

# K-divergence

kd.orig <- function(data, job, instance, ...) {
  res <- priv.kd(lpo = instance$lpo, lpm = instance$lpm, ...)
  -res
}

# differenced is commented out since it is only applied as conditional K-divergence 
# (and we do not use it)
#
# kd.diff <- function(data, job, instance, ...) {
#   res <- priv.kd(lpo = diff(instance$lpo), lpm = diff(instance$lpm), ...)
#   -res
# }

# KL divergence

kl.orig <- function(data, job, instance, ...) {
  res <- priv.kl.hist(lpo = instance$lpo, lpm = instance$lpm, ...)
  -res
}

kl.diff <- function(data, job, instance, ...) {
  res <- priv.kl.hist(lpo = diff(instance$lpo), lpm = diff(instance$lpm), ...)
  -res
}

# Load variance

lv.orig <- function(data, job, instance, ...) {
  res <- priv.lv(lpo = instance$lpo, lpm = instance$lpm, ...)
  res
}

# Mutual information

mi.orig <- function(data, job, instance, ...) {
  res <- priv.mi.hist(lpo = instance$lpo, lpm = instance$lpm, features = instance$features, ...)
  res
}

mi.diff <- function(data, job, instance, ...) {
  res <- priv.mi.hist(lpo = diff(instance$lpo), lpm = diff(instance$lpm), features = instance$features, ...)
  res
}

# Reconstruction

rc.orig <- function(data, job, instance, ...) {
  res <- priv.rc(lpo = instance$lpo, lpm = instance$lpm, ...)
  res
}

# coefficient of determination

dc.orig <- function(data, job, instance, ...) {
  res <- priv.dc(lpo = instance$lpo, lpm = instance$lpm, ...)
  res
}

dc.diff <- function(data, job, instance, ...) {
  res <- priv.dc(lpo = diff(instance$lpo), lpm = diff(instance$lpm), ...)
  res
}

# Total variation distance

tvd.orig <- function(data, job, instance, ...) {
  res <- priv.tvd(lpo = instance$lpo, lpm = instance$lpm, ...)
  -res
}


####

addAlgorithm(name = "cs.orig", fun = cs.orig)
addAlgorithm(name = "cs.diff", fun = cs.diff)

addAlgorithm(name = "ce.orig", fun = ce.orig)
addAlgorithm(name = "ce.diff", fun = ce.diff)

addAlgorithm(name = "er.diff", fun = er.diff)

addAlgorithm(name = "fm.diff", fun = fm.diff)

addAlgorithm(name = "kd.orig", fun = kd.orig)

addAlgorithm(name = "kl.orig", fun = kl.orig)
addAlgorithm(name = "kl.diff", fun = kl.diff)

addAlgorithm(name = "lv.orig", fun = lv.orig)

addAlgorithm(name = "mi.orig", fun = mi.orig)
addAlgorithm(name = "mi.diff", fun = mi.diff)

addAlgorithm(name = "rc.orig", fun = rc.orig)

addAlgorithm(name = "dc.orig", fun = dc.orig)
addAlgorithm(name = "dc.diff", fun = dc.diff)

addAlgorithm(name = "tvd.orig", fun = tvd.orig)



#### check what is there

reg$problems
reg$algorithms


######################
#### parameters (variants of measures, variants of BBLHs to use and storage characteristics - inputs of BBLHs)
######################

ades = list(
  cs.orig = data.table(),
  cs.diff = data.table(),
  ce.orig = data.table(),
  ce.diff = data.table(),
  er.diff = data.table(regime = c("zeros", "nozeros")),
  fm.diff = data.table(regime = c("fm", "rfm", "ed")),
  kd.orig = data.table(),
  kl.orig = data.table(),
  kl.diff = data.table(),
  lv.orig = data.table(),
  mi.orig = data.table(regime = c("iid", "ms", "mns")),
  mi.diff = data.table(regime = c("iid", "ms", "bin")),
  rc.orig = data.table(regime = c("v", "w")),
  dc.orig = data.table(regime = c("reg", "lpo")),
  dc.diff = data.table(regime = c("reg")),
  tvd.orig = data.table()
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



######################
#### run experiments
######################

submitJobs()
waitForJobs()



# list.of.packages <- c("batchtools")
# new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
# if(length(new.packages)) install.packages(new.packages)

library(bpriv)
library(batchtools)

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

modifylp <- function(data, job, algo.name, max.be, min.be, max.bc, max.bd, max.bp, be, ...){
  possible.args <- list(max.be = max.be, min.be = min.be, max.bc = max.bc, max.bd = max.bd, max.bp = max.bp, be = be)
  arg.names <- names(as.list(args(algo.name)))
  args.to.use <- possible.args[names(possible.args) %in% arg.names]

  #### check if better load profiles exist
  data$lpo[data$lpo < 0] <- 0
  
  d <- data[c("lpo", "interval")]
  res <- do.call(algo.name, append(d, args.to.use))

  list(lpo = data$lpo, lpm = res, features = data$features)
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

meas.er.knn.diff <- function(data, job, instance, ...) {
  res <- priv.er.knn(lpo = diff(instance$lpo), lpm = diff(instance$lpm), ...)
  res
}

# feature mass

meas.fm.diff <- function(data, job, instance, ...) {
  res <- priv.fm(lpo = diff(instance$lpo), lpm = diff(instance$lpm), ...)
  res
}

# K-divergence

meas.kd.orig <- function(data, job, instance, ...) {
  res <- priv.kd(lpo = instance$lpo, lpm = instance$lpm, ...)
  -res
}

# KL divergence

meas.kl.hist.orig <- function(data, job, instance, ...) {
  res <- priv.kl.hist(lpo = instance$lpo, lpm = instance$lpm, ...)
  -res
}

meas.kl.hist.diff <- function(data, job, instance, ...) {
  res <- priv.kl.hist(lpo = diff(instance$lpo), lpm = diff(instance$lpm), ...)
  -res
}

meas.kl.knn.orig <- function(data, job, instance, ...) {
  res <- priv.kl.knn(lpo = instance$lpo, lpm = instance$lpm, ...)
  -res
}

meas.kl.knn.diff <- function(data, job, instance, ...) {
  res <- priv.kl.knn(lpo = diff(instance$lpo), lpm = diff(instance$lpm), ...)
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

meas.mi.hist.orig.bin <- function(data, job, instance, ...) {
  res <- priv.mi.hist(lpo = instance$lpo, lpm = instance$lpm, 
                      features = instance$features, regime = "iid", num.bins = 2, ...)
  res
}

meas.mi.hist.diff <- function(data, job, instance, ...) {
  res <- priv.mi.hist(lpo = diff(instance$lpo), lpm = diff(instance$lpm), features = instance$features, ...)
  res
}

meas.mi.hist.diff.bin <- function(data, job, instance, ...) {
  res <- priv.mi.hist(lpo = diff(instance$lpo), lpm = diff(instance$lpm), 
                      features = instance$features, regime = "iid", num.bins = 2, ...)
  res
}

meas.mi.knn.orig <- function(data, job, instance, ...) {
  res <- priv.mi.knn(lpo = instance$lpo, lpm = instance$lpm, features = instance$features, ...)
  res
}

meas.mi.knn.diff <- function(data, job, instance, ...) {
  res <- priv.mi.knn(lpo = diff(instance$lpo), lpm = diff(instance$lpm), features = instance$features, ...)
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


addAlgorithm(name = "meas.cs.orig", fun = meas.cs.orig)
addAlgorithm(name = "meas.cs.diff", fun = meas.cs.diff)
addAlgorithm(name = "meas.ce.orig", fun = meas.ce.orig)
addAlgorithm(name = "meas.ce.diff", fun = meas.ce.diff)
addAlgorithm(name = "meas.er.hist.diff", fun = meas.er.hist.diff)
addAlgorithm(name = "meas.er.knn.diff", fun = meas.er.knn.diff)
addAlgorithm(name = "meas.fm.diff", fun = meas.fm.diff)
addAlgorithm(name = "meas.kd.orig", fun = meas.kd.orig)
addAlgorithm(name = "meas.kl.hist.orig", fun = meas.kl.hist.orig)
addAlgorithm(name = "meas.kl.hist.diff", fun = meas.kl.hist.diff)
addAlgorithm(name = "meas.kl.knn.orig", fun = meas.kl.knn.orig)
addAlgorithm(name = "meas.kl.knn.diff", fun = meas.kl.knn.diff)
addAlgorithm(name = "meas.lv.orig", fun = meas.lv.orig)
addAlgorithm(name = "meas.mi.hist.orig", fun = meas.mi.hist.orig)
addAlgorithm(name = "meas.mi.hist.diff", fun = meas.mi.hist.diff)
addAlgorithm(name = "meas.mi.knn.orig", fun = meas.mi.knn.orig)
addAlgorithm(name = "meas.mi.knn.diff", fun = meas.mi.knn.diff)
addAlgorithm(name = "meas.rc.orig", fun = meas.rc.orig)
addAlgorithm(name = "meas.dc.orig", fun = meas.dc.orig)
addAlgorithm(name = "meas.dc.diff", fun = meas.dc.diff)
addAlgorithm(name = "meas.tvd.orig", fun = meas.tvd.orig)
addAlgorithm(name = "meas.mi.hist.orig.bin", fun = meas.mi.hist.orig.bin)
addAlgorithm(name = "meas.mi.hist.diff.bin", fun = meas.mi.hist.diff.bin)



#### check what is there

reg$problems
reg$algorithms

#### parameters (names of BBLHs to use and storage characteristics - inputs of BBLHs)

tmp1 <- CJ(ind = 1:4,
  algo.name = c("best.effort.moderate", "lazy.charging.moderate", 
                "lazy.stepping.moderate", "nill.moderate",
                "random.charging.moderate"))
           
tmp2 <- data.table(ind = 1:4,
            max.be = c(1.2, 4.0, 2.4, 6.5),
            min.be = c(0.0, 0.0, 0.0, 0.0),
            max.bc = c(0.27, 2.0, 4.0, 5.0),
            max.bd = c(0.27, 2.0, 4.0, 5.0),
            max.bp = c(0.27, 2.0, 4.0, 5.0),
            be = c(0.5, 0.5, 0.5, 0.5))

pars <- merge(tmp1, tmp2, by = "ind")[, -1]

pdes = list()
for(i in 1:length(all.lps)){
  pdes <- append(pdes, list(pars))
}
names(pdes) <- names(all.lps)

ades = list(
  meas.cs.orig = data.table(),
  meas.cs.diff = data.table(),
  meas.ce.orig = data.table(),
  meas.ce.diff = data.table(),
  meas.er.hist.diff = data.table(regime = c("zeros", "nozeros")),
  meas.er.knn.diff = data.table(regime = c("zeros", "nozeros")),
  meas.fm.diff = data.table(regime = c("fm", "rfm", "ed")),
  meas.kd.orig = data.table(),
  meas.kl.hist.orig = data.table(),
  meas.kl.hist.diff = data.table(),
  meas.kl.knn.orig = data.table(),
  meas.kl.knn.diff = data.table(),
  meas.lv.orig = data.table(),
  meas.mi.hist.orig = data.table(regime = c("iid", "ms", "mns")),
  meas.mi.hist.diff = data.table(regime = c("iid", "ms")),
  meas.mi.knn.orig = data.table(regime = c("iid", "ms")),
  meas.mi.knn.diff = data.table(regime = c("iid", "ms")),
  meas.rc.orig = data.table(regime = c("v", "w")),
  meas.dc.orig = data.table(regime = c("reg", "lpo")),
  meas.dc.diff = data.table(regime = c("reg", "lpo")),
  meas.tvd.orig = data.table(),
  meas.mi.hist.orig.bin = data.table(),
  meas.mi.hist.diff.bin = data.table()
)

addExperiments(pdes, ades, repls = 1)


#### check what is where

summarizeExperiments()
unwrap(getJobPars())

#### run experiments

submitJobs()
waitForJobs()


#### results

reduce <- function(res) res
results = unwrap(reduceResultsDataTable(fun = reduce))
pars = unwrap(getJobPars())
results = ijoin(pars, results)

results[max.be == 1.2, battery := "battery1"]
results[max.be == 2.4, battery := "battery2"]
results[max.be == 4, battery := "battery3"]
results[max.be == 6.5, battery := "battery4"]


save(results, file = paste0(getwd(), "/results.RData"))


#### results processing

# for different batteries different algorithms:
results[problem == 'CER6' & algorithm == 'meas.ce.orig' & battery == 'battery1',]
results[problem == 'CER6' & algorithm == 'meas.ce.orig' & battery == 'battery2',]
results[problem == 'CER6' & algorithm == 'meas.ce.orig' & battery == 'battery3',]
results[problem == 'CER6' & algorithm == 'meas.ce.orig' & battery == 'battery4',]

results <- results[!grepl(".knn.", results$algorithm), ]
results[algorithm == "cs.diff" & result.1 == -1, result.1 := 0]
results$algorithm <- gsub(".hist", "", results$algorithm)
results$algorithm <- gsub("meas.", "", results$algorithm)

pivot <- unique(results[, c("algorithm", "regime")])
pivot$meas.id <- paste(pivot$algorithm, pivot$regime, sep = ".")
pivot$meas.id <- gsub(".NA", "", pivot$meas.id)


for(i in 1: nrow(pivot)){
  print(i)
  for(j in unique(results$problem)){
    for(k in unique(results$battery)){
      if(is.na(pivot$regime[i])){
        ids <- results[problem == j & battery == k & algorithm == pivot$algorithm[i], job.id]
      } else {
        ids <- results[problem == j & battery == k & algorithm == pivot$algorithm[i] & regime == pivot$regime[i], job.id]
      }
      if(length(ids) != 5) stop("length of the vector is not 5")
      results[job.id %in% ids, algo.rank := rank(results[job.id %in% ids, result.1])]
      results[job.id %in% ids, meas.id := pivot$meas.id[i]]
    }
  }
}

library(reshape)
aggregate(results$algo.rank, list(results$algo.name), mean)

d <- aggregate(results$algo.rank, list(results$algo.name, results$battery), mean)
names(d) <- c("algo.name", "battery", "value")
res <- cast(d, algo.name~battery, mean)
res <- apply(res, 2, rank)
write.csv(res, paste0(getwd(), "/res_bat.csv"))
# best effort on battery1 and battery2

d <- aggregate(results$algo.rank, list(results$algo.name, results$problem), mean)
names(d) <- c("algo.name", "load.profile", "value")
res <- cast(d, algo.name~load.profile, mean)
res <- apply(res, 2, rank)
write.csv(res, paste0(getwd(), "/res_lp.csv"))
# nill in CER2 and in SmartB2

d <- aggregate(results$algo.rank, list(results$algo.name, results$meas.id), mean)
names(d) <- c("algo.name", "measure", "value")
res <- cast(d, algo.name~measure, mean)
res <- apply(res, 2, rank)
write.csv(res, paste0(getwd(), "/res_meas.csv"))

write.csv(results, paste0(getwd(), "/res_raw.csv"))

d <- results[algo.name == "best.effort.moderate",]
barplot(table(d$algo.rank)/sum(table(d$algo.rank)))

d <- results[algo.name == "nill.moderate",]
barplot(table(d$algo.rank)/sum(table(d$algo.rank)))

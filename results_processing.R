#### results processing

library(batchtools)
library(reshape)

reg <- loadRegistry(file.dir = paste0(getwd(), "/registry"), work.dir = getwd(), writeable = TRUE)

#' make a flat table with the results. Add a column "battery"
#' specifying which of the four batteries was used in the experiment.
#' save the resulting table.

reduce <- function(res) res
results = unwrap(reduceResultsDataTable(fun = reduce))
pars = unwrap(getJobPars())
results = ijoin(pars, results)

results[max.be == 1.2, battery := "battery1"]
results[max.be == 2.4, battery := "battery2"]
results[max.be == 4, battery := "battery3"]
results[max.be == 6.5, battery := "battery4"]

save(results, file = paste0(getwd(), "/results.RData"))

#' load the result

load(paste0(getwd(), "/results.RData"))

#' delete unnecessary information from measure's names

# results[algorithm == "cs.diff" & result.1 == -1, result.1 := 0]
results$algorithm <- gsub(".hist", "", results$algorithm)
results$algorithm <- gsub("meas.", "", results$algorithm)

#' combine the measure names and its variant into a single string,
#' so-called measure ID

pivot <- unique(results[, c("algorithm", "regime")])
pivot$meas.id <- paste(pivot$algorithm, pivot$regime, sep = ".")
pivot$meas.id <- gsub(".NA", "", pivot$meas.id)

#' for each measure ID, for each dataset and for each battery
#' rank the BBLH and store the respective rank values in a separate column

for(i in 1: nrow(pivot)){
  print(paste(i, "/", nrow(pivot)))
  for(j in unique(results$problem)){
    for(k in unique(results$battery)){
      if(is.na(pivot$regime[i])){
        ids <- results[problem == j & battery == k & algorithm == pivot$algorithm[i], job.id]
      } else {
        ids <- results[problem == j & battery == k & algorithm == pivot$algorithm[i] & regime == pivot$regime[i], job.id]
      }
      if(length(ids) != 7) stop("length of the vector is not 7")
      results[job.id %in% ids, algo.rank := rank(results[job.id %in% ids, result.1])]
      results[job.id %in% ids, meas.id := pivot$meas.id[i]]
    }
  }
}

#' merge the columms "algo.name" and "mode" into a single string,
#' which will identyfy the BBLH used. The result is in "algo.name" field

results$algo.name <- paste0(results$algo.name, ".", results$mode)

#' average ranks of the BBLHs

aggregate(results$algo.rank, list(results$algo.name), mean)

#' average ranks of the BBLHs with respect to battery

d <- aggregate(results$algo.rank, list(results$algo.name, results$battery), mean)
names(d) <- c("algo.name", "battery", "value")
res <- cast(d, algo.name~battery, mean)
res <- t(apply(res, 2, rank))
write.csv(res, paste0(getwd(), "/results/res_bat.csv"))

#' average ranks of the BBLHs with respect to load profile

d <- aggregate(results$algo.rank, list(results$algo.name, results$problem), mean)
names(d) <- c("algo.name", "load.profile", "value")
res <- cast(d, algo.name~load.profile, mean)
res <- t(apply(res, 2, rank))
write.csv(res, paste0(getwd(), "/results/res_lp.csv"))

#' average ranks of the BBLHs with respect to privacy measure

d <- aggregate(results$algo.rank, list(results$algo.name, results$meas.id), mean)
names(d) <- c("algo.name", "measure", "value")
res <- cast(d, algo.name~measure, mean)
res <- t(apply(res, 2, rank))
write.csv(res, paste0(getwd(), "/results/res_meas.csv"))

#' some additional plots

d <- results[algo.name == "alg.be.first",]
barplot(table(d$algo.rank)/sum(table(d$algo.rank)))

d <- results[algo.name == "alg.nill.na",]
barplot(table(d$algo.rank)/sum(table(d$algo.rank)))

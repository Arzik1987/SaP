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

results[max.be == 1.2, battery := "AC Battery"]
results[max.be == 2.4, battery := "Myreserve Pack"]
results[max.be == 4, battery := "eco 8.0/4"]
results[max.be == 6.5, battery := "RESU 6.5"]

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
results[algo.name == "alg.be.first", algo.name := "BE1"]
results[algo.name == "alg.be.second", algo.name := "BE2"]
results[algo.name == "alg.nill.na", algo.name := "NILL"]
results[algo.name == "alg.stepping.lc", algo.name := "LC"]
results[algo.name == "alg.stepping.ls1", algo.name := "LS1"]
results[algo.name == "alg.stepping.ls2", algo.name := "LS2"]
results[algo.name == "alg.stepping.rc", algo.name := "RC"]

#' save and load the result

save(results, file = paste0(getwd(), "/results_ranked.RData"))
load(paste0(getwd(), "/results_ranked.RData"))

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
res <- as.data.frame(t(apply(res, 2, rank)))
res$measure <- rownames(res)
res <- data.table(res)

res$diff <- ifelse(grepl("diff", res$measure), "y", "n")
res$measure <- gsub(".diff", "", res$measure)
res$measure <- gsub(".orig", "", res$measure)
res[measure == "ce", measure := "H"]
res[measure == "cs", measure := "CS"]
res[measure == "dc.reg", measure := "R^2_2"]
res[measure == "dc.lpo", measure := "R^2_p"]
res[measure == "er.nozeros", measure := "ER_{nz}"]
res[measure == "er.zeros", measure := "Cs"]
res[measure == "fm.ed", measure := "FM_{ed}"]
res[measure == "fm.fm", measure := "FM"]
res[measure == "fm.rfm", measure := "FM_r"]
res[measure == "kd", measure := "K"]
res[measure == "kl", measure := "KL"]
res[measure == "lv", measure := "LV"]
res[measure == "mi.bin", measure := "MI_b"]
res[measure == "mi.iid", measure := "MI_i"]
res[measure == "mi.ms", measure := "MI_s"]
res[measure == "mi.mns", measure := "MI_v"]
res[measure == "rc.v", measure := "RU_v"]
res[measure == "rc.w", measure := "RU_w"]
res[measure == "tvd", measure := "TVD"]
res$measure <- paste0("$\\mathop{", res$measure, "}$")

write.csv(res, paste0(getwd(), "/results/res_meas.csv"), row.names = FALSE)

#' some additional plots

d <- results[algo.name == "BE1",]
barplot(table(d$algo.rank)/sum(table(d$algo.rank)))

d <- results[algo.name == "NILL",]
barplot(table(d$algo.rank)/sum(table(d$algo.rank)))

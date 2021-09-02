
## Code for the experiments comparing Battery-Based Load Hiding Algorithms

The code uses a fixed version of the R package [bpriv](https://github.com/Arzik1987/bpriv), implementing smart meter privacy measures and several BBLH algorithms. This is done for reproducibility.

### About
This repository contains the code to reproduce our experiments. The main components are
* `bpriv_0.0.1.tar.gz` - the package containing the implementation of privacy measures and BBLH algorithms;
* `data_preparation.R` - does data preprocessing to obtain load profiles we use. It finds the files containing required load profiles, extracts the data for the required period, inputs missing measurement values and aggregates to the specified sampling rate. Runtime ~ 10 min;
* `experiment_measures.R` - does the experiments.  It applies Battery-Based Load Hiding (BBLH) algorithms with the parameters defined by characteristics of different storage systems to different load profiles. After it, for each pair of "user load" (load profile before applying BBLH algorithm) and corresponding "grid load" (load profile altered with BBLH algorithm), the code computes the values of different measures. It creates the new directory "registry" to store the experimental parameters and the raw results there. Runtime ~ 2 hours on four cores;
* `results_processing.R`- processes the raw experimental results and outputs their summaries in the form of `.csv` files in the directory "results" (created automatically). It also saves the raw results in the form of a flat table in the file `results_ranked.RData`, so that one can construct other summaries. Runtime ~ 10 min;
* `main.R` - executes the three above source files in the correct order;
* `data_preparation.csv` - contains the descriptions of load profiles to be used - the information for data preparation. This includes:
   * *Name* - the name of the load profile how it will appear in the output. Should contain one of the dataset names (`ECO`, `REDD`, `Smart` or `CER`) as its part;
   * *Sample_rate_seconds* - the *required* frequency of measurements;
   * *Length* - information field. Number of observations the resulting (preprocessed) load profile;
   * *No_residents*- information field. Number of residents in the corresponding household;
   * *Location* - the code supports either 'Europe' or 'US'. It is assumed that the voltage in the European grid is 220V and in the US grid is 110V. This parameter is the input for the NILL BBLH algorithm;
   * *Time* - information field;
   * *Day_of_week* - information field;
   * *Employed* - information field;
   * *Source_unit* - information field;
   * *Pathpart* - specifies the path within the directory `\data\XXX\` to find the required file. Here `XXX` is either one of the values {`ECO`, `REDD`, `Smart`, `CER`};
   * *Start_date* - the earliest date in the required user load;
   * *End_date* - the latest date in the required user load;
   * *Warm_up_days* - number of days preceding *Start-date* which are used for BBLH algorithm, but excluded from privacy measure calculation (see the paper);
   * *Source_sample_rate_seconds* - information field. Sample rate in the original dataset.
   
### Requirements
The experiments are implemented in programming language [R](https://cran.r-project.org/). We used additional packages listed below. In case your R installation does not contain these packages, the code from this repository will install their latest versions automatically.
* anytime (version 0.3.7)
* batchtools (version 0.9.12)
* cluster (version 2.1.0)
* data.table (version 1.12.8)
* devtools (version 2.2.1)
* entropy (version 1.2.1)
* FNN (version 1.1.3)
* infotheo (version 1.2.0)
* lubridate (version 1.7.4)
* purrr (version 0.3.3)
* reshape (version 0.8.8)
* wavelets (version 0.3-0.1)

### Before you start
In our experiments, we use the commonly used datasets which we do not own. That is why we request you to download those datasets from their source and put them in the subfolders of the folder `data`. Each of subfolders contains the file `FILES.txt` specifying the link to the source in the first row and containing the list of files to be put in the respective subfolder.

### How to execute the code
We suggest to use [RStudio](https://www.rstudio.com/). Any of the above packages can be installed manually with the command install.packages("<package-name>")
To install our package (`bpriv_0.0.1.tar.gz`), execute `install.packages("<path-to-the-file>", repos = NULL, type = "source")`.
To repeat our experiments:
* open the file `experiments.Rproj`;
* open the file `main.R` and execute it.

The result will appear in the folder `results`.
To launch part of code, select it and press `Ctrl + Enter`.


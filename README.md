
## Code for the experiments comparing Battery-Based Load Hiding Algorithms

### About
This repository contains the code to reproduce our experiments. The main components are
 - `data_preparation.R` - does data preprocessing to obtain load profiles we use. It finds the files containing required load profiles, extracts the data for the required period, inputs missing measurement values and aggregates to the specified sampling rate. Runtime ~ 10 min.
 - `experiment_measures.R` - does the experiments.  It applies Battery-Based Load Hiding (BBLH) algorithms with the parameters defined by  characteristics of different storage systems to different load profiles. After it, for each pair of "user load" (load profile before applying BBLH algorithm) and corresponding "grid load" (load profile altered with BBLH algorithm), the code computes the values of different measures. It creates the new directory "registry" to store the experimental parameters and the raw results there. Runtime ~ 2 hour on 4 cores.
 - `results_processing.R`- processes the raw experimental results and outputs their summaries in the form of `.csv` files in the directory "results" (created automatically). It also saves the raw results in the form of a flat table in the file `results_ranked.RData`, so that one can construct other summaries. Runtime ~ 10 min.
 - `main.R` - executes the three above source files in the correct order.
 - `data_preparation.csv` - contains the descriptions of load profiles to be used - the information for data preparation. This includes:
 *--Name* - the name of the load profile how it will appear in the output. Should contain one of the dataset names (`ECO`, `REDD`, `Smart` or `CER`) as its part
*-- Sample_rate_seconds* - the *required* frequency of measurements
*--Length* - information field. Number of observations the resulting (pre-processed) load profile
*--No_residents*- information field. Number of residents in the corresponding household
*--Location* - either Europe or US are supported by the code. It is assumed that the voltage in European grid is 220V and in the US grid is 110V. This parameter is the input for NILL BBLH algorithm.
*--Time* - information field
*--Day_of_week* - information field
*--Employed* - information field
*--Source_unit* - information field
*--Pathpart* - specifies the path within the directory `\data\XXX\` to find the required file. Here `XXX` is either one of the values {`ECO`, `REDD`, `Smart`, `CER`}
*--Start_date* - the earliest date in the required user load
*--End_date* - the latest date in the required user load
*--Warm_up_days* - number of days preceiding *Start-date* which are used for BBLH algorithm, but excluded from privacy measure calculation (see the paper)
*--Source_sample_rate_seconds* - information field. Sample rate in the original dataset
### Before you start
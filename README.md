# Dissertation Experiments

Chapter 4 of my PhD dissertation contains a number of large-scale simulation studies 
run on the [Sol cluster computing environment](https://asurc.atlassian.net/wiki/spaces/RC/pages/1640103978/Sol+Supercomputer) 
at Arizona State University.

This repository contains the scripts used to design, run, and analyze those simulation studies. 

## Simulation Script

The core of the experiments is contained in the `scripts/simulation_study.R` script. This file:

1. Accepts 16 command line arguments governing the nature of the data to be simulated and the analysis to be run
2. Runs a parallel simulation study with all available cores via the [`doParallel`](https://cran.r-project.org/web/packages/doParallel/index.html) R package
3. Stores every iteration of the simulation study in a CSV file.

This file can be run on its own from the command line, provided R is installed along with the following R package dependencies:

* [`here`](https://cran.r-project.org/web/packages/here/index.html) 1.0.1
* [`foreach`](https://cran.r-project.org/web/packages/foreach/index.html) 1.5.2
* [`doParallel`](https://cran.r-project.org/web/packages/doParallel/index.html) 1.0.17
* [`dbarts`](https://cran.r-project.org/web/packages/dbarts/index.html) 0.9.23
* [`grf`](https://cran.r-project.org/web/packages/grf/index.html) 2.2.1
* [`xgboost`](https://cran.r-project.org/web/packages/xgboost/index.html) 1.7.3.1
* [`rpart`](https://cran.r-project.org/web/packages/rpart/index.html) 4.1.19
* My personal XBART fork, [`andrewherren/XBART`](https://github.com/andrewherren/XBART/tree/dissertation-simulations), using the version indexed by `dissertation-simulations` tag

Regarding the version of R itself, I have run these on R 4.1.0 and 4.2.2, and I have not investigated how far back in R history they will run properly. The script does not use many novel features of the base R language / ecosystem, but it's possible that some of the package dependencies do (most of them are "latest available" at the time I ran these experiments). 

With R and the resulting dependencies set, you can run the "default parameters" version of this script from the `scripts` folder at the command line (at least on unix and linux) via:

```{bash}
Rscript --verbose simulation_study.R 
```

or you can specify the parameters yourself via command line arguments. They must be provided in the following order:

1. `dgp_num`: Which of the 5 data-generating processes (DGPs) specified in the script to simulate
2. `sample_size`: Number of random observations drawn
3. `num_variables`: Number of covariates drawn
4. `kappa`: The "noise to signal" ratio used to scale the outcome variance
5. `number_simulations`: How many iterations of the simulation study to run
6. `ate_true`: The "true average treatment effect", around which the $\tau(X)$ conditional average treatment effect function is centered
7. `n_prop_submodels`: Number of "submodels" of the propensity function, $\pi(X)$, to include as covariates in the "multiple propensity" XBCF model
8. `on_remote`: Whether or not the simulations are running locally (0) or on a remote cluster (1)
9. `datestamp`: A string of the form `YYYYMMDDHHMM` indicating the time and day that the simulation study was run (which will label the snapshot folder in `outputs/snapshots` to which the results are exported)
10. `estimated_propensities`: Whether the propensity function is assumed known (0) or must be estimated from $(Z, X)$ pairs (1)
11. `residualize_xbcf`: Whether the XBCF models will be trained on the original outcome (0) or an outcome with an initial model-based estimate of $\hat{Y}$ subtracted out (1)
12. `project_pi_yhat`: Whether the propensity score is to be projected onto a model-based estimate of $\hat{Y}$ (1) or not (0)
13. `script_iteration`: Which "run" of a job is this (defaults to 1 if running as a single job from the command line)
14. `n_yhat_submodels`: Number of "submodels" of $\mathbb{E}(Y \mid X)$ to estimate on subsets of the variables in $X$.
15. `use_yhat_covariate`: Whether or not the marginal $\hat{Y}$ estimates should be used as covariates in the XBCF models (1) or not (0)
16. `grf_default`: Whether GRF's $\hat{Y}$ and $\hat{\pi}$ models are estimated using GRF's default procedures (1), or are supplied as XBART / xgboost estimates (0). [*The propensity score is only estimated in either case if estimated_propensities = 1.*]

## Job Submission Script

Using a cluster computing environment to run these simulations requires a bit of extra infrastructure. There are several workload management systems designed for large clusters, and ASU's Sol cluster uses [SLURM](https://en.wikipedia.org/wiki/Slurm_Workload_Manager). SLURM is directed by scripts with `#SBATCH` commands governing the resources requested, time allocated, and I/O details. 

These experiments rely on a series of python scripts which, based on an array of 
simulation parameters, create job files with `#SBATCH` scripts on the fly and submit them to the scheduler. I owe an enormous debt of gratitude to 
[Francisco Castillo](https://github.com/fjcasti1/), my friend and colleague whose original job submission python script forms the basis of the scripts in this repo.

The scripts are organized into "rounds" based on the layout of the dissertation.

1. `job_submission_remote_round1.py`: This runs the jobs and simulations underpinning the first round of simulation results. In this case, XBCF and GRF are both trained more or less using their default specifications, with both known and estimated propensities.
2. `job_submission_remote_round2.py`: This runs the jobs and simulations underpinning the second round of simulation results. In this case, GRF is trained using its default specification while the XBCF methods are augmented with $\hat{Y}$ covariates, with both known and estimated propensities.
3. `job_submission_remote_round3.py`: This runs the jobs and simulations underpinning the third round of simulation results. In this case, the XBCF methods are augmented with $\hat{Y}$ covariates and GRF is trained using XBART $\hat{Y}$ estimates and xgboost propensity estimates (except in the case where propensities are assumed known, which is one of the scenarios of this study).

They are run from the scripts folder on a remote cluster environment via:

```
python3 job_submission_remote_round1.py

python3 job_submission_remote_round2.py

python3 job_submission_remote_round3.py
```

Run at least one minute apart to ensure the datestamps are unique.

## Quarto Report Generation

For a quick view of the results, there are quarto files that generate reports for each of the job "rounds" in the previous section. The reports are parameterized with a datestamp, and can be run from the command line by setting a `DATESTAMP` environment variable:

```
DATESTAMP="202303280203"

quarto render simulation_analysis_round1.qmd --output "../reports/simulation_analysis_round1_$DATESTAMP.pdf" --to pdf -P datestamp="$DATESTAMP"
```

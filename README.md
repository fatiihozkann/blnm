# Bayes Factors for Conditional Independence in Latent Network Modeling (BLNM)

R code for the paper *Bayes factor for conditional independence in latent network modeling* by Fatih Ozkan (Baylor University).

OSF project: <https://osf.io/p7r2m>

The Bayesian latent network model (BLNM) combines a confirmatory factor measurement model with a Gaussian graphical model on the latent variables. It computes inclusion and exclusion Bayes factors for each latent partial correlation by enumerating network structures and comparing their marginal likelihoods, so conditional independence between latent variables can be tested directly at the latent level, with quantifiable evidence for absence as well as presence.

## Contents

The layout is flat on purpose: `bfi_analysis.R` and `study1_run.R` call `source("blnm_gen.R")` with a relative path, so run everything from the repository root.

- `blnm_gen.R` : Core implementation for an arbitrary number of latent variables M. Model configuration, structure enumeration, marginalized likelihood, positive definiteness truncation constants (exact and Monte Carlo), Laplace approximation, MCMC (adaptive random walk Metropolis and independence Metropolis with a multivariate t proposal), bridge sampling, and inclusion Bayes factors under uniform and beta binomial structure priors. Sourced by the two scripts below.
- `blnm_poc.R` : Self contained proof of concept with M = 3 factors on the Holzinger and Swineford (1939) data that ship with lavaan. Full enumeration of all 8 structures, comparison with a two step factor score approach, a likelihood ratio test, and a Savage-Dickey check. Optional slab width sensitivity analysis and synthetic data validation are controlled by flags near the top. Set the environment variable `BLNM_SMOKE=1` for a fast reduced run.
- `bfi_analysis.R` : Big Five empirical example on the `bfi` data that ship with psych (25 items, N = 2436 complete cases, M = 5, 1024 candidate structures). Uses the two stage scheme from the paper: a Laplace sweep over all structures, then bridge sampling refinement of the top set. Reports inclusion Bayes factors under both structure priors, posterior structure probabilities, and the model averaged network. Requires `blnm_gen.R` in the same folder.
- `study1_run.R` : The simulation study reported in the paper (18 conditions, 500 replications each). Runs in parallel on all available cores minus one, writes a checkpoint (`study1_checkpoint.rds`) as it goes, and resumes automatically if interrupted. Requires `blnm_gen.R`. The environment variable `STUDY1_NREP` overrides the replication count for smaller test runs.
- `study1_analysis.R` : Builds the summary tables and the three panel figure from the simulation output. Reads `study1_results.rds` (see below).
- `study1_results.csv` : Results of the reported 500 replication run, included so the tables and figure can be regenerated without rerunning the simulation.
- `figures/` : Output figures are written here; the published figure files can also be stored here.

## Requirements

R (version 4.0 or later) with these packages:

```r
install.packages(c("mvtnorm", "coda", "bridgesampling", "lavaan", "psych", "logspline"))
```

Optional: `Rcpp` (with a working C++ toolchain) speeds up the Monte Carlo positive definiteness constants used when M is 4 or more. If Rcpp compilation is unavailable, the code detects this and falls back to a pure R implementation automatically.

Both datasets (`HolzingerSwineford1939`, `bfi`) ship with lavaan and psych, so no data downloads are needed. The only data file in this repository is the simulation output.

## Running

From the repository root:

```sh
Rscript blnm_poc.R          # minutes; BLNM_SMOKE=1 Rscript blnm_poc.R for a quick check
Rscript bfi_analysis.R      # the Laplace sweep is fast, the bridge refinement stage takes hours
Rscript study1_run.R        # long; checkpointed, safe to interrupt and resume
Rscript study1_analysis.R   # seconds; needs study1_results.rds
```

`study1_analysis.R` expects `study1_results.rds`. To create it from the CSV shipped in this repository:

```sh
Rscript -e 'saveRDS(read.csv("study1_results.csv"), "study1_results.rds")'
```

For a reduced simulation run, for example 20 replications per condition:

```sh
STUDY1_NREP=20 Rscript study1_run.R
```

## Reproducibility

All random seeds are fixed in the scripts (each simulation replication uses seed `100000 * cell + replication`). Small Monte Carlo variation across platforms and package versions is expected in the bridge sampling estimates.

## Citation

If you use this code, please cite the accompanying paper. A preprint is available through the OSF project at <https://osf.io/p7r2m>; citation details will be updated here once the preprint DOI is live.

## License

MIT, see `LICENSE`.

## Acknowledgment


Contact: Fatih Ozkan, Baylor University, fatih_ozkan1@baylor.edu

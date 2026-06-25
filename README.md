# Multi-type Branching Inference on Contact Trees

Code, simulation pipeline, and data for the manuscript **"Multi-type branching inference on contact trees with application to COVID-19."** The framework infers epidemiological parameters from *contact-traced transmission trees* (reported who-infected-whom events), not from pathogen genome sequences. The mathematical ingredient is adapted from phylodynamics (BiSSE-style backward Kolmogorov ODEs), but the object being modelled is a contact tree rather than a phylogeny.

The central object is an augmented state space `(i, k)` in which `i` counts the secondary infections an individual has already produced and `k` is its contact degree. As `i` grows, the remaining transmission potential `(k - i)·β` declines, which is what allows the likelihood to read information out of the branching pattern of the tree. The estimand throughout is the basic reproduction number `R0 = k·β/λ` with `λ = μ + σ` fixed.

---

## What the repository does

1. **Simulate** epidemics on contact networks and extract the sampled subtree as an edge table (Python, `Tree Simulations/`).
2. **Validate** the backward ODEs for the extinction probability `E_i(t)` and the tip density `D^i_j(t)` against the simulation, for both fixed and random contact degree (R, `ODE vs Theory/`).
3. **Estimate** `R0` and `k` by maximum likelihood on simulated trees, including identifiability and sensitivity analyses (R, `MLE/`).
4. **Estimate** the same quantities (simulated trees) in a Bayesian setting with latent branching times and a partial-resolution study (R, `Bayesian/`).
5. **Apply** the method to the Karnataka COVID-19 contact-tracing data (R + data, `COVID-19-estimation/`).

---

## Requirements

**Python** (simulation): `numpy`, `pandas`, `tqdm`. Standard library: `random`, `bisect`, `math`, `os`.

**R** (everything else): `deSolve` (ODE integration), `splines`, `readxl` (Excel ingest in the COVID analysis).

No fixed seeds are baked into the estimation scripts beyond the `SEED` constants exposed at the top of each runner; true parameter values are never hardcoded into the likelihood and are recomputed dynamically.

---

## Repository layout

```
Multi-type-branching-inference-on-contact-trees/
├── Tree Simulations/                    # Python epidemic + subtree simulator
│   ├── simulation-epi-tree-structure-n-tip.py
│   └── simulation-epi-tree-structure-one-tip.py
│
├── ODE vs Theory/                       # ODE-vs-simulation validation (R)
│   ├── Ei-and-Dij-sim-vs-ode-fixed-degree.R
│   ├── Ei-and-Dij-sim-vs-ode-random-degree.R
│   ├── Simulated data/
│   │   ├── phylo-epi-sim-data-fixed.csv
│   │   └── phylo-epi-sim-data-Poisson.csv
│   └── Ei.* , Di0.* … Di4.* , Di0_random.* … Di4_random.*   (figures)
│
├── MLE/                                 # Maximum-likelihood estimation (R)
│   ├── parameter_estimate_mle.R
│   ├── parameter_estimate_mle_pobs.R
│   ├── mle_results_summary.R
│   ├── full_tree_edges.csv              # simulated edge table (input)
│   ├── mle_results_table.csv            # results (output)
│   ├── ro_vs_pobs_sensitivity_results.csv
│   └── ro_vs_k.* , pi.* , ro_vs_pobs.*  (figures)
│
├── Bayesian/                            # Bayesian / partial-resolution MCMC (R)
│   ├── parameter_estimate_bayesian_unknown_internals_unknown_k.R
│   ├── partial_resolved_tree.R
│   ├── run_partial_resolved_tree.R
│   ├── partial_tree_all_results.rds
│   └── mcmc_posterior_vs_R0.* , mcmc_post_vs_k.* ,
│       mcmc_ro_vs_pres.* , mcmc_k_vs_pres.*        (figures)
│
└── COVID-19-estimation/                 # Karnataka application (R + data)
    ├── covid_data_analysis.R
    ├── covid_data_prep.R
    ├── parameter_estimation_mle.R
    ├── covid_mcmc.R
    ├── run_covid_mcmc.R
    ├── covid_edges.csv                  # built edge table (output of prep)
    ├── covid_mcmc_chain.csv             # posterior draws (output)
    ├── covid_mcmc_results.rds           # full MCMC result object (output)
    ├── covid-mcmc-*.{eps,pdf}           (figures)
    ├── covid-19-data/
    │   ├── ann2.csv                     # full line list (71,068 records)
    │   ├── contacts.csv                 # directed transmission links
    │   ├── traced.csv                   # contact-traced subset
    │   ├── untraced.csv                 # surveillance-detected subset
    │   ├── 54_serial_interval_data.xlsx
    │   └── Delays_3 for histogram.xlsx
    └── Paper/
        ├── gupta.pdf                    # source study
        └── pone.0270789.s001.docx       # source supplementary material
```

Each figure is provided in both `.eps` (for LaTeX inclusion) and `.pdf` (for preview). Figures are saved manually from the interactive R plots; the scripts draw them on screen rather than writing image files.

---

## `Tree Simulations/` — the simulator

A discrete-event stochastic SIR-type simulator on a contact structure. A focal individual is tagged at a burn-in time, its clade is followed forward, and the *sampled* subtree (the observed, `σ`-detected tips and their ancestors) is extracted into an edge table whose rows feed directly into the likelihood.

The simulation tracks each individual's contact degree `k`, infection time, recovery type (unobserved `μ` removal versus observed `σ` detection), and the running state `phyState = (k, j)` where `j` is the number of downstream infections produced so far.

**`simulation-epi-tree-structure-n-tip.py`** — the production simulator used for the full-tree estimation work. It supports four contact-degree regimes via a `DegreeDistribution` class: `fixed` (every individual has exactly `k`), `poisson`, `negbinom` (mean and dispersion, variance `k + k²/dispersion`), and `geometric`. The root focal individual draws its degree once at infection; every newborn draws its own degree independently from the same distribution. For `fixed`, the behaviour collapses exactly to a constant `k`.

The `SampledSubtree` class walks from each observed tip up to the tagged root, segments each individual's lifetime at its within-individual transmission events, and emits one row per segment. The edge-table schema is:

```
rep_id, edge_type, lineage_type, par_id, chi_id, tau_a, tau_b, delta,
s, c, j_obs, k, branch_rate, ind_id, infector_id, seg_idx, n_segs,
n_tips, n_internal, root_state, k_root, t_root_abs
```

Here `tau_a`/`tau_b` are the absolute start/end times of a segment, `delta` the elapsed duration, `s`/`c` the state entering/leaving the segment, `j_obs` the observed state at a sampled tip (and `NaN` otherwise), `k` the *per-individual* degree, and `branch_rate = (k - s)·β` at internal branching events (`NaN` elsewhere). The `k` column varies row to row under random-degree distributions; under `fixed` it is constant.

The script also contains `verify_edge_table`, a battery of consistency checks (tip rows carry a finite `j_obs`, segment times are contiguous, newborns enter at state `0`, last and newborn segments carry no branch rate, `k ≥ j_obs`, and so on), and `log_likelihood_from_df`, a reference implementation of the fixed-`k` likelihood used to cross-check the R estimators.

Run it from the command line with a distribution name, for example `python simulation-epi-tree-structure-n-tip.py fixed`. The output is `full_tree_edges_<dist>.csv`. The committed `MLE/full_tree_edges.csv` (753 valid replicates at `k=4, β=1.5, μ=σ=0.5`, hence `R0=6`, `p_obs=0.5`) is the fixed-degree product of this script.

**`simulation-epi-tree-structure-one-tip.py`** — an earlier single-observation variant. Rather than emitting a full edge table, it records for each replicate the focal state and time, the single observed state `(j_obs, k_obs)`, and the clade lifespan, distinguishing the three terminal cases: tagged-and-observed, extinct, and never-tagged. It also carries scaffolding for forward, backward, and recursive contact tracing (`tracing` flags) that is not exercised in the main analysis. Its output is `phylo-epi-sim-data-<degree>.csv`, which is exactly the input consumed by the ODE-validation scripts in `ODE vs Theory/`.

---

## `ODE vs Theory/` — analytic ODEs versus simulation

These two scripts confirm that the backward ODE system reproduces the empirical distributions produced by the simulator. Both bin the simulated clade lifespans by time to estimate, empirically, the survival fraction `E_i(t)` and the tip-density `D^i_j(t)`, then overlay the ODE solution computed with `deSolve`.

**`Ei-and-Dij-sim-vs-ode-fixed-degree.R`** — the fixed-degree case. It reads `Simulated data/phylo-epi-sim-data-fixed.csv`, conditions on each focal degree `i_foc`, and solves the coupled `E`/`D` system

```
dE_i = μ − (μ + σ + (k−i)β) E_i + (k−i)β · E_1 · E_{i+1}
dD_i = − (μ + σ + (k−i)β) D_i + (k−i)β (E_1 D_{i+1} + E_{i+1} D_1)
```

with `D^{j}_{j}(0) = σ` and `E_i(0)=1`. The helper functions `Ei()`, `Dij()`, and `plot_all_Ei()` overlay ODE (line) and simulation (points). These produce `Ei.*` and `Di0.*` through `Di4.*`.

**`Ei-and-Dij-sim-vs-ode-random-degree.R`** — the random-degree (Poisson-mixed) case. It decouples the focal/latent degree `focal_k` from the Poisson mixing mean `lambda`, reads `Simulated data/phylo-epi-sim-data-Poisson.csv`, builds the newborn mixing weights `w_k = dpois(0:Kmax, lambda)` truncated at `Kmax`, and solves the corresponding mixed system. These produce the `Di0_random.*` through `Di4_random.*` figures.

The two simulated-data CSVs (`phylo-epi-sim-data-fixed.csv`, `phylo-epi-sim-data-Poisson.csv`, share the schema `i_foc, k_foc, t_foc, j_obs, k_obs, t_obs` and are the output of the one-tip simulator.

---

## `MLE/` — maximum-likelihood estimation

The likelihood is assembled bottom-up over the edge table. Each row contributes one of:

- a tip-density term `D^{s}_{j_obs}(δ)` when the segment is the last for a sampled individual,
- an extinction term `E_s(δ)` when the last segment belongs to an unobserved sub-clade,
- an interval-survival term `exp(−(λ + (k−s)β)·δ)` for non-terminal segments, and
- a branching-rate factor `(k−s)·β` at each internal branching event.

The root is conditioned on producing at least one observed tip via `log π_{i0} − log(1 − E_{i0}(Δ))`, where `π` is the equilibrium state frequency obtained from the characteristic equation. With `λ` fixed, the free parameters are `R0` (continuous) and `k` (integer grid); `β` is recovered as `β = R0·λ/k`.

**`parameter_estimate_mle.R`** — the core estimator. Solves the `E`/`D` ODEs with `lsoda`, precomputes interpolators for every `j_obs`, evaluates the per-replicate log-likelihood, and optimises `R0` by 1-D `optimize` over a grid of `k` values, reporting AIC and BIC per model and a profile-likelihood plot with a 95% interval. Reads `full_tree_edges.csv`.

**`parameter_estimate_mle_pobs.R`** — a sensitivity analysis. Re-estimates `R0` while fixing the sampling probability `p_obs` to each value on a grid from 0 to 1, using the same simulated data (true `p_obs = 0.5`). It demonstrates that the sampling probability cannot be read from tree shape alone and must be fixed externally; with `λ` fixed the profile in `p_obs` is monotone. Writes `ro_vs_pobs_sensitivity_results.csv` and produces `ro_vs_pobs.*`.

**`mle_results_summary.R`** — the reporting driver. Sources the estimator, runs it across `k = 1…12`, and produces the consolidated outputs: the results table with confidence intervals and AIC weights (`mle_results_table.csv`), the `R0`-versus-`k` plot with intervals (`ro_vs_k.*`), and the equilibrium-frequency comparison theory-versus-simulation (`pi.*`).

**Data and result files.** `full_tree_edges.csv` is the input edge table. `mle_results_table.csv` reports the fit at each `k`: the best model sits at `k=4` with `R0 ≈ 6.04` (true 6), and the near-constant product `k·β̂` across the grid is the visible signature of partial non-identifiability between `k` and `β`. `ro_vs_pobs_sensitivity_results.csv` holds the `p_obs` sweep.

---

## `Bayesian/` — partial-resolution MCMC

Applies the Bayesian method to simulated data in which internal branching times, their states, and the contact degree `k` are jointly latent. The sampler alternates three Metropolis-Hastings moves: a Gaussian random walk on `log R0`, a uniform proposal for one latent branching time, and a discrete proposal for `k`. States are analytically marginalised over `s ∈ {0,…,k}` at each latent node, so only continuous and integer unknowns are sampled. The backward ODE is solved at most twice per iteration and cached.

Priors: `R0 ~ LogNormal(log 5, 1)`, `k ~ DiscreteUniform{1,…,K_MAX}`, and each latent time uniform on its feasible interval.

**`parameter_estimate_bayesian_unknown_internals_unknown_k.R`** — the likelihood engine, sourced by the sampler. It is the Bayesian-side counterpart of `MLE/parameter_estimate_mle.R` and holds the shared machinery: the `E`/`D` ODE solver (`solve_ED`), the interpolator precomputation (`precompute_ED`), the equilibrium frequencies (`compute_pi`), the per-replicate log-likelihood (`log_lik_rep`), the aggregate negative log-likelihood (`neg_log_lik`), and the grid-over-`k` estimator and profile plot. Crucially it also exposes the two cached helpers the partial-resolution sampler relies on, `precompute_ED_for_R0` (one ODE solve per `(R0, k)` proposal, returned as an `ED_cache`) and `ll_from_ED` (evaluate a replicate against a cached solution). It documents the same estimand distinction as the COVID likelihood, `R0 = k·β/λ` versus the derived `R̄0 = (β/λ)·S_π(k) < R0`.

**`partial_resolved_tree.R`** — the sampler and the partial-resolution study. It sources the engine above and defines the edge-mutation helpers used to move a single latent time, `mutate_edge`, `build_edges`, and `ll_rep_cached`, which reuse the cached `ED` solution so that a time move costs no new ODE solve. A fraction `p_res` of internal branching times is held fixed at its true value and the remainder is left latent; sweeping `p_res` from 0 (all latent) to 1 (all fixed) quantifies how much knowing the branching times tightens the `R0` posterior. It produces the four MCMC figures: posterior of `R0` (`mcmc_posterior_vs_R0.*`), posterior of `k` (`mcmc_post_vs_k.*`), and the two `p_res` sweeps (`mcmc_ro_vs_pres.*`, `mcmc_k_vs_pres.*`).

**`run_partial_resolved_tree.R`** — the runner. Exposes a fast-test block and a production block (`N_REPS_USE`, `N_ITER`, `N_BURN`, `THIN`, `P_RES_VALS`, `ODE_DT`, `SEED`) and then sources the sampler. The aggregate result across `p_res` levels is saved to `partial_tree_all_results.rds`, with per-level chains optionally saved alongside. `partial_tree_all_results.rds` also lets you reload and replot the saved results without re-running the chains.

---

## `COVID-19-estimation/` — Karnataka application

Applies the method to early-outbreak (March to May 2020) Karnataka COVID-19 contact-tracing data. State `i` is the number of secondary transmissions an individual had made by detection; `j_obs = children_pri` is that individual's state at confirmation. Index cases (`parentid = 0`) are surveillance-detected and enter the `σ` tip-density term; contact-traced cases (`parentid > 0`) are not `σ`-detections and contribute extinction terms, with their confirmation time providing a feasibility bound on their parent's transmission time. Infection (birth) times are latent and proposed by the MCMC. Calendar time is rescaled by the mean serial interval (about 5.65 to 5.88 days per model unit, taken from the serial-interval data).

**`covid_data_analysis.R`** — reproduces the Section 5 descriptive numbers and the month filter from the raw line list and contact links, splits the data into traced and untraced subsets, and writes `covid-19-data/untraced.csv`. It also ingests the serial-interval and delay spreadsheets. Run this first, since it produces the `untraced.csv` consumed downstream.

**`covid_data_prep.R`** — builds the model edge table `covid_edges.csv` from `untraced.csv` and `contacts.csv`. It maps each individual to its observed state and detection time, sets the latent birth-time bounds, and emits the schema `rep_id, ind_id, par_id, is_tip, s, c, j_obs, tau_a, tau_b, delta, seg_idx, n_segs, branch_rate, n_tips, root_state, lo_born, hi_born`.

**`parameter_estimation_mle.R`** — the COVID-side likelihood, a sibling of `MLE/parameter_estimate_mle.R` adapted to latent birth times and an ED cache. Its header documents the distinction between the estimand `R0 = k·β/λ` and the derived equilibrium-weighted `R̄0 = (β/λ)·S_π(k) < R0`, which is reported but never substituted for `R0`.

**`covid_mcmc.R`** — the three-step sampler for the real data: random walk on `log R0`, uniform proposal for one latent time (a birth time or an internal branching time), and a random walk on `k` under a Negative-Binomial prior chosen to reflect COVID-19 superspreading overdispersion. It must be sourced from the runner, not run standalone. Saves `covid_mcmc_results.rds` (full result object including intervals, Geweke and ESS diagnostics, and the `k` mode) and `covid_mcmc_chain.csv` (the posterior draws).

**`run_covid_mcmc.R`** — the entry point that wires the pipeline together: build the edge table (if absent), load the likelihood, validate the edge table, take a quick MLE for the `R0` starting value, set the Negative-Binomial prior on `k` (`NBINOM_MU`, `NBINOM_PHI`, `K_MAX`, `K_INIT`, `STEP_K`) and the fixed `p_obs ≈ 0.745` (the empirical index-case fraction), then run the sampler. Has fast-test and production parameter blocks.

**Figures.** `covid-mcmc-ro.*` (R0 posterior), `covid-mcmc-joint-posterior.*` (joint R0/k), `covid-mcmc-mu-K-distribution.*` (k posterior), and `covid-mcmc-ACF-ro.*` (autocorrelation diagnostic) are saved manually from the sampler's plots.

### `covid-19-data/` and `Paper/`

The `covid-19-data/` folder (directly under `COVID-19-estimation/`) holds the raw and derived inputs:

- `ann2.csv` — the full line list (71,068 records) with confirmation date, parent id, primary/secondary children counts, cluster id, and category (`cat_4`).
- `contacts.csv` — directed transmission links (`from, to, Reason`).
- `traced.csv` / `untraced.csv` — the contact-traced and surveillance-detected subsets produced by the analysis script.
- `54_serial_interval_data.xlsx` and `Delays_3 for histogram.xlsx` — the serial-interval and reporting-delay data used to set the time scale.

The `Paper/` folder holds the source material: `gupta.pdf` (the Karnataka study) and `pone.0270789.s001.docx` (its supplementary material).

---

## Reproduction order

1. `Tree Simulations/` → generate `full_tree_edges.csv` (full-tree) and the `phylo-epi-sim-data-*.csv` files (one-tip).
2. `ODE vs Theory/` → confirm the ODEs match the simulation for fixed and random degree.
3. `MLE/` → `source("mle_results_summary.R")` for the full-tree MLE, intervals, and the `p_obs` sensitivity.
4. `Bayesian/` → `source("run_partial_resolved_tree.R")` for the partial-resolution study (it sources the likelihood engine automatically).
5. `COVID-19-estimation/` → `covid_data_analysis.R`, then `source("run_covid_mcmc.R")` (which calls `covid_data_prep.R` if needed).

---

## Data provenance and citation

The COVID-19 data derive from the Karnataka contact-tracing study (`Paper/gupta.pdf` and its supplementary material). Please cite that source for the data and the present manuscript for the inference method.

# Minority-group incubators and majority-group reservoirs for climate change adaptation

This repository contains code and supporting documentation for the agent-based model analyzed in our paper, "Minority-group incubators and majority-group reservoirs for climate change adaptation".


## Quick start

To get started, clone this repository, e.g., execute the following in the
terminal: 

```
git clone https://github.com/eehh-stanford/MinMaj-Adaptation-Diffusion.git
```

After cloning the repository, install all dependencies by first starting the
[Julia REPL](https://docs.julialang.org/en/v1/stdlib/REPL/) then run

```
julia> using Pkg; Pkg.activate("."); Pkg.instantiate()
```

## Unit tests

We developed our model using test-driven development, which uses small,
executable code snippets to confirm the model works as expected and to document
model mechanics; see
[`src/test/model.jl`](https://github.com/eehh-stanford/MinMaj-Adaptation-Diffusion/blob/main/src/test/model.jl) to view the test suite.

While still in the REPL, run the unit tests to make sure all is working well:

```
julia> include("src/test/model.jl")
```

This should print two "Test Summary" outputs where all tests are shown to pass.
The tests initialize specially-initialized models and checks that model outputs
are as expected. 

# Run the model and analyze results

## Model and computational experiments

The model is implemented in [`src/model.jl`](src/model.jl) and the computational experiments that run the model over all parameter settings for the desired number of trials and used by the Slurm scripts (below) is in [`src/experiment.jl`](src/experiment.jl).

## Run all simulations on Slurm cluster

To run simulations on a Slurm cluster, log in to the cluster then execute the following commands from the project directory, first
```
./scripts/slurm/main.sh
```
to run the main analyses, and
```
./scripts/slurm/supplement.sh
```
to run the supplemental analyses. This creates a fresh, distinct version of simulation results that can be analyzed as we explain below, using archived data of the simulations used to create our results in the submitted version of the paper.

## Analysis

With the addition of a network-theoretic modeling option there are two forms of
analysis: (1) processing and calculating success rates and time to fixation from
ABM output data, and (2) analyzing the frequency of network ties given different
parameter settings.

### Success rates and times to fixation

To make the jitterplot and means of success rate with _Start group_ condition
on the x-axis included in Figure 1 of the 2024 manuscript, use the
[`success_over_groups_jitter`](https://github.com/eehh-stanford/MinMaj-Adaptation-Diffusion/blob/main/scripts/plot.R#L113) function. 

Use `main_asymm_heatmaps` to create the main heatmap results of _success rate_
as a function of $h_\mathrm{min}$ and $h_\mathrm{maj}$, which can be found in
[`scripts/plot.R`](https://github.com/eehh-stanford/MinMaj-Adaptation-Diffusion/blob/main/scripts/plot.R#L72).
For creating the heatmaps of average time to model fixation, pass the keyword
argument `measure = "step"` to `main_asymm_heatmaps`. Similarly, to create
supplemental analyses use the `supp_asymm_heatmaps` function in
[`scripts/plot.R`](https://github.com/eehh-stanford/MinMaj-Adaptation-Diffusion/blob/main/scripts/plot.R#L15).

To create the heatmaps you need the output data from the simulations presented in our journal article, stored in the `data` folder in the root project directory. To get the data in the right place, first create a `data` directory, then download and unzip the two zip files in our OSF repository: https://osf.io/cd9hx/. 

To create time series of individual model runs, use the
`make_all_group_prevalence_comparisons` function in
[`scripts/analysis.jl`](https://github.com/eehh-stanford/MinMaj-Adaptation-Diffusion/blob/main/scripts/analysis.jl#L290).


### Network tie frequencies

To understand the minority-group incubator and majority-group reservoir effects,
we calculated the probability that a minority group member teaches
at least one majority group member (currently Equation 6) and the probability
that a majority group member teaches another majority group member (currently 
Equation 9). These probabilities are functions of minoirty group fraction, 
$m$, and majority group homophily, $h_{\mathrm{maj}}$. We inspect plots of
each partial function individually. Plots over $m$ are made with the
`plot_over_m` function, and plots over $h$ are made with the `plot_over_h`
function. Both functions are contained in the file 
[`scripts/graph_counting.R`](/scripts/graph_counting.R).

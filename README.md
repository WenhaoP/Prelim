# Research Prelim Exam, STAT 572 in Spring 2025 at UW

This repository stores the code and data for reproducing the numerical results in the paper "Extremes on river networks” written by Peiman Asadi, Anthony C. Davison, and Sebastian Engelke and published on The Annals of Applied Statistics in 2015.

## Structure

* `Codes/`
    * `Functions.R`: helpful functions constitue numerical algorithms in other code files.
    * `Main.Rmd`: code for validating Asadi et. al's results and generating the figures in the final report of STAT 572.
    * `MultivariateFitting.R`: code for fitting the max-stable process $\eta_\Gamma$ on the data.
    * `UnivariateAnalysis.R`: code for estimating the marginal transformations $U_j$'s and marginal GEVDs.
* `Plots/`: generated figures.
* `Data/`: Application data.
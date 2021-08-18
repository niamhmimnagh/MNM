## MNM Models - Bayesian Multi-Species N-Mixture Models for Unmarked Animal Communities


### Simulation Studies Files
The following files contain code which was used in the simulation studies, to assess the estimation of parameters in the MNM models.
* mnm.R - simulate data for, and fit the MNM model,
* mnm_ar.R - simulate data for, and fit the AR model which contains an autoregressive term in the abundance,
* mnm_hurdle.R - simulate data for, and fit the hurdle model, which allows us to deal with data that contains large numbers of zero counts,
* mnm_hurdle_ar.R - simulate data for, and fit the hurdle-AR model, which is a combination of the AR model and hurdle model.


### Case Study Files
The following files all pertain to the model, fitted to the North American Breeding Bird Survey data, which was deemed best fit by DIC value i.e. the Hurdle-AR model, which contains a lat/long response surface on abundance, and which has detection probability that varies by site and species.

* N.RData contains estimated latent abundances for each species, at each site and year,
* sites.RData contains latitude and longitude for the sites of interest in Alaska,
* p.RData contains the detection probabilities for each species, at each site,

* abundance_plots.Rmd can be used in conjunction with N.RData and sites.RData to generate a html file of maps of Alaska, which contain the abundances of each species at each site and year,

* detection_probability_plots.Rmd can be used in conjunction with p.RData and sites.RData to generate a html file, which contains the detection probability of each species at each site.

### Analytic Correlation Files
* analytic_correlations.R contains R code used to implement the analytic correlations in R.

## MNM Models - Bayesian Multi-Species N-Mixture Models for Unmarked Animal Communities

### Simulation Studies Files
mnm.R, mnm_ar.R, mnm_hurdle.R and mnm_hurdle_ar.R contain code used to in simulation studies. These can be used to simulate data for, and fit the MNM model, AR model, hurdle model and hurdle-AR model, respectively.

### Case Study Files
The following files all pertain to the model, fitted to the North American Breeding Bird Survey data, which was deemed best fit by DIC value i.e. the Hurdle-AR model, which contains a lat/long response surface on abundance, and which has detection probability that varies by site and species.

          
* abundance_plots.Rmd can be used in conjunction with N.RData and sites.RData to generate a html file of maps of Alaska, which contain the abundances of each species at each site and year

* detection_probability_plots.Rmd can be used in conjunction with p.RData and sites.RData to generate a html file, which contains the detection probability of each species at each site.



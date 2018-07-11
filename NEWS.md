# Jive 2.0

## Major changes

- This version is a major update of the Jive package. The biggest change is that is uses the 'OUwie.fixed' function of the OUwie package to fit the likelihood of the evolutionary models. This allows more flexibility with the different models (stationary root, sace tree to unit length).
- New models are available: white noise with multiple modes (WNM), Brownian Motion with multiple modes (BMM). 
- All models can now be selected for the evolution of the mean species trait as well.
- Addition of a function to calcuate the marginal likelihood of the models with thermodynamic integration, which are useful for Bayes factors calculations.
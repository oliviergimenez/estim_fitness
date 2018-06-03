Code for "Estimating individual fitness in the wild using capture-recapture data"
==================================================

The paper we wrote is available at https://link.springer.com/article/10.1007/s10144-017-0598-x.

You can get a PDF from [here](https://oliviergimenez.github.io/pubs/Gimenez&Gaillard2017PopEcol.pdf)

We provide here the code to estimate fitness while accounting for imperfect detection of individuals in the wild. We consider lifetime reproductive success (à la Clutton-Brock), individual growth rate (à la McGraw and Caswell) and lifetime individual contribution to population growth (à la Coulson). We illustrate our approach using individual data on a population of roe deer. This is an appendix for our paper "Estimating individual fitness in the wild using capture-recapture data" published in Population Ecology.

This repository contains the following files:

* `deerupdatedentree2y.txt`: the roe deer dataset, with 0 for 'the animal is not seen', 1 for 'the animal is seen without any fawn', 2 for 'the animal is seen with one fawn' and 3 for 'the animal is seen with two fawns'.
* `fit_model.r`: a R script that reads the roe deer dataset and fits a multistate capture-recapture state-space model to these data using JAGS (Table 2 in the paper).
* `calcul_all_fitness_metrics.r`: a R script that calculates all three fitness metrics considered in the paper (Figure 1 in the paper).
* `model_selection_in_esurge.zip`: a zip file that contains the E-SURGE output files of the 40 models fitted to the roe deer dataset (Table 1 in the paper).



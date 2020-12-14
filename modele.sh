#!/bin/bash


### Compute correlations of historical recombination profiles (lambda or rho) of 4 populations of landraces per chr*region with a heterocedastic mixed model
### Outputs are : variance-covariance matrix (*aleatoires*) ; residuals (*ajuste*) ; fixed effets (*fixes*)
## Inputs are : one estimate of lambda (or rho) per population and per interval. 



base=${1}


source ${base}

module purge

module load system/R-3.4.3_bis




Rscript ${r_scripts}modele_lambda.R ${r_intervals}intervals.txt ${r_modele}chrregion_lambda_random.txt ${r_modele}NA.txt ${r_modele}NA.txt

Rscript ${r_scripts}modele_rho.R ${r_intervals}intervals.txt ${r_modele}chrregion_rho_random.txt ${r_modele}NA.txt ${r_modele}NA.txt


rm ${r_modele}NA.txt

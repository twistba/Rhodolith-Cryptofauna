# Rhodolith-Cryptofauna

This repository contains data used in Méndez Trejo et al. Cryptofauna associated with rhodoliths: diversity is species-specific and influenced by habitat.

Data are separated into three spreadsheets each with a corresponding R script, outlined below.


(1)Rhodolith_dimensions.csv --> Rhodolith_dimensions.R

This spreadsheet contains individual rhodolith dimension.

id= Unique identification for each rhodolith collected in this study
site= location at which sample was collected
species= rhodolith species
site_species= combined site species variable
transect= transect ID at which the sample was taken
rhodolith_number= unique rhodolith number for a given transect

buoyant_weight= buoyant weight
volume= volume(cm3)of rhodolith calculated using buoyant weights
volume_free_space= volume fo space between branches calculated using rhodolith volume and volume of sphere using length_x width_y and height_z measurements
tip_number= number of rhodolith branch tips - averaged using 5 times 1cm2 quadrats



(2)cryptofauna_richness.csv --> Cryptofauna_Richness.R

This spreadsheet contains the species richness(number of species present)of cryptofauna for each individual rhodolith.

id= Unique identification for each rhodolith collected in this study
site= location at which sample was collected
species= rhodolith species
site_species= combined site species variable
transect= transect ID at which the sample was taken
rhodolith_number= unique rhodolith number for a given transect

richness= species richness of all cryptofauna
richness_motile= species richness of motile cryptofauna
richness_sessile= species richness of sessile cryptofauna


(3) cryptofauna_species_abundance.csv --> Cryptofauna_Multivariate.R

This spreadsheet contains abundances for all the cryptofaunal species identified on each individual rhodolith

id= Unique identification for each rhodolith collected in this study
site= location at which sample was collected
species= rhodolith species
site_species= combined site species variable
transect= transect ID at which the sample was taken
rhodolith_number= unique rhodolith number for a given transect

Amphithoidae:Trichomusculus= Names of cryptofauna and there assosiated abundances


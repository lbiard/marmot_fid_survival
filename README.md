# marmot_fid_survival

This repository hosts data and R code for Blumstein D. T., Sanchez M., Philson C. S., Bliard L. (2023). Is flight initiation distance associated with longer-term survival in yellow-bellied marmots (*Marmota flaviventer*)? Animal Behaviour. https://doi.org/10.1016/j.anbehav.2023.05.013

Data and code also available on OSF https://doi.org/10.17605/OSF.IO/3WY58 

## GENERAL INFORMATION

1. Title: "Data and code from: Is flight initiation distance associated with longer-term survival in yellow-bellied marmots (Marmota flaviventer)?"

2. Author information:
       
       A.  Name: Dan Blumstein
		   Institution: UCLA
		   Email: marmots@ucla.edu

       B.  Name: McKenna Sanchez
		   Institution: Texas A&M University
		   Email: mckennagsanchez@gmail.com

       C.  Name: Conner Philson
		   Institution: UCLA
		   Email: cphilson@g.ucla.edu

       D.  Name: Louis Bliard
		   Institution: University of Zurich
		   Email: Louis.bliard@evobio.eu


## DATA & FILE OVERVIEW

1. File List: 

- `fid_data_osf.csv` Dataset containing the FID data.
- `summer_data_osf.csv` Dataset containing the summer survival data.
- `winter_data_osf.csv` Dataset containing the winter survival data.


- `bivariate_model_summer.R` Code to look at the among-individual covariance between fid and summer survival.
- `bivariate_model_winter.R` Code to look at the among-individual covariance between fid and winter survival.


## METHODOLOGICAL INFORMATION


1. Methods for processing the data: All data were processed in R.

2. Software-specific information needed to interpret the data:
R statistical software, version 4.1.2. 
Main R packages: cmdstanr, version 0.4.0. brms, version 2.16.3.


### DATA-SPECIFIC INFORMATION FOR: `fid_data_osf.csv`

1. Number of variables: 13

2. Number of rows: 1109

3. Variable List: 

- year: year the fid trial was performed.
- jdate: Julian date of the trial.
- col_area: name of the colony.
- SD: starting distance (in meters).
- AD: alert distance (in meters).
- FID: flight initiation distance (in meters).
- dist_burrow: distance from the nearest burrow  (in meters).
- sex: male (M) or female (F).
- age: age in years.
- agelclass: adult or yearling.
- valley_location: Up valley or Down valley.
- trial_number: number of trial within the season for this individual.
- ID_Code: unique individual identifier.

4. Missing data codes: NA



### DATA-SPECIFIC INFORMATION FOR: `summer_data_osf.csv`

1. Number of variables: 13

2. Number of rows: 1119

3. Variable List: 

- year: year of the observation.
- col_area: name of the colony.
- valley_location: Up valley or Down valley.
- sex: male (M) or female (F).
- summer_survival: whether the individual survived the summer season ("YES"/"NO").
- combined_colony: overall colonies (some colonies in close proximity are merged)
- predator_index: level of predation risk ("low"/"high").
- n_predator_observations: number of predator observartions at the given location.
- total_obslength_min: number of minutes of observation at the given location
- age: age in years.
- agelclass: adult or yearling.
- massjun: imputed mass of the individual in June, in grams.
- ID_Code: unique individual identifier.

4. Missing data codes: NA



### DATA-SPECIFIC INFORMATION FOR: `winter_data_osf.csv`

1. Number of variables: 11

2. Number of rows: 4227

3. Variable List: 

- year: year of the observation.
- Winter_Survival: whether the individual survived the winter season (yes=1/no=0).
- col_area: name of the colony.
- massaug: mass of the individual in August (in grams).
- nb_mass: number of time the mass of the individual was taken during the season.
- age: age in years.
- agelclass: adult or yearling.
- sex: male (M) or female (F).
- valley_location: Up valley or Down valley.
- date_of_melt: Julian day of the year of snowmelt.
- ID_Code: unique individual identifier.

4. Missing data codes: NA




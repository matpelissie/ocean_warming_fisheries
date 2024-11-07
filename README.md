# Ocean warming drives abrupt declines in fish productivity at global scale 

This repository contains the minimal data and R scripts to reproduce the analyses and figures associated with the manuscript entitled "Ocean warming drives abrupt declines in fish productivity at global scale" (submitted to PNAS).

## R

The analyses were run using R version 4.3.3.  
The 'R' directory contains the functions used for performing the analyses and the figures.

## Load R packages

`00_Rpackages.R` loads R packages. Some may not be installed, if so the command to install them is indicated. 

## Download datasets

`01_datasets.R` downloads and prepares datasets.

## Run time series classification

`02_classification.R` defines the fish stocks suitable for the trajectory classification, computes surplus production, and performs the classification for each time series. This step may take about one hour to be completed.

## Perform the analyses and make the figures

`03_analyses.R` produces the figures and additional analyses.

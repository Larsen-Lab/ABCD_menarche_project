# ABCD menarche project

This repository contains the analysis code for "Menarche is an inflection point for mental health and brain development"

## Data Availability

The following datasets are required to run the analysis:

* ABCD Data (Adolescent Brain Cognitive Development Study):
	* Data is available on the ABCD Data Repository through the NIH Brain Development Cohort (NBCD) Data Hub (https://abcdstudy.org).
Data from the latest release (v6.0) was used: https://dx.doi.rog/10.15154/z563-zd24
* ABCD ARMS
	* The ABCD dataframe is split into to two matched groups (group 1 and group 2) to perfom split-half reliability sensitvity analyses. Find more information here: https://osf.io/7xn4f/overview
* ABCD General exposome
	* Socioeconomic status (SES) was accounted for in all analyses using the general exposome.
	The general expsome represents the totality of co-occurring environmental exposures and experiences 
	that an individual encounters throughout the lifespan, in contrast to examining single environmental features 
	(e.g. household income) in isolation.
	Find more information here: https://doi.org/10.1016/j.dcn.2024.101370

## System Requirements

* **R Version:** 4.4.0 (or later)
* **R Studio Version:** 2024-04-24 (or later)
* **Dependencies:** `tidyverse` `mgcv` `gratia` `lme4` `parallel` `ggplot2` `ggpattern`

## Scripts to Plot Results (referred to as Tools)

The following scripts need to be *sourced* in the analysis code in order to plot the results

* `partial_residual_plot_new.R`
* `GAM_posterior_stats.R`
* `GAMM_plotting2_with_2nd_derv_new.R`
* `GAMM_firstderv_table.R`

## Analysis Code

* `ABCD_6.0_brief_communication_2026-04-13.Rmd`


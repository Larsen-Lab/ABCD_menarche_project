# ABCD menarche project

This repository contains the analysis code for "Menarche is an inflection point for mental health and brain development"

## Data availability

The following datasets are required to run the analysis:

* ABCD Data (Adolescent Brain Cognitive Development Study):
** Data is available on the ABCD Data Repository through the NIH Brain Development Cohort (NBCD) Data Hub (https://abcdstudy.org).
Data from the latest release (v6.0) was used: https://dx.doi.rog/10.15154/z563-zd24
* ABCD ARMS
** 
* ABCD General exposome
** 

## System Requirements

* **R Version:** 4.4.0 (or later)
* **R Studio Version:** 2024-04-24 (or later)
* **Dependencies:** `tidyverse` `mgcv` `gratia` `lme4` `parallel` `ggplot2` `ggpattern`

## Scripts that are needed to plot results (referred to as Tools)

* `partial_residual_plot_new.R`
* `GAM_posterior_stats.R`
* `GAMM_plotting2_with_2nd_derv_new.R`
* `GAMM_firstderv_table.R`

## Analysis code

* `ABCD_6.0_brief_communication_2026-04-13.Rmd`


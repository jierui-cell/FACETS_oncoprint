## Overview

**Generates FACETS integer copy number and genomic alterations for MSK-IMPACT patients in oncoprint format simultaneously.** 

This repository contains example code for:
1.	Generating alteration tables (mutations, fusions, copy-number events) from cBioPortal-style files.
2.	Creating arm-level CNA clustering visualizations.
3.	Creating oncoprint plots that display per-gene alterations (e.g., mutations, CNAs, fusions).

The workflow:
1.	Load mutation/fusion/CNA data from local text files (cBioPortal or similar format).
2.	Generate a comprehensive alterations table for each sample and each gene of interest (should take around 4-8 minutes depending on your laptop).
3.	Filter or subset the data to a cohort-of-interest by inputing tumor sample ID matching cbioportal (e.g., by tumor type).
4.	Create arm-level CNA clustering plots and oncoprint plots using a combined function.

Please find `example.ipynb` for details of usage.
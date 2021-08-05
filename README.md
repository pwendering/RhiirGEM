# Reconstruction and validation of a genome-scale metabolic model for _Rhizophagus irregularis_

## Description
This repository contains all code and result files that were generated to analyse the _R. irregularis_ genome-scale metabolic model iRi1574.

## Requirements
* IBM CPLEX solver 12.9
* Gurobi solver 9.1.1 (only myristate uptake analysis)
* MATLAB (tested with MATLAB 2017b and 2020b)
* COBRA toolbox for MATLAB (v3.0) [https://github.com/opencobra/cobratoolbox](https://github.com/opencobra/cobratoolbox)
* R (3.6.3) with additional packages:  wesanderson, plotrix, scales, pheatmap, ape
* Python (3.6)
* all simulations were run under both a Fedora operating system and Windows 10 

## Usage
* make sure software requirements are met
* add all MATLAB code files to the MATLAB search path
* change into top-level directory (RhiirGEM/) or update topDir variable in the analysis scripts

## Repository structure
* RhiirGEM/code/analysis
    + MATLAB code files that were used for the analysis of the iRi1574 model
* RhiirGEM/code/kcats
    + get-kcats.sh calls three python scripts that were used to retrieve information on enzyme turnover numbers and full organism lineages
* RhiirGEM/code/plotting
    + contains one R script per figure/subfigure
* RhiirGEM/data
    + corrected-EC-numbers.csv is a two-column file containing updates from old to new EC standard
    + uniprot.tab contains information on all enzymes of _R. irregularis_ that are present in UniProt
* RhiirGEM/data/kcats
    + kcat-reference-data.tsv contains all turnover numbers retrieved from BRENDA and SABIO-RK
    + kcats_model.txt contains a list of turnover numbers per (irreversible) reaction in the iRi1574 model
* RhiirGEM/data/transcriptomic-data
    + contains averaged expression data of three replicates from the three developmental stages (ERM, IRM, ARB)
    + the data were published by Zeng et al. (2018) [https://doi.org/10.1111/tpj.13908](https://doi.org/10.1111/tpj.13908)
* RhiirGEM/memote-report
    + contains the html-formatted file that was obtained by running the [memote](https://memote.io/) test suite on the iRi1574 model
* RhiirGEM/model
    + contains the iRi1574 model as sbml, xlsx, and mat file
* RhiirGEM/results/carbon-sources
    + results from growth predictions and flux sampling of the iRi1574 model with different media conditions
* RhiirGEM/results/developmental-stages
    + results from growth predictions and flux sampling of the iRi1574 model at three simulated developmental stages
* RhiirGEM/results/figures
    + contains all figures shown in the manuscript
* RhiirGEM/results/fungal-models
    + fungal model files that were used for comparisons (references given in Suppl. Table 2)
* RhiirGEM/results/stats
    + model statistics for summarizing the iRi1574 model

## To reproduce the results:
Run the following MATLAB scripts to reproduce the results of our study
(RhiirGEM/code/analysis/)
* Create files for plotting model statistics:
    + generate\_model\_stat_files.m (several minutes)
* Comparison to published fungal models:
    + compare\_subsystems\_fungal\_models.m (about 30 minutes)
* Growth simulations on single carbon sources:
    + growth\_on\_single\_carbon\_sources.m (about a minute)
* Growth prediction on myristate: 
    + growth\_on\_myristate\_fba.m (seconds)
    + growth\_on\_myristate\_emoment.m (several minutes)
* Simulation of growth on different media:
    + analysis\_carbon\_source\_concentration.m (1-2 hours)
* Growth prediction at developmental stages:
    + analysis\_fungal\_structures.m (about 1 hour)

Output files of these scripts are stored in the respective subdirectories of RhiirGEM/results/.

All figures were generated using R scripts (RhiirGEM/code/plotting/), which are named according to their output figure, which is stored at RhiirGEM/results/figures/.

## Reference
**to be added**

 







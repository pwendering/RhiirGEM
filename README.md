# Reconstruction and validation of a genome-scale metabolic model for _Rhizophagus irregularis_

## Description
This repository contains all code and result files that were generated to analyse the _R. irregularis_ genome-scale metabolic model iRi1572.

## Requirements
* IBM CPLEX solver 12.9
* MATLAB (tested with MATLAB 2017b)
* COBRA toolbox for MATLAB (v3.0) [https://github.com/opencobra/cobratoolbox](https://github.com/opencobra/cobratoolbox)
* R (3.6.3) with additional packages:  wesanderson, plotrix, scales, pheatmap, ape
* Python (3.6)
* all simulations were run under a Fedora operating system

## Usage
* make sure software requirements are met
* add all MATLAB code files to the MATLAB search path
* change into top-level directory (iRi1572/) or update topDir variable in analysis scripts

## Repository structure
* iRi1572/code/analysis
    + MATLAB code files that were used for the analysis of the iRi1572 model
* iRi1572/code/kcats
    + get-kcats.sh calls three python scripts that were used to retrieve information on enzyme turnover numbers and full organism lineages
* iRi1572/code/plotting
    + contains one R script per figure/subfigure
* iRi1572/data
    + corrected-EC-numbers.csv is a two-column file containing updates from old to new EC standard
    + uniprot.tab contains information on all enzymes of _R. irregularis_ that are present in UniProt
* iRi1572/data/kcats
    + kcat-reference-data.tsv contains all turnover numbers retrieved from BRENDA and SABIO-RK
    + kcats_model.txt contains a list of turnover numbers per (irreversible) reaction in the iRi1572 model
* iRi1572/data/transcriptomic-data
    + contains averaged expression data of three replicates from the three developmental stages (ERM, IRM, ARB)
    + the data were published by Zeng et al. (2018) [https://doi.org/10.1111/tpj.13908](https://doi.org/10.1111/tpj.13908)
* iRi1572/memote-report
    + contains the html-formatted file that was obtained by running the [memote](https://memote.io/) test suite on the iRi1572 model
* iRi1572/model
    + contains the iRi1572 model as sbml, xlsx, and mat file
* iRi1572/results/carbon-sources
    + results from growth predictions and flux sampling of the iRi1572 model with different media conditions
* iRi1572/results/developmental-stages
    + results from growth predictions and flux sampling of the iRi1572 model at three simulated developmental stages
* iRi1572/results/figures
    + contains all figures shown in the manuscript
* iRi1572/results/fungal-models
    + 18S rRNA sequence of compared fungal species
    + generated tree description file from 18S rRNA sequences from the [MAFFT online service](https://mafft.cbrc.jp/alignment/server/)
    + E.C. numbers contained in the compared models
    + Jaccard distances between E.C. classes of compared models
* iRi1572/results/stats
    + model statistics for summarizing the iRi1572 model

## To reproduce the results:
Run the following MATLAB scripts to reproduce the results of our study
(iRi1572/code/analysis/)
* Create files for plotting model statistics:
    + generate\_model\_stat_files.m (several minutes)
* Comparison to published fungal models:
    + compare\_to\_fungal\_models.m (several minutes)
* Growth prediction on myristate: 
    + growth\_on\_myristate.m (seconds)
* Simulation of growth on different media:
    + analysis\_carbon\_source\_concentration.m (1-2 hours)
* Growth prediction at developmental stages:
    + analysis\_developmental\_stages.m (about 1 hour)

Output files of these scripts are stored in the respective subdirectories of iRi1572/results/.

All figures were generated using R scripts (iRi1572/code/plotting/), which are named according to their output figure, which is stored at iRi1572/results/figures/.

## Reference
**to be added**

 







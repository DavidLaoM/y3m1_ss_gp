# y3m1_ss_gp

Descrition: Glycolysis revisited: from steady-state growth to glucose pulses.

Citation: TBA.

Keywords: Sacchamoyces cerevisiae; glycolysis; central carbon metabolism; kinetic modeling; parameter estimation; steady state; growing cell; glucose perturbation response 

Description: 
This repository contains the work behind our research and the required code and data to reproduce the figures in the manuscript and working simulations with the model developed. This work consisted of the development of a parameter estimation toolbox suited for our specific case study (data, model complexity), model deveolopment and validation and subsequent analysis of model simulations. 
The model consists of a series of ordinary differential equations (ODE) representing the mass balances involved in each reaction. Reactions covered in yeast (Saccharomyces cerevisiae) glycolysis, trehalose cycle, glycerol branch, and simplified expressions for sink reactions, respiration and fermentation.
The data/conditions studied are steady states (dilution rates between 0.025-0.375 h-1) and a 110 mM glucose perturbation at a dilution rate of 0.1 h-1.

Required software:
A functional matlab installation (MATLAB 9.3 or higher).

Reproducibility:
The entire computational work, including article figures and supplementary figures, can be reproduced. For this we recommend running the files 'F0_...' to 'F8_...' in the root folder. A simple simulation is generated by running the file 'F1_reference_simulations.m'.
An SBML version of this model will be uploaded to the Biomodels and JWS online databases. A running python version will be published in a separate github repo. 

Folder structure: TBA

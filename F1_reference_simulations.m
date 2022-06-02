% % F0_REFERENCE_SIMULATIONS.m
% This code runs and plots the steady state (ss) and glucose perturbation
% (gp) simulations.
% David Lao-Martil, 2022-05-23

% initializing
set_paths;
colorpalette;
[legenda] = legendaFull;
loadData_Y3M1;
Y3M1manus_labels;

% setup simulations and plotting
setup.runSScanelas = 1; % simulate steady states 0.025-0.375 h-1 (Canelas et al., 2011)
setup.runGSvanHeerden = 1; % simulate 110 mM glucose perturbation at 0.1 h-1 (van Heerden et al., 2014)
setup.plotResultsMode = 99; % 99 == plot results. 0 == do not plot.
% setup model structure
setup.biomass = 1; % 1 == include biomass. Else not.
setup.clamp10.Pi = 0; % 0 == cytsolic inorganic phosphate as dynamic variable. 1 == fixed at 10 mM.
setup.clamp10.TRE = 0; % 0 == trehalose cycle as dynamic variable. 1 == fixed at experimental data values. 2 == fixed only at steady state data values
setup.clamp10.NADX = 0; % 0 == mitochondrial NADH recycle. 1 == no mitochondiral recycle + NAD(H) concentrations fixed. 
setup.clamp10.AXP = 0; % 0 == AXP cofactor metabolism modelled. 1 == not modelled + A(X)P concentrations fixed. 
setup.clamp10.IXP = 0; % 0 == inosine salvage pathway modelled. 1 == not modelled + I(X)P concentrations fixed.  
% % % % % setup.plotBackColor = [.94 .94 .94]; % backgroun
% setup parameter set
load('x108c.mat'),
x = x_3rd;

% simulation 
sim_Y3M1 = simulateY3M1_separateGPSS(x, canelas_SS, data, dataset, setup);


% % saving the plots
% savefig(1,'results\F1_SSmets.fig')
% savefig(2,'results\F1_SSrates.fig')
% savefig(3,'results\F1_GPmets.fig')
% savefig(4,'results\F1_GPrates.fig')
% savefig(5,'results\F1_SSfysio.fig')




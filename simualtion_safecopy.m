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
setup.runSScanelas = 1; % SS sims cancelled for now. Only GP to be seen
setup.runGSvanHeerden = 1;
setup.plotResultsMode = 99; %0; %30; %20; %10;
% setup model structure
setup.biomass = 1; % biomass
setup.clamp10.Pi = 0;%1; 
setup.clamp10.TRE = 0;%1; 
setup.clamp10.NADX = 0;%1; 
setup.clamp10.AXP = 0;%1; 
setup.clamp10.IXP = 0; %1; 
setup.plotBackColor = [.94 .94 .94]; % standard grey
% setup parameter set
load('x108c.mat'),
x = x_3rd;

% setup to simplify % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% % % % % setup.conditionsSS = 0;
% % % % % setup.conditionsGP = 0;
% % % % % setup.latestTopology = 1;
% % % % % setup.lastErrorsCorrect = 1;
% % % % % setup.Katpase_fermIncrease = 1;
% % % % % setup.Katpase_fermIncrease2 = 0;
% % % % % setup.startATPgp = 1;
% % % % % setup.importPiON = 1;
% % % % % setup.treCycleON = 1;
% % % % % setup.startParamsOFF = 1;
% % % % % setup.IXPdevelopmentOFF = 0; 
% % % % % setup.adjust_mito_atpase = 1;
% % % % % setup.Amd_IXP = 0;
% % % % % setup.clamp_vALDgp = 0;
% % % % % setup.NADHrecycle_syncPYR = 0;
% % % % % setup.latest_IXP_change = 1;
% % % % % setup.pHchangeON = 0;
% % % % % setup.GPcorrectStart = 1;
% % recallWalther2010data;
% % setup.IXPstartConcs_w2010 = 1;
% % % % % setup.ATPgpissue = 1;
% % % % % setup.separateParsGPSS = 1;
% % % % % setup.final_Vmito_implementation = 1;
% % % % % setup.sameGLTHXK4all = 1;
% % % % % setup.GLT_km_same_extra = 1;
% % % % % setup.VmGLT_mu_dependent = 1; 
% % % % % setup.GLT_km_same = 1;
% % % % % setup.VmGLT_mu_dependent = 1; % HXT growth rate dependency 
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 

% simulation 
sim_Y3M1 = simulateY3M1_separateGPSS(x, canelas_SS, data, dataset, setup);


% % saving the plots
% savefig(1,'results\F1_SSmets.fig')
% savefig(2,'results\F1_SSrates.fig')
% savefig(3,'results\F1_GPmets.fig')
% savefig(4,'results\F1_GPrates.fig')
% savefig(5,'results\F1_SSfysio.fig')




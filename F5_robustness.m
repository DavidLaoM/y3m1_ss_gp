% % F5_ROBUSTNESS.m
% This figure reproduces the simulations with model ensembles.
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
setup.plotResultsMode = 30; % 99 == plot results. 0 == do not plot.
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
x_2nd = x_3rd;

% % % % % % initializing
% % % % % set_paths;
% % % % % colorpalette;
% % % % % [legenda] = legendaFull;
% % % % % loadData_Y3M1;
% % % % % Y3M1manus_labels;
% % % % % 
% % % % % % setup simulations and plotting
% % % % % setup.runSScanelas = 1; % SS sims cancelled for now. Only GP to be seen
% % % % % setup.runGSvanHeerden = 1;
% % % % % setup.plotResultsMode = 30; %0; %30; %20; %10;
% % % % % % setup model structure
% % % % % setup.biomass = 1; % biomass
% % % % % setup.clamp10.Pi = 0;%1; 
% % % % % setup.clamp10.TRE = 0;%1; 
% % % % % setup.clamp10.NADX = 0;%1; 
% % % % % setup.clamp10.AXP = 0;%1; 
% % % % % setup.clamp10.IXP = 0; %1; 
% % % % % setup.plotBackColor = [.94 .94 .94]; % standard grey
% % % % % % setup parameter set
% % % % % load('x108c.mat'),
% % % % % x = x_3rd;
% % % % % x_2nd = x_3rd;
% % % % % 
% % % % % % setup to simplify % % % % % % % % % % % % % % % % % % % % % % % % % % % 
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
% % % % % % % recallWalther2010data;
% % % % % % % setup.IXPstartConcs_w2010 = 1;
% % % % % setup.ATPgpissue = 1;
% % % % % setup.separateParsGPSS = 1;
% % % % % setup.final_Vmito_implementation = 1;
% % % % % setup.sameGLTHXK4all = 1;
% % % % % setup.GLT_km_same_extra = 1;
% % % % % setup.VmGLT_mu_dependent = 1;
% % % % % setup.GLT_km_same = 1;
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 



%% Reference simulation
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
setup.plotResultsMode = 0;
[refSim] = simulateY3M1_separateGPSS(x_2nd, canelas_SS, data, dataset, setup);
setup.plotResultsMode = 0;

% %% plotting
setup.plotResultsMode = 0;
% namesHits = cell(1,2); namesHits{1} = 'start'; namesHits{2} = 'x108res';
namesHits = cell(1,1); namesHits{1} = 'x108res'; %namesHits{2} = 'x108res';
setup.namesHits = namesHits;
% simRes = [simRes_reference, refSim];
simRes = refSim;
% %%
[~] = plotAll_Y3M1(simRes,legenda,canelas_SS,data,setup,x_2nd);
setup.plotResultsMode = 0;


%% Adjusted simulation
setup.GLT_km_same_extra = 1;
setup.VmGLT_mu_dependent = 1; % to re run
% 
x_3rd = x_2nd; x_3rd(35) = x_2nd(38); x_3rd(132) = x_2nd(132) + 0.1;
% 
setup.plotResultsMode = 0;
[adjSim] = simulateY3M1_separateGPSS(x_3rd, canelas_SS, data, dataset, setup);
setup.plotResultsMode = 0;

% %% plotting
setup.plotResultsMode = 30;
namesHits = cell(1,2); namesHits{1} = 'refSim'; namesHits{2} = 'adjSim';
setup.namesHits = namesHits;
simRes = [refSim, adjSim];
% simRes = refSim;
% %%
% [~] = plotAll_Y3M1(simRes,legenda,canelas_SS,data,setup,x_2nd);
[~] = plotAll_Y3M1(simRes,legenda,canelas_SS,data,setup,[x_2nd; x_3rd]);
setup.plotResultsMode = 0;
% % 
% save('x108c.mat','x_3rd')


%% creating random noise array
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
N = 1000; %test multiple combinations
% 10 % noise
rng(1), rand_dist = rand(N,length(x_2nd)) * (log10(1.1)-log10(0.9)) + log10(0.9); % 10% noise
% 20 + (150-20) 
% figure, histogram(rand_dist)

% select only some paraemters to change
parsGLT = 35:36;
parsHXK = 28:34;
parsPFK = 43:56;
parsPGI = 57:60;
parsFBA = 11:15;

parsTPI = 79:82;
parsPGM = 67:70;
parsENO = 16:19;
parsPYK = 71:77;
parsTDH = 21:27;

parsPGK = 61:66;
parsGPD = 96:104;
parsPDC = 39:42;
parsADH = 1:10;

% additions rest
parsHOR2 = 105:107;

% addition tre cycle
parsPGM1 = 83:86;
parsTPS1 = 124:128;
parsTPS2 = 119:121;
parsNTH1 = 122:123;
parsTRE = [parsPGM1, parsTPS1, parsTPS2, parsNTH1];

% addition transporters
parsTransporters = [95, 108];

% additions cofactors
parsMitoNADH = 92:93;
parsMitoATP = 109:111;
parsATPase = 129;
parsADK1 = 130:131;
parsVac = 37;
parsCofactors = [parsMitoNADH, parsMitoATP, parsATPase, parsADK1, parsVac];

% addition inosine
parsInosine = 112:118;

% additions cofactors, growth rate-dependent
parsMitoNADH_mu = 132;
parsMitoATP_mu = 133;
parsATPase_mu = 134;
parsADK1_mu = 135:136;
parsCofactors_mu = [parsMitoNADH_mu, parsMitoATP_mu, parsATPase_mu, parsADK1_mu];

% creating selPars
selPars = [parsGLT, parsHXK, parsPFK, parsPGI, parsFBA,...
    parsTPI, parsPGM, parsENO, parsPYK, parsTDH,...
    parsPGK, parsGPD, parsPDC, parsADH];
% selPars_HOR2 = [selPars, parsHOR2];
% selPars_TRE = [selPars, parsTRE];
% selPars_Transporters = [selPars, parsTransporters];
% selPars_Cofactors = [selPars, parsCofactors];
% selPars_Inosine = [selPars, parsInosine];
% selPars_Cofactors_mu = [selPars, parsCofactors_mu];
% selPars_all = selPars_Cofactors_mu;
selPars = [selPars, parsHOR2, parsTRE, parsTransporters, parsCofactors, parsInosine, parsCofactors_mu];

% 
xRobustness = zeros(N,length(x_2nd));
for i = 1:N
%     xRobustness(i,:) = x_2nd + rand_dist(i,:);
    xRobustness(i,:) = x_2nd;
    xRobustness(i,selPars) = x_2nd(selPars) + rand_dist(i,selPars);
end


%% Other simulations in the sensitivity analysis
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 

% % % loop to simulate the robustness check
% % simRobust_r1 = cell(1,250);
% % cluster = parcluster('local');
% % pool = parpool(cluster,16);
% % % pool = parpool(cluster,4);
% % parfor i = 1:250
% %     disp(i)
% %     [robSim] = simulateY3M1_separateGPSS(xRobustness(i,:), canelas_SS, data, dataset, setup);
% %     % ss
% %     temp_Y = zeros(8,32);
% %         temp_Y(1,:) = robSim{1}.ss.Y.ss2can{1}(end,:);
% %         temp_Y(2,:) = robSim{1}.ss.Y.ss2can{2}(end,:);
% %         temp_Y(3,:) = robSim{1}.ss.Y.ss2can{3}(end,:);
% %         temp_Y(4,:) = robSim{1}.ss.Y.ss2can{4}(end,:);
% %         temp_Y(5,:) = robSim{1}.ss.Y.ss2can{5}(end,:);
% %         temp_Y(6,:) = robSim{1}.ss.Y.ss2can{6}(end,:);
% %         temp_Y(7,:) = robSim{1}.ss.Y.ss2can{7}(end,:);
% %         temp_Y(8,:) = robSim{1}.ss.Y.ss2can{8}(end,:);
% %         simRobust_r1{i}.ss.Y = temp_Y;
% %     temp_V = zeros(8,42);
% %         temp_V(1,:) = robSim{1}.ss.V.ss2can{1}(end,:);
% %         temp_V(2,:) = robSim{1}.ss.V.ss2can{2}(end,:);
% %         temp_V(3,:) = robSim{1}.ss.V.ss2can{3}(end,:);
% %         temp_V(4,:) = robSim{1}.ss.V.ss2can{4}(end,:);
% %         temp_V(5,:) = robSim{1}.ss.V.ss2can{5}(end,:);
% %         temp_V(6,:) = robSim{1}.ss.V.ss2can{6}(end,:);
% %         temp_V(7,:) = robSim{1}.ss.V.ss2can{7}(end,:);
% %         temp_V(8,:) = robSim{1}.ss.V.ss2can{8}(end,:);
% %         simRobust_r1{i}.ss.V = temp_V;
% %     % gs
% %     simRobust_r1{i}.gs.T = robSim{1}.gs.T(3002:end); 
% %     simRobust_r1{i}.gs.Y = robSim{1}.gs.Y(3002:end,:); 
% %     simRobust_r1{i}.gs.V = robSim{1}.gs.V(3002:end,:); 
% % end
% % % save('simRes_Y3M1_mF5_robust_core_enzymes_round1_cluster.mat','refSim','simRobust_r1')
% % % save('simRes_Y3M1_mF5_robust_core_enzymes_round1_local.mat','refSim','simRobust_r1')
% % save('simRes_Y3M1_mF5_robust_all_enzymes_round1.mat','refSim','simRobust_r1')
% % delete(pool)
% % 
% % % loop to simulate the robustness check
% % simRobust_r2 = cell(1,250);
% % cluster = parcluster('local');
% % pool = parpool(cluster,16);
% % % pool = parpool(cluster,4);
% % parfor i = 251:500
% %     disp(i)
% %     [robSim] = simulateY3M1_separateGPSS(xRobustness(i,:), canelas_SS, data, dataset, setup);
% %     % ss
% %     temp_Y = zeros(8,32);
% %         temp_Y(1,:) = robSim{1}.ss.Y.ss2can{1}(end,:);
% %         temp_Y(2,:) = robSim{1}.ss.Y.ss2can{2}(end,:);
% %         temp_Y(3,:) = robSim{1}.ss.Y.ss2can{3}(end,:);
% %         temp_Y(4,:) = robSim{1}.ss.Y.ss2can{4}(end,:);
% %         temp_Y(5,:) = robSim{1}.ss.Y.ss2can{5}(end,:);
% %         temp_Y(6,:) = robSim{1}.ss.Y.ss2can{6}(end,:);
% %         temp_Y(7,:) = robSim{1}.ss.Y.ss2can{7}(end,:);
% %         temp_Y(8,:) = robSim{1}.ss.Y.ss2can{8}(end,:);
% %         simRobust_r2{i-250}.ss.Y = temp_Y;
% %     temp_V = zeros(8,42);
% %         temp_V(1,:) = robSim{1}.ss.V.ss2can{1}(end,:);
% %         temp_V(2,:) = robSim{1}.ss.V.ss2can{2}(end,:);
% %         temp_V(3,:) = robSim{1}.ss.V.ss2can{3}(end,:);
% %         temp_V(4,:) = robSim{1}.ss.V.ss2can{4}(end,:);
% %         temp_V(5,:) = robSim{1}.ss.V.ss2can{5}(end,:);
% %         temp_V(6,:) = robSim{1}.ss.V.ss2can{6}(end,:);
% %         temp_V(7,:) = robSim{1}.ss.V.ss2can{7}(end,:);
% %         temp_V(8,:) = robSim{1}.ss.V.ss2can{8}(end,:);
% %         simRobust_r2{i-250}.ss.V = temp_V;
% %     % gs
% %     simRobust_r2{i-250}.gs.T = robSim{1}.gs.T(3002:end); 
% %     simRobust_r2{i-250}.gs.Y = robSim{1}.gs.Y(3002:end,:); 
% %     simRobust_r2{i-250}.gs.V = robSim{1}.gs.V(3002:end,:); 
% % end
% % % save('simRes_Y3M1_mF5_robust_core_enzymes_round2_cluster.mat','refSim','simRobust_r2')
% % % save('simRes_Y3M1_mF5_robust_core_enzymes_round2_local.mat','refSim','simRobust_r2')
% % save('simRes_Y3M1_mF5_robust_all_enzymes_round2.mat','simRobust_r2')
% % delete(pool)
% % 
% % % loop to simulate the robustness check
% % simRobust_r3 = cell(1,250);
% % cluster = parcluster('local');
% % pool = parpool(cluster,16);
% % % pool = parpool(cluster,4);
% % parfor i = 501:750
% %     disp(i)
% %     [robSim] = simulateY3M1_separateGPSS(xRobustness(i,:), canelas_SS, data, dataset, setup);
% %     % ss
% %     temp_Y = zeros(8,32);
% %         temp_Y(1,:) = robSim{1}.ss.Y.ss2can{1}(end,:);
% %         temp_Y(2,:) = robSim{1}.ss.Y.ss2can{2}(end,:);
% %         temp_Y(3,:) = robSim{1}.ss.Y.ss2can{3}(end,:);
% %         temp_Y(4,:) = robSim{1}.ss.Y.ss2can{4}(end,:);
% %         temp_Y(5,:) = robSim{1}.ss.Y.ss2can{5}(end,:);
% %         temp_Y(6,:) = robSim{1}.ss.Y.ss2can{6}(end,:);
% %         temp_Y(7,:) = robSim{1}.ss.Y.ss2can{7}(end,:);
% %         temp_Y(8,:) = robSim{1}.ss.Y.ss2can{8}(end,:);
% %         simRobust_r3{i-500}.ss.Y = temp_Y;
% %     temp_V = zeros(8,42);
% %         temp_V(1,:) = robSim{1}.ss.V.ss2can{1}(end,:);
% %         temp_V(2,:) = robSim{1}.ss.V.ss2can{2}(end,:);
% %         temp_V(3,:) = robSim{1}.ss.V.ss2can{3}(end,:);
% %         temp_V(4,:) = robSim{1}.ss.V.ss2can{4}(end,:);
% %         temp_V(5,:) = robSim{1}.ss.V.ss2can{5}(end,:);
% %         temp_V(6,:) = robSim{1}.ss.V.ss2can{6}(end,:);
% %         temp_V(7,:) = robSim{1}.ss.V.ss2can{7}(end,:);
% %         temp_V(8,:) = robSim{1}.ss.V.ss2can{8}(end,:);
% %         simRobust_r3{i-500}.ss.V = temp_V;
% %     % gs
% %     simRobust_r3{i-500}.gs.T = robSim{1}.gs.T(3002:end); 
% %     simRobust_r3{i-500}.gs.Y = robSim{1}.gs.Y(3002:end,:); 
% %     simRobust_r3{i-500}.gs.V = robSim{1}.gs.V(3002:end,:); 
% % end
% % % save('simRes_Y3M1_mF5_robust_round3.mat','refSim','simRobust_r3')
% % save('simRes_Y3M1_mF5_robust_all_enzymes_round3.mat','simRobust_r3')
% % delete(pool)
% % 
% % % loop to simulate the robustness check
% % simRobust_r4 = cell(1,250);
% % cluster = parcluster('local');
% % pool = parpool(cluster,16);
% % % pool = parpool(cluster,4);
% % parfor i = 751:1000
% %     disp(i)
% %     [robSim] = simulateY3M1_separateGPSS(xRobustness(i,:), canelas_SS, data, dataset, setup);
% %     % ss
% %     temp_Y = zeros(8,32);
% %         temp_Y(1,:) = robSim{1}.ss.Y.ss2can{1}(end,:);
% %         temp_Y(2,:) = robSim{1}.ss.Y.ss2can{2}(end,:);
% %         temp_Y(3,:) = robSim{1}.ss.Y.ss2can{3}(end,:);
% %         temp_Y(4,:) = robSim{1}.ss.Y.ss2can{4}(end,:);
% %         temp_Y(5,:) = robSim{1}.ss.Y.ss2can{5}(end,:);
% %         temp_Y(6,:) = robSim{1}.ss.Y.ss2can{6}(end,:);
% %         temp_Y(7,:) = robSim{1}.ss.Y.ss2can{7}(end,:);
% %         temp_Y(8,:) = robSim{1}.ss.Y.ss2can{8}(end,:);
% %         simRobust_r4{i-750}.ss.Y = temp_Y;
% %     temp_V = zeros(8,42);
% %         temp_V(1,:) = robSim{1}.ss.V.ss2can{1}(end,:);
% %         temp_V(2,:) = robSim{1}.ss.V.ss2can{2}(end,:);
% %         temp_V(3,:) = robSim{1}.ss.V.ss2can{3}(end,:);
% %         temp_V(4,:) = robSim{1}.ss.V.ss2can{4}(end,:);
% %         temp_V(5,:) = robSim{1}.ss.V.ss2can{5}(end,:);
% %         temp_V(6,:) = robSim{1}.ss.V.ss2can{6}(end,:);
% %         temp_V(7,:) = robSim{1}.ss.V.ss2can{7}(end,:);
% %         temp_V(8,:) = robSim{1}.ss.V.ss2can{8}(end,:);
% %         simRobust_r4{i-750}.ss.V = temp_V;
% %     % gs
% %     simRobust_r4{i-750}.gs.T = robSim{1}.gs.T(3002:end); 
% %     simRobust_r4{i-750}.gs.Y = robSim{1}.gs.Y(3002:end,:); 
% %     simRobust_r4{i-750}.gs.V = robSim{1}.gs.V(3002:end,:); 
% % end
% % % save('simRes_Y3M1_mF5_robust_round4.mat','refSim','simRobust_r4')
% % save('simRes_Y3M1_mF5_robust_all_enzymes_round4.mat','simRobust_r4')
% % delete(pool)


%% load and better split results array
% % 
% load('simRes_Y3M1_mF5_robust_round1.mat','refSim','simRobust_r1')
% load('simRes_Y3M1_mF5_robust_round2.mat','simRobust_r2')
% load('simRes_Y3M1_mF5_robust_round3.mat','simRobust_r3')
% load('simRes_Y3M1_mF5_robust_round4.mat','simRobust_r4')
% simRobust = [simRobust_r1, simRobust_r2, simRobust_r3, simRobust_r4];
% clear simRobust_r1 simRobust_r2 simRobust_r3 simRobust_r4

% % check, local
% load('simRes_Y3M1_mF5_robust_core_enzymes_round1_local.mat','refSim','simRobust_r1')
% load('simRes_Y3M1_mF5_robust_core_enzymes_round2_local.mat','simRobust_r2')
% simRobust = [simRobust_r1, simRobust_r2];
% clear simRobust_r1 simRobust_r2

% % check, cluster
% load('simRes_Y3M1_mF5_robust_core_enzymes_round1_cluster.mat','refSim','simRobust_r1')
% load('simRes_Y3M1_mF5_robust_core_enzymes_round2_cluster.mat','simRobust_r2')
% simRobust = [simRobust_r1, simRobust_r2];
% clear simRobust_r1 simRobust_r2

% 
% load('safecopy_simRes_Y3M1_mF5_robust_round1.mat','refSim','simRobust_r1')
% load('safecopy_simRes_Y3M1_mF5_robust_round2.mat','simRobust_r2')
% load('safecopy_simRes_Y3M1_mF5_robust_round3.mat','simRobust_r3')
% load('safecopy_simRes_Y3M1_mF5_robust_round4.mat','simRobust_r4')
% simRobust = [simRobust_r1, simRobust_r2, simRobust_r3, simRobust_r4];
% clear simRobust_r1 simRobust_r2 simRobust_r3 simRobust_r4

% 
load('simRes_Y3M1_mF5_robust_all_enzymes_round1.mat','refSim','simRobust_r1')
load('simRes_Y3M1_mF5_robust_all_enzymes_round2.mat','simRobust_r2')
load('simRes_Y3M1_mF5_robust_all_enzymes_round3.mat','simRobust_r3')
load('simRes_Y3M1_mF5_robust_all_enzymes_round4.mat','simRobust_r4')
simRobust = [simRobust_r1, simRobust_r2, simRobust_r3, simRobust_r4];
clear simRobust_r1 simRobust_r2 simRobust_r3 simRobust_r4


%% Plotting
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 

% % dummy plotting
% simRobust_safecopy = simRobust;
% ntemp = 1000;
% simRobust_dummy = cell(1,ntemp);
% for i = 1:ntemp
%     temp.gs.T = refSim{1}.gs.T(3002:end,:);
%     temp.gs.Y = refSim{1}.gs.Y(3002:end,:);
%     temp.gs.V = refSim{1}.gs.V(3002:end,:);
%     simRobust_dummy{i} = temp;
% end
% simRobust = simRobust_dummy;
% % 
% plot_robustnessSimulations_Y3M1
% % 
% simRobust = simRobust_safecopy;
plot_robustnessSimulations_Y3M1


% % saving the plots
% savefig(2001,'results\F5_robust_ss_mets.fig')
% savefig(2002,'results\F5_robust_ss_rates.fig')
% savefig(2003,'results\F5_robust_gp_mets.fig')
% savefig(2004,'results\F5_robust_gp_rates.fig')






%% Starting setup
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% % % % % 
% % % % clear, %close all
% % % % set_paths;
% % % % path(path, strcat(folder,'/tempResults'))
% % % % dbstop if error
% % % % colorpalette;
% % % % for starting_setup_1 = 1
% % % %     % get the starting parameters, datasets and legend
% % % %     [legenda] = legendaFull;
% % % %     loadData_Y3M1;
% % % %     % simulation setup
% % % %     setup.runSScanelas = 1; % SS sims cancelled for now. Only GP to be seen
% % % %     setup.runGSvanHeerden = 1;
% % % %     setup.plotResultsMode = 0; %30; %20; %10;
% % % %     setup.conditionsSS = 0;
% % % %     setup.conditionsGP = 0;
% % % %     % model structure
% % % %     setup.biomass = 1; % biomass
% % % %     setup.clamp10.Pi = 0;%1; 
% % % %     setup.clamp10.TRE = 0;%1; 
% % % %     setup.clamp10.NADX = 0;%1; 
% % % %     setup.clamp10.AXP = 0;%1; 
% % % %     setup.clamp10.IXP = 0; %1; 
% % % %     % setup updated to last topology
% % % %     setup.latestTopology = 1;
% % % %     setup.lastErrorsCorrect = 1;
% % % %     % prior parameter adjustments
% % % %     for prior = 1
% % % %         load('x_step3_3TREreg.mat','x_step3_3TREreg'); pset2 = x_step3_3TREreg;
% % % %         load('x_step3_1highMU_selected.mat'); pset2(137:146) = x_highMU;
% % % %         pars2zeros = [...%38 144,...
% % % %             28 29 30 31 139 140,...
% % % %             46 49 50 51 87,...
% % % %             12,14,...
% % % %             80,81,...
% % % %             71,76,...
% % % %             ];
% % % %         pset2(pars2zeros) = zeros;
% % % %         x = pset2;
% % % %             % 2021-02-22 changed: lowering the intensity of glucose entry via GLT
% % % %             x(35) = -0.6250;
% % % %             % 2021-02-25 changed: increased trehalose cycle uptake (from second test)
% % % %             x([83    86   124   125   126   127   128   119   120    33]) = [...
% % % %                 1.7707   -0.8062   -0.3905   -0.7950    0.7822   -0.2262   -0.5143    0.1574   -0.3881   -1.5375,...
% % % %                 ]; % option B
% % % %         %     x([83    86   124   125   126   127   128   119   120    33]) = [...
% % % %         %         1.8444    1.9588   -0.3952   -0.8152    0.8085   -0.3085   -0.5936    0.1164   -0.3871   -1.5606,...
% % % %         %         ]; % option C
% % % %             % 2021-02-25 changed: adjusted nth1 kcat
% % % %             x(123) = -0.9528 - 0.6;
% % % %             x1 = x;
% % % %             % 2021-03-04 changed: adjustment on pgi-pfk-fba-tpi
% % % %             idxs_change =   [11    13    15    43    44,...    
% % % %                 45    47    48    52    53,...
% % % %                 54    55    56    57    58,...    
% % % %                 59    60    79    82];
% % % %             x(idxs_change) = [...
% % % %                 -2.9923    0.4644   -0.0912   -1.3570   -0.0428,...
% % % %                 0.9930   -0.0624   -1.3087    1.2157    0.9114,...
% % % %                 0.0475   -0.4054    1.8658    0.4785    1.3032,...
% % % %                 1.5172   -0.7006    0.0404   -0.1552];
% % % %             x2 = x;
% % % %             % 2021-03-09 changed: IXP influence lowered for correct ATP start +
% % % %             % ATPase
% % % %                 % ATPase possibility for effect
% % % %                 setup.Katpase_fermIncrease = 1;
% % % %                 setup.Katpase_fermIncrease2 = 0;
% % % %                 % ATPase parameter number 
% % % %                 % -> x109
% % % %                 % Start ATP
% % % %                 setup.startATPgp = 1;
% % % %                 % lower intake by IXP
% % % %                 x(113) = .9;
% % % %             % 2021-03-11 changed: ATP profile fitted by pEst (78b, line720)
% % % %             x([43 44 47 48 49 53, 130 131, 113]) = [...
% % % %                 -1.1530   -0.0000   -0.1086   -1.1062,...
% % % %                 0.1840    0.6468   -0.3382    0.0761, 1.4793];
% % % %             % 2021-03-12 changed: start glk to fit the SS rates (pE79a, line820)
% % % %             x([11 12 13 14 15 79 80 81 82 140]) = [...
% % % %                 -2.4033    0.2485    0.1829    0.2502    0.1458...
% % % %                  1.3006    0.1197    0.2267   -0.6856    1.0000];
% % % %                 % x140 could also be set to zero, but then the rates should be
% % % %                 % better fit in another way.
% % % %     end
% % % %     load('pset_SS_GP_20210406.mat','x_SSsel','x_GPsel');
% % % %     load('parameters_step1.mat');
% % % % 
% % % %     % visualization color
% % % %     % setup.plotBackColor = [.5 .5 1]; % blue
% % % %     % setup.plotBackColor = [.5 1 .5]; % green
% % % %     % setup.plotBackColor = [1 .5 .5]; % red
% % % %     setup.plotBackColor = [.94 .94 .94]; % standard grey
% % % %     % latest changes (2020 - 08 - 13)
% % % %     setup.importPiON = 1;
% % % %     setup.treCycleON = 1;
% % % %     setup.startParamsOFF = 1;
% % % %     setup.IXPdevelopmentOFF = 0; 
% % % %     setup.adjust_mito_atpase = 1;
% % % %     setup.Amd_IXP = 0;
% % % %     setup.clamp_vALDgp = 0;
% % % %     setup.NADHrecycle_syncPYR = 0;
% % % %     setup.latest_IXP_change = 1;
% % % %     setup.pHchangeON = 0;
% % % %     % latest changes (2020 - 11 - 20)
% % % %     setup.GPcorrectStart = 1;
% % % %     % latest changes (2021 - 04 - 16)
% % % %     recallWalther2010data;
% % % %     setup.IXPstartConcs_w2010 = 1;
% % % %     setup.ATPgpissue = 1;
% % % %     % save('pset_SS_GP_20210406.mat','x_SSsel','x_GPsel');
% % % %     load('pset_SS_GP_20210406.mat','x_SSsel','x_GPsel');
% % % %     % saved variables: x_SSsel, x_GPsel
% % % %     x = x_SSsel;
% % % %     %     x([149 148 150 151]) = [1.0005    1.4666    1.5084    1.5043]; % 59    52    79    22
% % % %     %     x(147) = x(11);
% % % %         x(147:151) = [-1.9033    1.4666    1.5084    1.5043    1.0005];
% % % %         % xAll1_pC28(4,[11 52 59 79 22])
% % % %         % xAll1_pC28(4,[147 148 149 150 151])
% % % %         % 
% % % %         % ans =
% % % %         % 
% % % %         %    -1.9033    1.4666    1.5084    1.5043    1.0005
% % % % 
% % % %     %     x(147) = x(11);     % fba  (11:15)
% % % %     %     x(148) = x(52);     % pfk  (43:56)
% % % %     %     x(149) = x(59);     % pgi  (57:60)
% % % %     %     x(150) = x(79);     % tpi  (79:82)
% % % %     %     x(151) = x(22);     % gapdh(21:27)
% % % %     %     x(147) = x_GPsel(11);     % fba  (11:15)
% % % %     %     x(148) = x_GPsel(52);     % pfk  (43:56)
% % % %     %     x(149) = x_GPsel(59);     % pgi  (57:60)
% % % %     %     x(150) = x_GPsel(79);     % tpi  (79:82)
% % % %     %     x(151) = x_GPsel(22);     % gapdh(21:27)
% % % %     setup.separateParsGPSS = 1;
% % % %     x(57:60) = [0.5680    1.3728    0.8920   -0.7635];
% % % %     % (2021 05 11)
% % % %     x(83) = 2.0549;
% % % %     % (2021 05 18) latest change on glt to fit GP concentrations better
% % % %     load('tempParameterSave_20210512.mat')
% % % %     x = x97b_up;
% % % %     % (2021 05 20) 
% % % %     setup.final_Vmito_implementation = 1;
% % % %     % (2021 05 24)
% % % %     x(109:111) = [0.0042    0.0252   -0.0762];
% % % %     % x(152) = 0; x2 = x; x2(152) = 100;% IXP is off before pulse
% % % %     x(152) = 100; % IXP is off before pulse
% % % % 
% % % %     % (2021 09? 27) finally looking good:
% % % %     load('parSet_104a.mat','xFinal');
% % % %     x = xFinal;
% % % %     x(71:77) = zeros;
% % % % 
% % % %     % (2021 10 25) regularization UG
% % % %     load('parSet_105b.mat','x105b');
% % % %     x([43:56, 148]) = x105b([43:56, 148]); % parsPFK
% % % %     x(73) = x105b(73); % parsPYK (kmfbp)
% % % % 
% % % %     % 
% % % %     Y3M1manus_labels
% % % %     % %% Reference simulation
% % % %     setup.GLT_km_same = 1;
% % % %     % %% (1) Recall x108a estimate + simple simulation
% % % %     setup.sameGLTHXK4all = 1;
% % % %     
% % % % end
% % % % % load parameters
% % % % load('x108a.mat','x_1st','x_2nd')



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
setup.conditionsSS = 0;
setup.conditionsGP = 0;
setup.latestTopology = 1;
setup.lastErrorsCorrect = 1;
setup.Katpase_fermIncrease = 1;
setup.Katpase_fermIncrease2 = 0;
setup.startATPgp = 1;
setup.importPiON = 1;
setup.treCycleON = 1;
setup.startParamsOFF = 1;
setup.IXPdevelopmentOFF = 0; 
setup.adjust_mito_atpase = 1;
setup.Amd_IXP = 0;
setup.clamp_vALDgp = 0;
setup.NADHrecycle_syncPYR = 0;
setup.latest_IXP_change = 1;
setup.pHchangeON = 0;
setup.GPcorrectStart = 1;
% % recallWalther2010data;
% % setup.IXPstartConcs_w2010 = 1;
setup.ATPgpissue = 1;
setup.separateParsGPSS = 1;
setup.final_Vmito_implementation = 1;
setup.sameGLTHXK4all = 1;
setup.GLT_km_same_extra = 1;
setup.VmGLT_mu_dependent = 1;
setup.GLT_km_same = 1;
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 

setup.save = 0;



x_2nd = x_3rd; 
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


%% making final figure
% disp('making final figure')


%% recall plot to put in the background
% %%
I = imread('fig1_model_diagram.jpg');
f6001_h = figure('Name','fig6_robust_SSconc');
imshow(I)
% % f6002_h = figure('Name','fig6_robust_SSrate');
% % imshow(I)
% % f6003_h = figure('Name','fig6_robust_GPconc');
% % imshow(I)
% % f6004_h = figure('Name','fig6_robust_GPrate');
% % imshow(I)
% 
f6001_h.Position = [1921 257 1536 747];
% % f6002_h.Position = [1921 257 1536 747];
% % f6003_h.Position = [1921 257 1536 747];
% % f6004_h.Position = [1921 257 1536 747];

% f1_h -> fig_h_c_ss
    % f1_h.Children(end+1-i)
% f2_h -> fig_h_r_ss
% f3_h -> fig_h_c_gs
% f4_h -> fig_h_r_gs

%%
for i = 3
    new_f1_h = copyobj(fig_h_c_ss.Children(end+1-i),f6001_h);
%     [0.1300 0.7441 0.0759 0.1512]
end





%% 
% %% 
new_f1_h = copyobj([f1_h.Children(1),f1_h.Children(2)],3);
% %% 
new_f1_h(1).Position = [0.5095 0.4099 0.0854 0.0934]; %[0.5068 0.7220 0.0854 0.0934];
new_f1_h(2).Position = [0.3890 0.3934 0.2183 0.2437]; %[0.4154 0.7176 0.1811 0.2162];


%% in the background
% recall
% correct label the y-axis (probably move it a bit)
% correct label the x-axis (if needed)
% adjust colormap according to data type

% % % % % % % % % % % % % % % % 
% % SS concentrations (f2001_h)
% G6P
% F6P
% FBP
% GAP
% DHAP
% GLYCEROL3P

% P3G
% P2G
% PEP
% PYR

% GLYCEROL
% ETOH

% ATP
% ADP
% AMP
% or Energy charge

% % % % % % % % % % % % % % % % 
% % SS rates (f2002_h)

% vHXT

% vPYK

% vPGM1
% vTPS1
% vTPS2
% vNTH1


% % % % % % % % % % % % % % % % 
% % GP concentrations (f2003_h)
% G6P
% F6P
% FBP
% GAP
% DHAP
% GLYCEROL3P

% P3G
% P2G
% PEP
% PYR

% G1P
% T6P
% TRE

% GLYCEROL
% ETOH

% ATP
% ADP
% AMP
% or Energy charge

% % % % % % % % % % % % % % % % 
% % GP rates (f20042_h)

% vHXT
% vHXK
% vPGI
% vPFK
% vALD

% vPYK

% vPGM1
% vTPS1
% vTPS2
% vNTH1

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% %% saving 2001 SSconc
% savefig_loc = 'E:\OneDrive - TU Eindhoven\Documents\ch4_Y3M1\results\figures\manuscript\';
% saveName = [savefig_loc, 'fig6sup_SSconc']; savefig(2001, saveName);
% figure(2001)
% 
% % savefig(1, 'tempfig21.fig')
% % specs printing (method 3)
% set(gcf,'Units','inches');
% screenposition = get(gcf,'Position');
% set(gcf,...
%     'PaperPosition',[0 0 screenposition(3:4)],...
%     'PaperSize',[screenposition(3:4)]);
% print -dpdf -painters fig6sup_SSconc
% print -dpng -painters fig6sup_SSconc
% 
% 
% % %% saving 2002 SSrate
% savefig_loc = 'E:\OneDrive - TU Eindhoven\Documents\ch4_Y3M1\results\figures\manuscript\';
% saveName = [savefig_loc, 'fig6sup_SSrate']; savefig(2002, saveName);
% figure(2002)
% 
% % savefig(1, 'tempfig21.fig')
% % specs printing (method 3)
% set(gcf,'Units','inches');
% screenposition = get(gcf,'Position');
% set(gcf,...
%     'PaperPosition',[0 0 screenposition(3:4)],...
%     'PaperSize',[screenposition(3:4)]);
% print -dpdf -painters fig6sup_SSrate
% print -dpng -painters fig6sup_SSrate
% 
% 
% % %% saving 2003 GPconc
% savefig_loc = 'E:\OneDrive - TU Eindhoven\Documents\ch4_Y3M1\results\figures\manuscript\';
% saveName = [savefig_loc, 'fig6sup_GPconc']; savefig(2003, saveName);
% figure(2003)
% 
% % savefig(1, 'tempfig21.fig')
% % specs printing (method 3)
% set(gcf,'Units','inches');
% screenposition = get(gcf,'Position');
% set(gcf,...
%     'PaperPosition',[0 0 screenposition(3:4)],...
%     'PaperSize',[screenposition(3:4)]);
% print -dpdf -painters fig6sup_GPconc
% print -dpng -painters fig6sup_GPconc
% 
% 
% % %% saving 2004 GPrate
% savefig_loc = 'E:\OneDrive - TU Eindhoven\Documents\ch4_Y3M1\results\figures\manuscript\';
% saveName = [savefig_loc, 'fig6sup_GPrate']; savefig(2004, saveName);
% figure(2004)
% 
% % savefig(1, 'tempfig21.fig')
% % specs printing (method 3)
% set(gcf,'Units','inches');
% screenposition = get(gcf,'Position');
% set(gcf,...
%     'PaperPosition',[0 0 screenposition(3:4)],...
%     'PaperSize',[screenposition(3:4)]);
% print -dpdf -painters fig6sup_GPrate
% print -dpng -painters fig6sup_GPrate



% % %% examine why the lack of NADH
% % % check when it goes down
% % reacts_to_plot = 1;
% % % 
% % figure(1102)
% % for i = 1:1000
% % % for i = 1:500
% %     plot(canelas_SS.mtlD.D, simRobust{i}.ss.V(:,reacts_to_plot),...
% %         'linewidth',1,'color',[.75 .75 .75])
% %     if i == 1
% %         hold on
% %     end
% % end
% % temp_plot = zeros(1,8);
% %     temp_plot(1) = refSim{1}.ss.V.ss2can{1}(end,reacts_to_plot);
% %     temp_plot(2) = refSim{1}.ss.V.ss2can{2}(end,reacts_to_plot);
% %     temp_plot(3) = refSim{1}.ss.V.ss2can{3}(end,reacts_to_plot);
% %     temp_plot(4) = refSim{1}.ss.V.ss2can{4}(end,reacts_to_plot);
% %     temp_plot(5) = refSim{1}.ss.V.ss2can{5}(end,reacts_to_plot);
% %     temp_plot(6) = refSim{1}.ss.V.ss2can{6}(end,reacts_to_plot);
% %     temp_plot(7) = refSim{1}.ss.V.ss2can{7}(end,reacts_to_plot);
% %     temp_plot(8) = refSim{1}.ss.V.ss2can{8}(end,reacts_to_plot);
% % plot(canelas_SS.mtlD.D, temp_plot,...
% %     'linewidth',1,'color','black')
% % hold off
% % %% How many fell down there? select them
% % % only 27 cases out of 1000, 2.7%.
% % idxs_lessNADHmito = [];
% % for i = 1:1000
% %     % 
% %     val2scan = simRobust{i}.ss.V(8,reacts_to_plot);
% %     if val2scan < 1
% %         idxs_lessNADHmito = [idxs_lessNADHmito; i];
% %     end
% % end
% % length(idxs_lessNADHmito)
% % 
% % %% Plot parameters, like scatter, and check which is uniquely down
% % % 
% % idxs_remain = [1:1000]';
% % % idxs_remain(idxs_remain == idxs_lessNADHmito) = [];
% % idxs_remain(ismember(idxs_remain,idxs_lessNADHmito)) = [];
% % % 
% % fig78_h = figure(78);
% % fig78_h.Position = [1921 257 1536 747];
% % for i = 1:50
% %     subplot(5,10,i)
% %     histogram(xRobustness(idxs_remain,i))
% %     hold on
% %     histogram(xRobustness(idxs_lessNADHmito,i))
% % end
% % legend('remain','lessNADHmito')
% % title('1 to 50')
% % % 
% % fig79_h = figure(79);
% % fig79_h.Position = [1921 257 1536 747];
% % for i = 51:100
% %     subplot(5,10,i-50)
% %     histogram(xRobustness(idxs_remain,i))
% %     hold on
% %     histogram(xRobustness(idxs_lessNADHmito,i))
% % end
% % legend('remain','lessNADHmito')
% % title('51 to 100')
% % % 
% % fig80_h = figure(80);
% % fig80_h.Position = [1921 257 1536 747];
% % for i = 101:150
% %     subplot(5,10,i-100)
% %     histogram(xRobustness(idxs_remain,i))
% %     hold on
% %     histogram(xRobustness(idxs_lessNADHmito,i))
% % end
% % legend('remain','lessNADHmito')
% % title('101 to 150')
% % %%
% % % 
% % xRobustness_remain = xRobustness(idxs_remain,:);
% % xRobustness_lessNADHmito = xRobustness(idxs_lessNADHmito,:);
% % % 
% % mean_remain = mean(xRobustness_remain);
% % mean_lessNADHmito = mean(xRobustness_lessNADHmito);
% % std_remain = std(xRobustness_remain);
% % std_lessNADHmito = std(xRobustness_lessNADHmito);
% % % 
% % fh41 = figure(41); plot(mean_remain - mean_lessNADHmito,'o-')
% % fh42 = figure(42); errorbar(mean_remain - mean_lessNADHmito, std_remain - std_lessNADHmito, 'o-')
% % % main parameters that clearly differ are x92 and x132
% % % mean_remain([92 132]) = [-0.6207    1.9980]
% % % mean_lessNADHmito([92 132]) = [-0.6615    1.9727] -> both are smaller in
% % % the cases of lessNADHmito
% % % p.mitoNADHVmax  = 1     * 10 .^ x(92);%100;
% % % p.mitoNADHKm    = 0.1   * 10 .^ x(93);
% % %     % temporary to account for the growth rate dependent mitoNADH
% % %     setup.ratioNADH = p.mitoNADHVmax;
% % %     setup.x132 = x(132);
% % % Basically, this means less capacity for NADH mitochondrial, regeneration,
% % % both regariding Vmax and the growth rate dependent icrease






%% safecopy
% % % % 
% % % % 
% % % % %% via PSA, set it higher in the reference model
% % % % % 
% % % % 
% % % % 
% % % % %% testing
% % % % % % first
% % % % % met_to_plot = 18; %17 %1
% % % % % met_to_plot = 8; %17 %1
% % % % % met_to_plot = 26;%4;%3;
% % % % % met_to_plot = 4;
% % % % % met_to_plot = 5;
% % % % % second
% % % % % met_to_plot = 26;
% % % % % met_to_plot = 20;
% % % % met_to_plot = 5;%5;
% % % % 
% % % % figure(1111)
% % % % for i = 1:1000
% % % % % for i = 1:500
% % % %     plot(simRobust{i}.gs.T, simRobust{i}.gs.Y(:,met_to_plot),...
% % % %         'linewidth',1,'color',[.75 .75 .75])
% % % %     if i == 1
% % % %         hold on
% % % %     end
% % % % end
% % % % plot(refSim{1}.gs.T(3002:end), refSim{1}.gs.Y(3002:end,met_to_plot),...
% % % %     'linewidth',1,'color','black')
% % % % hold off
% % % % %%
% % % % % met_to_plot = 6;
% % % % % met_to_plot = 17;
% % % % % met_to_plot = 14;
% % % % met_to_plot = 5;
% % % % % 
% % % % figure(1112)
% % % % for i = 1:1000
% % % % % for i = 1:500
% % % %     plot(canelas_SS.mtlD.D, simRobust{i}.ss.Y(:,met_to_plot),...
% % % %         'linewidth',1,'color',[.75 .75 .75])
% % % %     if i == 1
% % % %         hold on
% % % %     end
% % % % end
% % % % temp_plot = zeros(1,8);
% % % %     temp_plot(1) = refSim{1}.ss.Y.ss2can{1}(end,met_to_plot);
% % % %     temp_plot(2) = refSim{1}.ss.Y.ss2can{2}(end,met_to_plot);
% % % %     temp_plot(3) = refSim{1}.ss.Y.ss2can{3}(end,met_to_plot);
% % % %     temp_plot(4) = refSim{1}.ss.Y.ss2can{4}(end,met_to_plot);
% % % %     temp_plot(5) = refSim{1}.ss.Y.ss2can{5}(end,met_to_plot);
% % % %     temp_plot(6) = refSim{1}.ss.Y.ss2can{6}(end,met_to_plot);
% % % %     temp_plot(7) = refSim{1}.ss.Y.ss2can{7}(end,met_to_plot);
% % % %     temp_plot(8) = refSim{1}.ss.Y.ss2can{8}(end,met_to_plot);
% % % % plot(canelas_SS.mtlD.D, temp_plot,...
% % % %     'linewidth',1,'color','black')
% % % % hold off
% % % % 
% % % % 
% % % % %% gradually adding parameters in the change, to detect the issue
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% % % % N = 1000; %test multiple combinations
% % % % % 10 % noise
% % % % rng(1), rand_dist = rand(N,length(x_2nd)) * (log10(1.1)-log10(0.9)) + log10(0.9); % 10% noise
% % % % % 20 + (150-20) 
% % % % % figure, histogram(rand_dist)
% % % % 
% % % % % select only some paraemters to change
% % % % parsGLT = 35:36;
% % % % parsHXK = 28:34;
% % % % parsPFK = 43:56;
% % % % parsPGI = 57:60;
% % % % parsFBA = 11:15;
% % % % 
% % % % parsTPI = 79:82;
% % % % parsPGM = 67:70;
% % % % parsENO = 16:19;
% % % % parsPYK = 71:77;
% % % % parsTDH = 21:27;
% % % % 
% % % % parsPGK = 61:66;
% % % % parsGPD = 96:104;
% % % % parsPDC = 39:42;
% % % % parsADH = 1:10;
% % % % 
% % % % % additions rest
% % % % parsHOR2 = 105:107;
% % % % 
% % % % % addition tre cycle
% % % % parsPGM1 = 83:86;
% % % % parsTPS1 = 124:128;
% % % % parsTPS2 = 119:121;
% % % % parsNTH1 = 122:123;
% % % % parsTRE = [parsPGM1, parsTPS1, parsTPS2, parsNTH1];
% % % % 
% % % % % addition transporters
% % % % parsTransporters = [95, 108];
% % % % 
% % % % % additions cofactors
% % % % parsMitoNADH = 92:93;
% % % % parsMitoATP = 109:111;
% % % % parsATPase = 129;
% % % % parsADK1 = 130:131;
% % % % parsVac = 37;
% % % % parsCofactors = [parsMitoNADH, parsMitoATP, parsATPase, parsADK1, parsVac];
% % % % 
% % % % % addition inosine
% % % % parsInosine = 112:118;
% % % % 
% % % % % additions cofactors, growth rate-dependent
% % % % parsMitoNADH_mu = 132;
% % % % parsMitoATP_mu = 133;
% % % % parsATPase_mu = 134;
% % % % parsADK1_mu = 135:136;
% % % % parsCofactors_mu = [parsMitoNADH_mu, parsMitoATP_mu, parsATPase_mu, parsADK1_mu];
% % % % 
% % % % % creating selPars
% % % % selPars = [parsGLT, parsHXK, parsPFK, parsPGI, parsFBA,...
% % % %     parsTPI, parsPGM, parsENO, parsPYK, parsTDH,...
% % % %     parsPGK, parsGPD, parsPDC, parsADH];
% % % % selPars_HOR2 = [selPars, parsHOR2];
% % % % selPars_TRE = [selPars, parsTRE];
% % % % selPars_Transporters = [selPars, parsTransporters];
% % % % selPars_Cofactors = [selPars, parsCofactors];
% % % % selPars_Inosine = [selPars, parsInosine];
% % % % selPars_Cofactors_mu = [selPars, parsCofactors_mu];
% % % % 
% % % % % selPars
% % % % x_selPars = x_2nd; x_selPars(selPars) = x_selPars(selPars)+  rand_dist(1,selPars);
% % % % x_selPars_HOR2 = x_2nd; x_selPars_HOR2(selPars_HOR2) = x_selPars(selPars_HOR2)+  rand_dist(1,selPars_HOR2);
% % % % x_selPars_TRE = x_2nd; x_selPars_TRE(selPars_TRE) = x_selPars(selPars_TRE)+  rand_dist(1,selPars_TRE);
% % % % x_selPars_Transporters = x_2nd; x_selPars_Transporters(selPars_Transporters) = x_selPars(selPars_Transporters)+  rand_dist(1,selPars_Transporters);
% % % % x_selPars_Cofactors = x_2nd; x_selPars_Cofactors(selPars_Cofactors) = x_selPars(selPars_Cofactors)+  rand_dist(1,selPars_Cofactors);
% % % % x_selPars_Inosine = x_2nd; x_selPars_Inosine(selPars_Inosine) = x_selPars(selPars_Inosine)+  rand_dist(1,selPars_Inosine);
% % % % x_selPars_Cofactors_mu = x_2nd; x_selPars_Cofactors_mu(selPars_Cofactors_mu) = x_selPars(selPars_Cofactors_mu)+  rand_dist(1,selPars_Cofactors_mu);
% % % % 
% % % % % simulations
% % % % setup.plotResultsMode = 0;
% % % % [refSim] = simulateY3M1_separateGPSS(x_2nd, canelas_SS, data, dataset, setup);
% % % % [refSim_selPars] = simulateY3M1_separateGPSS(x_selPars, canelas_SS, data, dataset, setup);
% % % % [refSim_selPars_HOR2] = simulateY3M1_separateGPSS(x_selPars_HOR2, canelas_SS, data, dataset, setup);
% % % % [refSim_selPars_TRE] = simulateY3M1_separateGPSS(x_selPars_TRE, canelas_SS, data, dataset, setup);
% % % % [refSim_selPars_Transporters] = simulateY3M1_separateGPSS(x_selPars_Transporters, canelas_SS, data, dataset, setup);
% % % % [refSim_selPars_Cofactors] = simulateY3M1_separateGPSS(x_selPars_Cofactors, canelas_SS, data, dataset, setup);
% % % % [refSim_selPars_Inosine] = simulateY3M1_separateGPSS(x_selPars_Inosine, canelas_SS, data, dataset, setup);
% % % % [refSim_selPars_Cofactors_mu] = simulateY3M1_separateGPSS(x_selPars_Cofactors_mu, canelas_SS, data, dataset, setup);
% % % % setup.plotResultsMode = 0;
% % % % 
% % % % %% plotting
% % % % 
% % % % % %% selPars
% % % % setup.plotResultsMode = 30;
% % % % namesHits = cell(1,2); namesHits{1} = 'refSim'; namesHits{2} = 'refSim,selPars';
% % % % setup.namesHits = namesHits;
% % % % simRes = [refSim, refSim_selPars];
% % % % [~] = plotAll_Y3M1(simRes,legenda,canelas_SS,data,setup,[x_2nd;x_2nd]);
% % % % setup.plotResultsMode = 0;
% % % % 
% % % % % %% selPars_HOR2
% % % % setup.plotResultsMode = 30;
% % % % namesHits = cell(1,2); namesHits{1} = 'refSim'; namesHits{2} = 'refSim,selPars,HOR2';
% % % % setup.namesHits = namesHits;
% % % % simRes = [refSim, refSim_selPars_HOR2];
% % % % [~] = plotAll_Y3M1(simRes,legenda,canelas_SS,data,setup,[x_2nd;x_2nd]);
% % % % setup.plotResultsMode = 0;
% % % % 
% % % % % %% selPars_TRE
% % % % setup.plotResultsMode = 30;
% % % % namesHits = cell(1,2); namesHits{1} = 'refSim'; namesHits{2} = 'refSim,selPars,TRE';
% % % % setup.namesHits = namesHits;
% % % % simRes = [refSim, refSim_selPars_TRE];
% % % % [~] = plotAll_Y3M1(simRes,legenda,canelas_SS,data,setup,[x_2nd;x_2nd]);
% % % % setup.plotResultsMode = 0;
% % % % 
% % % % % %% selPars_Transporters
% % % % setup.plotResultsMode = 30;
% % % % namesHits = cell(1,2); namesHits{1} = 'refSim'; namesHits{2} = 'refSim,selPars,Transporters';
% % % % setup.namesHits = namesHits;
% % % % simRes = [refSim, refSim_selPars_Transporters];
% % % % [~] = plotAll_Y3M1(simRes,legenda,canelas_SS,data,setup,[x_2nd;x_2nd]);
% % % % setup.plotResultsMode = 0;
% % % % 
% % % % % %% selPars_Cofactors
% % % % setup.plotResultsMode = 30;
% % % % namesHits = cell(1,2); namesHits{1} = 'refSim'; namesHits{2} = 'refSim,selPars,Cofactors';
% % % % setup.namesHits = namesHits;
% % % % simRes = [refSim, refSim_selPars_Cofactors];
% % % % [~] = plotAll_Y3M1(simRes,legenda,canelas_SS,data,setup,[x_2nd;x_2nd]);
% % % % setup.plotResultsMode = 0;
% % % % 
% % % % % %% selPars_Inosine
% % % % setup.plotResultsMode = 30;
% % % % namesHits = cell(1,2); namesHits{1} = 'refSim'; namesHits{2} = 'refSim,selPars,Inosine';
% % % % setup.namesHits = namesHits;
% % % % simRes = [refSim, refSim_selPars_Inosine];
% % % % [~] = plotAll_Y3M1(simRes,legenda,canelas_SS,data,setup,[x_2nd;x_2nd]);
% % % % setup.plotResultsMode = 0;
% % % % 
% % % % % %% selPars_Cofactors_mu
% % % % setup.plotResultsMode = 30;
% % % % namesHits = cell(1,2); namesHits{1} = 'refSim'; namesHits{2} = 'refSim,selPars,Cofactors,mu';
% % % % setup.namesHits = namesHits;
% % % % simRes = [refSim, refSim_selPars_Cofactors_mu];
% % % % [~] = plotAll_Y3M1(simRes,legenda,canelas_SS,data,setup,[x_2nd;x_2nd]);
% % % % setup.plotResultsMode = 0;
% % % % 
% % % % 
% % % % %%





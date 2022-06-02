
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
x_2nd = x;

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


% % % % % % PEFULL_R90A_HARMONIZEAXP
% % % % % Get the right parameter set.
% % % % % Change starting AXP_GP to experimental data GP and not SS.
% % % % % Harmonize ADP and AMP profiles
% % % % clear, %close all
% % % % set_paths;
% % % % path(path, strcat(folder,'/tempResults'))
% % % % dbstop if error
% % % % colorpalette;
% % % % % get the starting parameters, datasets and legend
% % % % [legenda] = legendaFull;
% % % % loadData_Y3M1;
% % % % % simulation setup
% % % % setup.runSScanelas = 1; % SS sims cancelled for now. Only GP to be seen
% % % % setup.runGSvanHeerden = 1;
% % % % setup.plotResultsMode = 0; %30; %20; %10;
% % % % setup.conditionsSS = 0;
% % % % setup.conditionsGP = 0;
% % % % % model structure
% % % % setup.biomass = 1; % biomass
% % % % setup.clamp10.Pi = 0;%1; 
% % % % setup.clamp10.TRE = 0;%1; 
% % % % setup.clamp10.NADX = 0;%1; 
% % % % setup.clamp10.AXP = 0;%1; 
% % % % setup.clamp10.IXP = 0; %1; 
% % % % % setup updated to last topology
% % % % setup.latestTopology = 1;
% % % % setup.lastErrorsCorrect = 1;
% % % % % prior parameter adjustments
% % % % for prior = 1
% % % %     load('x_step3_3TREreg.mat','x_step3_3TREreg'); pset2 = x_step3_3TREreg;
% % % %     load('x_step3_1highMU_selected.mat'); pset2(137:146) = x_highMU;
% % % %     pars2zeros = [...%38 144,...
% % % %         28 29 30 31 139 140,...
% % % %         46 49 50 51 87,...
% % % %         12,14,...
% % % %         80,81,...
% % % %         71,76,...
% % % %         ];
% % % %     pset2(pars2zeros) = zeros;
% % % %     x = pset2;
% % % %         % 2021-02-22 changed: lowering the intensity of glucose entry via GLT
% % % %         x(35) = -0.6250;
% % % %         % 2021-02-25 changed: increased trehalose cycle uptake (from second test)
% % % %         x([83    86   124   125   126   127   128   119   120    33]) = [...
% % % %             1.7707   -0.8062   -0.3905   -0.7950    0.7822   -0.2262   -0.5143    0.1574   -0.3881   -1.5375,...
% % % %             ]; % option B
% % % %     %     x([83    86   124   125   126   127   128   119   120    33]) = [...
% % % %     %         1.8444    1.9588   -0.3952   -0.8152    0.8085   -0.3085   -0.5936    0.1164   -0.3871   -1.5606,...
% % % %     %         ]; % option C
% % % %         % 2021-02-25 changed: adjusted nth1 kcat
% % % %         x(123) = -0.9528 - 0.6;
% % % %         x1 = x;
% % % %         % 2021-03-04 changed: adjustment on pgi-pfk-fba-tpi
% % % %         idxs_change =   [11    13    15    43    44,...    
% % % %             45    47    48    52    53,...
% % % %             54    55    56    57    58,...    
% % % %             59    60    79    82];
% % % %         x(idxs_change) = [...
% % % %             -2.9923    0.4644   -0.0912   -1.3570   -0.0428,...
% % % %             0.9930   -0.0624   -1.3087    1.2157    0.9114,...
% % % %             0.0475   -0.4054    1.8658    0.4785    1.3032,...
% % % %             1.5172   -0.7006    0.0404   -0.1552];
% % % %         x2 = x;
% % % %         % 2021-03-09 changed: IXP influence lowered for correct ATP start +
% % % %         % ATPase
% % % %             % ATPase possibility for effect
% % % %             setup.Katpase_fermIncrease = 1;
% % % %             setup.Katpase_fermIncrease2 = 0;
% % % %             % ATPase parameter number 
% % % %             % -> x109
% % % %             % Start ATP
% % % %             setup.startATPgp = 1;
% % % %             % lower intake by IXP
% % % %             x(113) = .9;
% % % %         % 2021-03-11 changed: ATP profile fitted by pEst (78b, line720)
% % % %         x([43 44 47 48 49 53, 130 131, 113]) = [...
% % % %             -1.1530   -0.0000   -0.1086   -1.1062,...
% % % %             0.1840    0.6468   -0.3382    0.0761, 1.4793];
% % % %         % 2021-03-12 changed: start glk to fit the SS rates (pE79a, line820)
% % % %         x([11 12 13 14 15 79 80 81 82 140]) = [...
% % % %             -2.4033    0.2485    0.1829    0.2502    0.1458...
% % % %              1.3006    0.1197    0.2267   -0.6856    1.0000];
% % % %             % x140 could also be set to zero, but then the rates should be
% % % %             % better fit in another way.
% % % % end
% % % % load('pset_SS_GP_20210406.mat','x_SSsel','x_GPsel');
% % % % load('parameters_step1.mat');
% % % % 
% % % % % visualization color
% % % % % setup.plotBackColor = [.5 .5 1]; % blue
% % % % % setup.plotBackColor = [.5 1 .5]; % green
% % % % % setup.plotBackColor = [1 .5 .5]; % red
% % % % setup.plotBackColor = [.94 .94 .94]; % standard grey
% % % % % latest changes (2020 - 08 - 13)
% % % % setup.importPiON = 1;
% % % % setup.treCycleON = 1;
% % % % setup.startParamsOFF = 1;
% % % % setup.IXPdevelopmentOFF = 0; 
% % % % setup.adjust_mito_atpase = 1;
% % % % setup.Amd_IXP = 0;
% % % % setup.clamp_vALDgp = 0;
% % % % setup.NADHrecycle_syncPYR = 0;
% % % % setup.latest_IXP_change = 1;
% % % % setup.pHchangeON = 0;
% % % % % latest changes (2020 - 11 - 20)
% % % % setup.GPcorrectStart = 1;
% % % % % latest changes (2021 - 04 - 16)
% % % % recallWalther2010data;
% % % % setup.IXPstartConcs_w2010 = 1;
% % % % setup.ATPgpissue = 1;
% % % % % save('pset_SS_GP_20210406.mat','x_SSsel','x_GPsel');
% % % % load('pset_SS_GP_20210406.mat','x_SSsel','x_GPsel');
% % % % % saved variables: x_SSsel, x_GPsel
% % % % x = x_SSsel;
% % % % %     x([149 148 150 151]) = [1.0005    1.4666    1.5084    1.5043]; % 59    52    79    22
% % % % %     x(147) = x(11);
% % % %     x(147:151) = [-1.9033    1.4666    1.5084    1.5043    1.0005];
% % % %     % xAll1_pC28(4,[11 52 59 79 22])
% % % %     % xAll1_pC28(4,[147 148 149 150 151])
% % % %     % 
% % % %     % ans =
% % % %     % 
% % % %     %    -1.9033    1.4666    1.5084    1.5043    1.0005
% % % % 
% % % % %     x(147) = x(11);     % fba  (11:15)
% % % % %     x(148) = x(52);     % pfk  (43:56)
% % % % %     x(149) = x(59);     % pgi  (57:60)
% % % % %     x(150) = x(79);     % tpi  (79:82)
% % % % %     x(151) = x(22);     % gapdh(21:27)
% % % % %     x(147) = x_GPsel(11);     % fba  (11:15)
% % % % %     x(148) = x_GPsel(52);     % pfk  (43:56)
% % % % %     x(149) = x_GPsel(59);     % pgi  (57:60)
% % % % %     x(150) = x_GPsel(79);     % tpi  (79:82)
% % % % %     x(151) = x_GPsel(22);     % gapdh(21:27)
% % % % setup.separateParsGPSS = 1;
% % % % x(57:60) = [0.5680    1.3728    0.8920   -0.7635];
% % % % % (2021 05 11)
% % % % x(83) = 2.0549;
% % % % % (2021 05 18) latest change on glt to fit GP concentrations better
% % % % load('tempParameterSave_20210512.mat')
% % % % x = x97b_up;
% % % % % (2021 05 20) 
% % % % setup.final_Vmito_implementation = 1;
% % % % % (2021 05 24)
% % % % x(109:111) = [0.0042    0.0252   -0.0762];
% % % % % x(152) = 0; x2 = x; x2(152) = 100;% IXP is off before pulse
% % % % x(152) = 100; % IXP is off before pulse
% % % % 
% % % % % (2021 09? 27) finally looking good:
% % % % load('parSet_104a.mat','xFinal');
% % % % x = xFinal;
% % % % x(71:77) = zeros;
% % % % 
% % % % % (2021 10 25) regularization UG
% % % % load('parSet_105b.mat','x105b');
% % % % x([43:56, 148]) = x105b([43:56, 148]); % parsPFK
% % % % x(73) = x105b(73); % parsPYK (kmfbp)
% % % % % 
% % % % Y3M1manus_labels
% % % % 
% % % % %% Reference simulation
% % % % setup.GLT_km_same = 1;




% 
setup.plotResultsMode = 0;
[simRes_reference] = simulateY3M1_separateGPSS(x, canelas_SS, data, dataset, setup);
setup.simGPdata = simRes_reference{1}.gs;
setup.simSSdata = simRes_reference{1}.ss;
setup.litParams = zeros(size(x));
setup.plotResultsMode = 0;


% % % % %% (1) Recall x108a estimate + simple simulation
% % % % setup.sameGLTHXK4all = 1;
% % % % load('x108a.mat','x_1st','x_2nd')
% % % % 
% % % % % %% late simulation
% % % % setup.plotResultsMode = 0;
% % % % [refSim] = simulateY3M1_separateGPSS(x_2nd, canelas_SS, data, dataset, setup);
% % % % setup.plotResultsMode = 0;



% % %% plotting
% setup.plotResultsMode = 30;
% namesHits = cell(1,2); namesHits{1} = 'start'; namesHits{2} = 'x108res';
% setup.namesHits = namesHits;
% simRes = [simRes_reference, refSim];
% [~] = plotAll_Y3M1(simRes,legenda,canelas_SS,data,setup,[x; x_2nd]);
% setup.plotResultsMode = 0;


% %% early PSA to determine if focusing on FBP (PFK) or LG mets (PYK)
% % prepare parametrs
% x73_psa = [x_2nd;   x_2nd; x_2nd; x_2nd;   x_2nd; x_2nd; x_2nd];
% x73_psa(2:end,73) = [-2 -1 -.5 .5 1 2];
% 
% % simulation and simulation
% setup.plotResultsMode = 30;
% namesHits = cell(1,7); namesHits{1} = 'x108res';
%     namesHits{2} = '-2'; namesHits{3} = '-1'; namesHits{4} = '-0.5';
%     namesHits{5} = '0.5'; namesHits{6} = '1'; namesHits{7} = '2';
% [sim_x73_psa] = simulateY3M1_separateGPSS(x73_psa, canelas_SS, data, dataset, setup);
% setup.plotResultsMode = 0;


% %% check Pi producting in the inositol cycle
% setup.IXPcycle_correctPi_production = 0;
% % simulations
% setup.plotResultsMode = 0;
% [sim_IXP_corrPi_0] = simulateY3M1_separateGPSS(x_2nd, canelas_SS, data, dataset, setup);
% setup.IXPcycle_correctPi_production = 1;
% [sim_IXP_corrPi_1] = simulateY3M1_separateGPSS(x_2nd, canelas_SS, data, dataset, setup);
% % plotting
% setup.plotResultsMode = 30;
% namesHits = cell(1,2); namesHits{1} = 'ixp.corrPi.0'; namesHits{2} = 'ixp.corrPi.1';
% setup.namesHits = namesHits;
% simRes = [sim_IXP_corrPi_0, sim_IXP_corrPi_1];
% [~] = plotAll_Y3M1(simRes,legenda,canelas_SS,data,setup,[x; x_2nd]);
% setup.plotResultsMode = 0;
% % 
% setup.IXPcycle_correctPi_production = 0;


%% parameter estimation
% blank and constant setup
setup.plotResultsMode = 0;
% blankWeight = zeros(1,63);
blankWeight = zeros(1,140);
lambdalist = 0;
setup.parEst.lambda = lambdalist(1); lam = setup.parEst.lambda;
% % % options = optimoptions('lsqnonlin','Display','iter','OutputFcn',{@saveIterationsMain},'PlotFcn',{@optimplotx});
% % options = optimoptions('lsqnonlin','Display','iter','OutputFcn',{@saveIterationsMain},'PlotFcn',{@optimplotx},...
% %     'FunctionTolerance', 1e-4, 'OptimalityTolerance', 1e-4);
% options = optimoptions('lsqnonlin','Display','iter',...
%     'OutputFcn',{@saveIterationsMain},'PlotFcn',{@optimplotx},...
%     'FiniteDifferenceStepSize',0.2);
options = optimoptions('lsqnonlin','Display','iter',...
    'FunctionTolerance', 1e-4, 'OptimalityTolerance', 1e-4);
%     'OutputFcn',{@saveIterationsMain},...
setup.parEst.costfun = 10015;

            
%% (5) Saving names.
ntests = 13;
names_cell = cell(1,ntests);
for i = 1:ntests
    names_cell{i} = sprintf('mF3_w%d.mat',i);
end
setup.litParams = x_2nd; 


%% ESTIMATION 
% 
selPars = 73;
% end
plength = length(selPars);
lb = -3*ones(1,plength);
ub = 3*ones(1,plength);
x_temp = x_2nd(selPars);

% select warray (first fitting PEP in SS and gradually GP)
w_PEP_SS = 9;
w_PEP_GP = 39;
w_arr = [1E-2 2E-2 5E-2 1E-1 2E-1 5E-1 1E0 2E0 5E0 1E1 2E1 5E1 1E2];
w_arr_flip = flip(w_arr); 
w_mat = zeros(ntests,length(blankWeight));
for i = 1:ntests
    w_mat(i,:) = blankWeight;
    w_mat(i,w_PEP_SS) = w_arr_flip(i);
    w_mat(i,w_PEP_GP) = w_arr(i);
end


%% loop lambda
% cluster = parcluster('local');
% pool = parpool(cluster,13);
% parfor j = 1:ntests
% 
%     % select weights
%     lam = 0;
%     w = w_mat(j,:);
% 
%     % select name
%     saveName = names_cell{j};
% 
%     % estimation
%     sprintf('enzyme = %d, lambda = %d', 1, j)
%     tic
%     [xres,resnorm,residual,exitflag] = lsqnonlin(@costfunSystemY3M1_temp_separateParsGPSS3,x_temp,lb,ub,options,canelas_SS,setup,x_2nd,data,dataset,w,lam,selPars);
%     t = toc; 
%     disp(xres);
%     
%     % simulation
%     xres_full = x_2nd;
%     xres_full(selPars) = xres;
%     [simRes_temp] = simulateY3M1_separateGPSS(xres_full, canelas_SS, data, dataset, setup);
% 
%     % error simulation
%     x2 = x_2nd; x2(selPars) = xres;
%     lam2 = 0;
%     [error]=costfunSystemY3M1_temp_separateParsGPSS3(xres,canelas_SS,setup,x2,data,dataset,w,lam2,selPars);
%     errorSim = sum(abs(error));
%     % error parameters
%     errorLam = sum(abs(xres));
% 
%     % saving result
% %         parsave_Y3M2_cluster_3(saveName,xres,resnorm,residual,exitflag,t,w,lam,simRes_temp)
%     parsave_Y3M2_cluster_4(saveName,xres,resnorm,residual,exitflag,t,w,lam,simRes_temp, error, errorSim, errorLam)
% 
% end
% delete(pool)


%% DATA RECALL (not the multistart cases yet)

% recalling data
nEnzs = 1;
xAll = cell(nEnzs,1+ntests);
namesHits = cell(nEnzs,1+ntests);
simRes = cell(nEnzs,1+ntests);
% recalling errors
errorAll = cell(nEnzs,1+ntests);
errorSim_mat = zeros(nEnzs,1+ntests);
errorPars_mat = zeros(nEnzs,1+ntests);
% 
for i = 1:nEnzs
%     tempLoadName = names_cell{i}; % tempLoadName -> enzyme
    for j = 1:(1+ntests)
        if j == 1 % if first case
            % 
            xAll{i,j} = x_2nd; % assign x
            namesHits{i,j} = 'initial'; % assign name
            simRes{i,j} = simRes_reference; % assign simRes
            % 
            errorAll{i,j} = 0;
            errorSim(i,j) = 0;
            errorPars(i,j) = 0;
        else
%             loadName = sprintf(tempLoadName,j-1); % LoadName -> lambdalist (j-1)
            loadName = names_cell{j-1};
            if exist(loadName,'file') == 2 % if it exists
                load(loadName);
                % 
%                 xtemp = x; xtemp(selParsStruct{i}) = xres; 
                xtemp = x_2nd; xtemp(selPars) = xres; 
                xAll{i,j} = xtemp; % assign x
                namesHits{i,j} = loadName; % assign name
                simRes{i,j} = simRes_temp{1}; % assign simRes
                % 
                errorAll{i,j} = error;
                errorSim_mat(i,j) = errorSim;
                errorPars_mat(i,j) = errorLam;
            end
        end
    end
end

%% visualization 3: actualy simulations
xDummy = [x; x; x; x; x; ...
            x; x; x; x; x; ...
            x; x; x; x];
for i = 1:nEnzs
    simRes{i,1} = simRes{i,1}{1};
end
% % % % From here on it will probably change


%% 108A
% idxsSel = 1:23; % all of them. None too bacd anyways
% idxsSel = 1:10; % until lambdalist(9), which seem to be the best case
% idxsSel = 1:12; % 13 and 14 already weird values 
% idxsSel = 1:12; % 13 and 14 already weird values 
setup.plotResultsMode = 30;
setup.namesHits = namesHits(1,:);
[~] = plotAll_Y3M1(simRes(1,:),legenda,canelas_SS,data,setup,xDummy);


%% analysze error via errorAll
errorAll_focus = errorAll(2:14);
len = length(errorAll_focus);
% find([1 0 1 1])
for i = 1:len
    find(errorAll_focus{i})
end
%% 
idxs_PEPss = [65    66    67    68    69    70    71    72];
idxs_PEPgp = [325   326   327   328   329   330   331   332   333   334   335];
%%
% 
error_PEPss = zeros(1,len);
error_PEPgp = zeros(1,len);
% 
for i = 1:len
    error_PEPss(i) = sum(abs(errorAll_focus{i}(idxs_PEPss)));
    error_PEPgp(i) = sum(abs(errorAll_focus{i}(idxs_PEPgp)));
end







% % % % %% 
% % % % simRes_mF3 = simRes;
% % % % save('mF3sims.mat','error_PEPss','error_PEPgp','simRes_mF3','legenda','canelas_SS','data','setup','xDummy','namesHits')


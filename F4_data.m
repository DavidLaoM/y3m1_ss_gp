% % F4_DATA.m
% This code reproduces the plots that highlight the relevance of including
% different types of datasets for more accurate parameter estimates and
% simulations.
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
% [refSim] = simulateY3M1_separateGPSS(x_2nd, canelas_SS, data, dataset, setup);
[refSim] = simulateY3M1_separateGPSS(x_3rd, canelas_SS, data, dataset, setup);
setup.plotResultsMode = 0;

% % %% plotting
% setup.plotResultsMode = 30;
% % namesHits = cell(1,2); namesHits{1} = 'start'; namesHits{2} = 'x108res';
% namesHits = cell(1,1); namesHits{1} = 'x108res'; %namesHits{2} = 'x108res';
% setup.namesHits = namesHits;
% % simRes = [simRes_reference, refSim];
% simRes = refSim;
% % %%
% [~] = plotAll_Y3M1(simRes,legenda,canelas_SS,data,setup,x_2nd);
% setup.plotResultsMode = 0;


% Intention on take-aways from this section:
% (1) Missing information if only one data type is missing is negative.
% (2) Better parameter solution space is achieved.
% For which
% (1) Compare different error in simulations and in bar plots
% (2) Compare error and confidence interval change
% In this file
% (1) Contribution proteomics: turn it on and off
% (2) Comparison parameter fits SS AND/OR GP: fit parameters + Error + CI.
% (3) Comparison parameter fits with AND/OR without TRE: fit parameters + Error + CI.
% (4) Comparison parameter fits metabolomics AND/OR fluxomics: fit parameters + Error + CI
% (5) Develop diagram
    % For now, HXK as a case study, since high influence under all the
    % cases described
% Parameters estimated in the following way:
    % (option) - Start adding sme noise to estimate. More as global sampling.
    % - Estimate parameters, using changes in cost functions.
    % - Covariance method to find confidence intervals.


%% (1) Contribution proteomics: turn it on and off
% simulate
setup.vanHoek_off = 1;
setup.plotResultsMode = 0;
% [PvHoek_off_Sim] = simulateY3M1_separateGPSS(x_2nd, canelas_SS, data, dataset, setup);
[PvHoek_off_Sim] = simulateY3M1_separateGPSS(x_3rd, canelas_SS, data, dataset, setup);
setup.plotResultsMode = 0;
setup.vanHoek_off = 0;

%% plot simulations
setup.plotResultsMode = 30;
namesHits = cell(1,2); namesHits{1} = 'reference'; namesHits{2} = 'PvHoek off';
simRes = [refSim, PvHoek_off_Sim];
%
% [~] = plotAll_Y3M1(simRes,legenda,canelas_SS,data,setup,[x_2nd; x_2nd]);
[~] = plotAll_Y3M1(simRes,legenda,canelas_SS,data,setup,[x_3rd; x_3rd]);
setup.plotResultsMode = 0;


%%
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% % % % % % % % PLOTS PROTEOMICS  % % % % % % % % % % % % % % % % % % % % % 
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
%% selecting relevant data
% checkup
if exist('f1001_h','var')
    clf(f1001_h)
end
% recall info
f2_h = figure(2);
f2_r_h = get(f2_h,'children');
mu = f2_r_h(42).Children(1).XData;
HXK_exp = f2_r_h(42).Children(1).YData;
HXK_sim = f2_r_h(42).Children(3).YData;
HXK_sim_protOFF = f2_r_h(42).Children(2).YData;
% calculate error (normalized)
HXK_sim_error_full = (HXK_sim - HXK_exp);
HXK_sim_protOFF_error_full = (HXK_sim_protOFF - HXK_exp);
HXK_sim_error = sum(abs(HXK_sim_error_full ./ HXK_exp));
HXK_sim_protOFF_error = sum(abs(HXK_sim_protOFF_error_full ./ HXK_exp));

%% making the plot (subplots just needed for the final one)

% figure
f1001_h = figure(1001);
% f1001_h.Position = [1978 275 1227 710];
f1001_h.Position = [2425 275 780 710];
clf(1001)

% simulated vs experimental
sp2 = subplot(2,2,1);
% 
line([0 2.25],[0 2.25],'color','k','LineStyle','--')
hold on
% 
sc1 = scatter(HXK_exp,HXK_sim_protOFF);
sc1.LineWidth = 0.6;
sc1.MarkerEdgeColor = 'k';
sc1.MarkerFaceColor = [.7 .7 .7];%'r';
hold on
% 
sc2 = scatter(HXK_exp,HXK_sim);
sc2.LineWidth = 0.6;
sc2.MarkerEdgeColor = 'k';
sc2.MarkerFaceColor = 'k';%'b';
% 
sp2.XLim = [0 2.25];
sp2.YLim = [0 2.25];
% %     % 
% %     line(sp2.XLim,sp2.YLim,'color','k','LineStyle','--')
% 
xlabel('Experimental Reaction Rate (mmol L-1 s-1)')
ylabel('Simulated Reaction Rate (mmol L-1 s-1)')
title('Experimental vs Simulation plot')
% %     legend('simulation, no proteomics','simulation proteomics',...
% %         'location','southeast')
box on

% error plots 2
% c2 = categorical({'error sim NO prot','error sim'});
c2 = categorical({'sim + Prot','sim - prot'});
sp4 = subplot(2,2,3);
bp2 = bar(c2,[sum(abs(HXK_sim_error_full)), sum(abs(HXK_sim_protOFF_error_full))],'FaceColor','flat');
sp4.YLim = [0 3.5];
bp2.CData(1,:) = [0 0 0];%[0 0 .7];
bp2.CData(2,:) = [.7 .7 .7];%[1 0 0];
% legend('error sim NO prot','error sim','location','northwest')
ylabel('error')
title('error plots, only part')

% simulated vs experimental + error plots
a1 = subplot(2,2,2);
copyobj(get(sp2,'children'),a1);
% 
xlabel('Experimental Reaction Rate (mmol L-1 s-1)')
ylabel('Simulated Reaction Rate (mmol L-1 s-1)')
%     title('proposed final plot','color','red')
%     legend('simulation, no proteomics','simulation proteomics',...
%         'location','southeast')
% 
a1.XLim = sp2.XLim;
a1.YLim = sp2.YLim;
box on
% addition
%     handaxes1 = axes('Position', [0.695 0.8125 a1.Position(3)/3.2 a1.Position(4)/3.2]);
handaxes1 = axes('Position', [0.81 0.6 a1.Position(3)/2.5 a1.Position(4)/2.5]);
copyobj(get(sp4,'children'),handaxes1);
handaxes1.XLim = [0 3];
handaxes1.YLim = sp4.YLim; %[0 4];
set(gca,'XTickLabel',[],'YTickLabel',[]);
box on
txt1 = text(0.2,2.75,'error','FontSize',12);
% txt2 = text(-3.5,7,'B','FontSize',20);

    % 
% % % % print(gcf,'-dpdf','E:\OneDrive - TU Eindhoven\Documents\ch4_Y3M1\results\figures\manuscript\F4B');


%%
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% % % % % % % % PLOTS CONFIDENCE INTERVALS  % % % % % % % % % % % % % % % % 
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% selecting relevant data
for intermediate_stuff = 1
% blank and constant setup
setup.plotResultsMode = 0;
blankWeight = zeros(1,140);
lambdalist = 0;
setup.parEst.lambda = lambdalist(1); lam = setup.parEst.lambda;
options = optimoptions('lsqnonlin','Display','iter',...
    'FunctionTolerance', 1e-4, 'OptimalityTolerance', 1e-4);
setup.parEst.costfun = 10015;
% 
ntests = 12;

% parameters 
parsHXK = 28:34;
% 
parComb_cell = cell(1,ntests); 
for i = 1:ntests
    parComb_cell{i} = parsHXK;
end

% saving names
names_cell = cell(1,ntests);
% 
names_cell{1} = 'mF6_2ssgp_1control.mat';
names_cell{2} = 'mF6_2ssgp_2ss.mat';
names_cell{3} = 'mF6_2ssgp_3gp.mat';
names_cell{4} = 'mF6_2ssgp_4ssgp.mat';
% 
names_cell{5} = 'mF6_3tre_1control.mat';
names_cell{6} = 'mF6_3tre_2core.mat';
names_cell{7} = 'mF6_3tre_3tre.mat';
names_cell{8} = 'mF6_3tre_4coretre.mat';
% 
names_cell{9} =  'mF6_4metrate_1control.mat';
names_cell{10} = 'mF6_4metrate_2met.mat';
names_cell{11} = 'mF6_4metrate_3rate.mat';
names_cell{12} = 'mF6_4metrate_4metrate.mat';
% 
setup.litParams = zeros(size(x)); 

% weight array
w2ss = 15; % ss HXK
w2gp = 48; % gp HXK
w3core = 37; % gp G6P
w3tre = 45; % gp T6P
w4met = 37; % gp G6P
w4rate = 48; % gp HXK
% 
warray_cell = cell(1,ntests);
for i = 1:ntests
    warray_cell{i} = blankWeight;
end
% 
warray_cell{2}(w2ss) = ones;
warray_cell{3}(w2gp) = ones;
warray_cell{4}([w2ss, w2gp]) = ones;
% 
warray_cell{6}(w3core) = ones;
warray_cell{7}(w3tre) = ones;
warray_cell{8}([w3core, w3tre]) = ones;
% 
warray_cell{10}(w4met) = ones;
warray_cell{11}(w4rate) = ones;
warray_cell{12}([w4met, w4rate]) = ones;

%% DATA RECALL (not the multistart cases yet)
ntests = 12;
nEnzs = 1;

% recalling data
xAll = cell(nEnzs,1+ntests);
namesHits = cell(nEnzs,1+ntests);
simRes = cell(nEnzs,1+ntests);
% recalling errors
errorAll = cell(nEnzs,1+ntests);
errorSim_mat = zeros(nEnzs,1+ntests);
errorPars_mat = zeros(nEnzs,1+ntests);
% 
residual_cell = cell(nEnzs,1+ntests);
output_cell = cell(nEnzs,1+ntests);
tempLam_cell = cell(nEnzs,1+ntests);
jacobian_cell = cell(nEnzs,1+ntests);

% 
for i = 1:nEnzs
%     tempLoadName = names_cell{i}; % tempLoadName -> enzyme
    for j = 1:(1+ntests)
        if j == 1 % if first case
            % 
            xAll{i,j} = x; % assign x
            namesHits{i,j} = 'initial'; % assign name
            simRes{i,j} = refSim{1}; % assign simRes
            % 
            errorAll{i,j} = 0;
            errorSim(i,j) = 0;
            errorPars(i,j) = 0;
        else
%             loadName = sprintf(tempLoadName,j-1); % LoadName -> lambdalist (j-1)
            loadName = names_cell{j-1}; % LoadName -> lambdalist (j-1)
            if exist(loadName,'file') == 2 % if it exists
                load(loadName);
                % 
                xtemp = x; xtemp(parComb_cell{j-1}) = xres; 
                % xtemp = x; xtemp(pars_all) = xres; 
                xAll{i,j} = xtemp; % assign x
                namesHits{i,j} = loadName; % assign name
                simRes{i,j} = simRes_temp{1}; % assign simRes
                % 
                errorAll{i,j} = error;
                errorSim_mat(i,j) = errorSim;
                errorPars_mat(i,j) = errorLam;
                %
                output_cell{i,j} = output;
                residual_cell{i,j} = residual;
                tempLam_cell{i,j} = tempLam;
                jacobian_cell{i,j} = jacobian;
            end
        end
    end
end
case1_xAll = xAll;

%% visualization 3: actualy simulations
xDummy = [x; x; x; x; x; ...
            x; x; x; x; ...
            x; x; x; x; ...
            ];
% simRes{1} = simRes{1}{1};
% % % % for i = 1:nEnzs
% % % %     simRes{i,1}{1} = simRes{i,1};
% % % % end
% % % % From here on it will probably change


%% confidence intervals
% (1) 
% option to use nlparci (https://nl.mathworks.com/help/stats/nlparci.html),
% but this needs more info (also: https://nl.mathworks.com/matlabcentral/answers/171048-how-can-i-get-the-cv-and-confidence-interval-for-my-parameter-estimates-from-output-of-lsqnonlin-fu).
% (2) Also the second method of direct calculation, to check.
%  

% %  for now i need to re run estimations and request more outputs.
% 2:5
temp_idx = 5;
parameters = xAll{temp_idx}(parsHXK);
residual = residual_cell{temp_idx};
Jacobian = jacobian_cell{temp_idx};
% 
ci = nlparci(parameters,residual,'Jacobian',Jacobian);

% 
nS = 12;
ci_cell = cell(1,nS);
for i = 1:nS
    % 
    temp_idx = 1 + i;
    parameters = xAll{temp_idx}(parsHXK);
    residual = residual_cell{temp_idx};
    Jacobian = jacobian_cell{temp_idx};
    % 
    ci_cell{i} = nlparci(parameters,residual,'Jacobian',Jacobian);
end

% % Narrative: can use the first story with HXK pretty well. Now for the
% others it does not look too good. Test in a second run with TPS1 the
% second case and GLT the third.

% % from
% warray_cell{2}(w2ss) = ones;
% % to
% warray_cell{4}([w2ss, w2gp]) = ones;
abs([ci_cell{1}(:,1) - ci_cell{1}(:,2),...
    ci_cell{2}(:,1) - ci_cell{2}(:,2),...
    ci_cell{3}(:,1) - ci_cell{3}(:,2),...
    ci_cell{4}(:,1) - ci_cell{4}(:,2)])
% 
% % ideally from w3 core to w3tre
% warray_cell{6}(w3core) = ones;
% warray_cell{7}(w3tre) = ones;
% warray_cell{8}([w3core, w3tre]) = ones;
% % tre to core+tre, but of course explosive
% % could try additionally TPS1
% 
% % 
% warray_cell{10}(w4met) = ones;
% warray_cell{11}(w4rate) = ones;
% warray_cell{12}([w4met, w4rate]) = ones;
% %  limited to hexokinase for simplicity of story
% % Could try additionally GLT


%% recall daata from v2
for recallData_2 = 1
    % parameters 
    parsHXK = 28:34;
    parsTPS1 = 124:128;
    parsGLT = [36,38];
    % 
    parComb_cell = cell(1,ntests); 
    % for i = 1:4
    %     parComb_cell{i} = parsHXK;
    % end
    % for i = 5:8
    %     parComb_cell{i} = parsTPS1;
    % end
    % for i = 9:12
    %     parComb_cell{i} = parsGLT;
    % end
    for i = 1:3
        parComb_cell{i} = parsTPS1;
    end
    for i = 4:6
        parComb_cell{i} = parsGLT;
    end

    % saving names
    names_cell = cell(1,ntests);
    % 
    % names_cell{1} = 'mF6_v2_2ssgp_1control.mat';
    % names_cell{2} = 'mF6_v2_2ssgp_2ss.mat';
    % names_cell{3} = 'mF6_v2_2ssgp_3gp.mat';
    % names_cell{4} = 'mF6_v2_2ssgp_4ssgp.mat';
    % 
    % names_cell{5} = 'mF6_v2_3tre_1control.mat';
    % names_cell{6} = 'mF6_v2_3tre_2core.mat';
    % names_cell{7} = 'mF6_v2_3tre_3tre.mat';
    % names_cell{8} = 'mF6_v2_3tre_4coretre.mat';
    names_cell{1} = 'mF6_v2_3tre_2core.mat';
    names_cell{2} = 'mF6_v2_3tre_3tre.mat';
    names_cell{3} = 'mF6_v2_3tre_4coretre.mat';
    % 
    % names_cell{9} =  'mF6_v2_4metrate_1control.mat';
    names_cell{4} = 'mF6_v2_4metrate_2met.mat';
    names_cell{5} = 'mF6_v2_4metrate_3rate.mat';
    names_cell{6} = 'mF6_v2_4metrate_4metrate.mat';
    % 
    setup.litParams = zeros(size(x)); 

    % weight array
    w2ss = 15; % ss HXK
    w2gp = 48; % gp HXK
    w3core = 37; % gp G6P
    w3tre = 45; % gp T6P
    w4met = 37; % gp G6P
    w4rate = 48; % gp HXK
    % 
    warray_cell = cell(1,ntests);
    for i = 1:ntests
        warray_cell{i} = blankWeight;
    end
    % 
    % warray_cell{2}(w2ss) = ones;
    % warray_cell{3}(w2gp) = ones;
    % warray_cell{4}([w2ss, w2gp]) = ones;
    % 
    warray_cell{1}(w3core) = ones;
    warray_cell{2}(w3tre) = ones;
    warray_cell{3}([w3core, w3tre]) = ones;
    % 
    warray_cell{4}(w4met) = ones;
    warray_cell{5}(w4rate) = ones;
    warray_cell{6}([w4met, w4rate]) = ones;



    % %%
    ntests = 6;
    nEnzs = 1;

    % recalling data
    xAll = cell(nEnzs,1+ntests);
    namesHits = cell(nEnzs,1+ntests);
    simRes = cell(nEnzs,1+ntests);
    % recalling errors
    errorAll = cell(nEnzs,1+ntests);
    errorSim_mat = zeros(nEnzs,1+ntests);
    errorPars_mat = zeros(nEnzs,1+ntests);
    % 
    residual_cell = cell(nEnzs,1+ntests);
    output_cell = cell(nEnzs,1+ntests);
    tempLam_cell = cell(nEnzs,1+ntests);
    jacobian_cell = cell(nEnzs,1+ntests);

    % 
    for i = 1:nEnzs
    %     tempLoadName = names_cell{i}; % tempLoadName -> enzyme
        for j = 1:(1+ntests)
            if j == 1 % if first case
                % 
                xAll{i,j} = x; % assign x
                namesHits{i,j} = 'initial'; % assign name
                simRes{i,j} = refSim{1}; % assign simRes
                % 
                errorAll{i,j} = 0;
                errorSim(i,j) = 0;
                errorPars(i,j) = 0;
            else
    %             loadName = sprintf(tempLoadName,j-1); % LoadName -> lambdalist (j-1)
                loadName = names_cell{j-1}; % LoadName -> lambdalist (j-1)
                if exist(loadName,'file') == 2 % if it exists
                    load(loadName);
                    % 
                    xtemp = x; xtemp(parComb_cell{j-1}) = xres; 
                    % xtemp = x; xtemp(pars_all) = xres; 
                    xAll{i,j} = xtemp; % assign x
                    namesHits{i,j} = loadName; % assign name
                    simRes{i,j} = simRes_temp{1}; % assign simRes
                    % 
                    errorAll{i,j} = error;
                    errorSim_mat(i,j) = errorSim;
                    errorPars_mat(i,j) = errorLam;
                    %
                    output_cell{i,j} = output;
                    residual_cell{i,j} = residual;
                    tempLam_cell{i,j} = tempLam;
                    jacobian_cell{i,j} = jacobian;
                end
            end
        end
    end


    % %% visualization 3: actualy simulations
    xDummy = [x; ...
                x; x; x; ...
                x; x; x; ...
                ];
    % simRes{1} = simRes{1}{1};
    % % % % for i = 1:nEnzs
    % % % %     simRes{i,1}{1} = simRes{i,1};
    % % % % end
    % % % % From here on it will probably change


    % %% 110A
    % % idxsSel = 1; % just reference simulation
    % % idxsSel = [1:17]; %
    % idxsSel_test1 = [1,2:4];
    % idxsSel_test2 = [1,5:7];
    % 
    % %% TPS1
    % close all
    % idxsSel_temp = idxsSel_test1;
    % setup.plotResultsMode = 30;
    % setup.namesHits = namesHits(idxsSel_temp);
    % [~] = plotAll_Y3M1(simRes(idxsSel_temp),legenda,canelas_SS,data,setup,xDummy(idxsSel_temp,:));
    % % looks as expected
    % 
    % %% GLT
    % close all
    % idxsSel_temp = idxsSel_test2;
    % setup.plotResultsMode = 30;
    % setup.namesHits = namesHits(idxsSel_temp);
    % [~] = plotAll_Y3M1(simRes(idxsSel_temp),legenda,canelas_SS,data,setup,xDummy(idxsSel_temp,:));
    % % looks as expected

    % %% confidence intervals
    % (1) 
    % option to use nlparci (https://nl.mathworks.com/help/stats/nlparci.html),
    % but this needs more info (also: https://nl.mathworks.com/matlabcentral/answers/171048-how-can-i-get-the-cv-and-confidence-interval-for-my-parameter-estimates-from-output-of-lsqnonlin-fu).
    % (2) Also the second method of direct calculation, to check.

    % %  draft
    % % 2:5
    % temp_idx = 5;
    % parameters = xAll{temp_idx}(parsHXK);
    % residual = residual_cell{temp_idx};
    % Jacobian = jacobian_cell{temp_idx};
    % % 
    % ci = nlparci(parameters,residual,'Jacobian',Jacobian)

    % 
    nS = 6;
    ci_cell2 = cell(1,nS);
    for i = 1:nS
        % 
        temp_idx = 1 + i;
        if i < 4
            parameters = xAll{temp_idx}(parsTPS1);
        else
            parameters = xAll{temp_idx}(parsGLT);
        end
        residual = residual_cell{temp_idx};
        Jacobian = jacobian_cell{temp_idx};
        % 
        ci_cell2{i} = nlparci(parameters,residual,'Jacobian',Jacobian);
    end
end

% %% nice trend for TPS1
abs([ci_cell2{1}(:,1) - ci_cell2{1}(:,2),...
    ci_cell2{2}(:,1) - ci_cell2{2}(:,2),...
    ci_cell2{3}(:,1) - ci_cell2{3}(:,2)])

% %% can be only narrated form rate to rate+met
abs([ci_cell2{4}(:,1) - ci_cell2{4}(:,2),...
    ci_cell2{5}(:,1) - ci_cell2{5}(:,2),...
    ci_cell2{6}(:,1) - ci_cell2{6}(:,2)])


% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
%% making the final plot

% 
name_cell_selected = cell(1,7);
    % 
    name_cell_selected{1} = 'TPS1. wG6P'; % classic data
    name_cell_selected{2} = 'TPS1. wTRE'; % new trehalose cycle data
    name_cell_selected{3} = 'TPS1. wG6P+wTRE';
    %
    name_cell_selected{4} = 'GLT. wHXK'; % rate data
    name_cell_selected{5} = 'GLT. wHXK+wG6P'; % combined rate and concentration data
    %
    name_cell_selected{6} = 'HXK. wHXK,SS'; % SS data
    name_cell_selected{7} = 'HXK. wHXK,SS+GP'; % SS+GP data

% 
ci_cell_selected = cell(1,7);
    % 
    ci_cell_selected{1} = ci_cell2{1}; % classic data
    ci_cell_selected{2} = ci_cell2{2}; % new trehalose cycle data
    ci_cell_selected{3} = ci_cell2{3};
    %
    ci_cell_selected{4} = ci_cell2{5}; % rate data
    ci_cell_selected{5} = ci_cell2{6}; % combined rate and concentration data
    %
    ci_cell_selected{6} = ci_cell{2}; % SS data
    ci_cell_selected{7} = ci_cell{4}; % SS+GP data
% 
errbar_cell_selected = cell(1,7);
for i = 1:7
    errbar_cell_selected{i} =abs(ci_cell_selected{i}(:,1) - ci_cell_selected{i}(:,2)); % classic data
end

%
pars_cell_selected = cell(1,7);
    % 
    pars_cell_selected{1} = xAll{2}(parsTPS1);
    pars_cell_selected{2} = xAll{3}(parsTPS1);
    pars_cell_selected{3} = xAll{4}(parsTPS1);
    % 
    pars_cell_selected{4} = xAll{6}(parsGLT);
    pars_cell_selected{5} = xAll{7}(parsGLT);
    % 
    pars_cell_selected{6} = case1_xAll{3}(parsHXK);
    pars_cell_selected{7} = case1_xAll{5}(parsHXK);
end


%% Making the plots
% checkup
if exist('fh_6002','var')
    clf(fh_6002)
end
% clf(6002)
fh_6002 = figure(6002);
% fh_6002.Position = [1986 593 1375 326];
fh_6002.Position = [1986 307 657 612];
% %% 
dp1 = subplot(2,2,1);
% % errorbar(pars_cell_selected{1}(1:4),[1:4],errbar_cell_selected{1}(1:4),...
% %     'horizontal','o','MarkerFaceColor',[0 0.4470 0.7410],'color',[0 0.4470 0.7410],...
% %     'MarkerEdgeColor',[0 0.4470 0.7410],'CapSize',0,'LineWidth',1)
c_grayish = [.6 .6 .6];
errorbar(pars_cell_selected{1}(1:4),[1:4],errbar_cell_selected{1}(1:4),...
    'horizontal','o','MarkerFaceColor',c_grayish,'color',c_grayish,...
    'MarkerEdgeColor',c_grayish,'CapSize',0,'LineWidth',1)
hold on
% % errorbar(pars_cell_selected{2}(1:4),[1:4]-1/3,errbar_cell_selected{2}(1:4),'horizontal','o','MarkerFaceColor',[0.8500    0.3250    0.0980],'color',[0.8500    0.3250    0.0980],'MarkerEdgeColor',[0.8500    0.3250    0.0980],'CapSize',0,'LineWidth',1)
% errorbar(pars_cell_selected{3}(1:4),[1:4]-2/3,errbar_cell_selected{3}(1:4),...
%     'horizontal','o','MarkerFaceColor',[0.9290    0.6940    0.1250],'color',[0.9290    0.6940    0.1250],...
%     'MarkerEdgeColor',[0.9290    0.6940    0.1250],'CapSize',0,'LineWidth',1)
errorbar(pars_cell_selected{3}(1:4),[1:4]-1/2,errbar_cell_selected{3}(1:4),...
    'horizontal','o','MarkerFaceColor','k','color','k',...
    'MarkerEdgeColor','k','CapSize',0,'LineWidth',1)
% title('TPS1. Using new measured metabolites (trehalose)')
tempXlim = [-8 8];
xlim(tempXlim)
% ylim([-0.5 5.5])
% ylim([-1 5])
ylim([-1 6])

% yticks([1 2 3 4 5])
yticks([0+1/3/2+0.5 1+1/3/2+0.5 2+1/3/2+0.5 3+1/3/2+0.5])
% yticklabels({'k_{m,G6P}','k_{m,UDPGLC}','k_{cat}','k_{m,Pi}'}) % ,'k_{m,F6P}' taken out (1:4)
yticklabels({'k_{m,G6P}','k_{m,UDPGLC}','v_{max}','k_{m,Pi}'}) % ,'k_{m,F6P}' taken out (1:4)

% % leg1 = legend('w_{G6P}','w_{TRE}','w_{G6P+TRE}');
leg1 = legend('w_{G6P}','w_{G6P+TRE}');
% leg1.Postion = [0.2579 0.6677 0.0756 0.2178];
legend('boxoff')
legend('orientation','horizontal')
% legend('location','north')
legend('location','southeast')
box on
% 
temp_xloc = dp1.XLim(1) + (dp1.XLim(2) - dp1.XLim(1)) * 0.1;
temp_xloc2 = dp1.XLim(1) + (dp1.XLim(2) - dp1.XLim(1)) * 0.8;
temp_yloc = dp1.YLim(1) + (dp1.YLim(2) - dp1.YLim(1)) * 0.9;
temp_txt = text(temp_xloc, temp_yloc, 'C','FontSize',12);
temp_txt2 = text(temp_xloc2, temp_yloc, 'TPS1','FontSize',12);

line(tempXlim,[0+1/3 0+1/3],'LineStyle',':','color','k')
line(tempXlim,[1+1/3 1+1/3],'LineStyle',':','color','k')
line(tempXlim,[2+1/3 2+1/3],'LineStyle',':','color','k')
line(tempXlim,[3+1/3 3+1/3],'LineStyle',':','color','k')
line(tempXlim,[4+1/3 4+1/3],'LineStyle',':','color','k')
% % leg1.String(8) = [];
leg1.String(7) = [];
leg1.String(6) = [];
leg1.String(5) = [];
leg1.String(4) = [];
leg1.String(3) = [];
% 
set(gca,'XTickLabel',[]);
% 
leg1.FontSize = 11;
dp1.FontSize = 11;
hold off
    % add lines horizontal -- separating
    % style errobar plot

% %% 
sp2 = subplot(2,2,2);
% 
% % errorbar(pars_cell_selected{4},[1:2],errbar_cell_selected{4},...
% %     'horizontal','o','MarkerFaceColor',[0 0.4470 0.7410],...
% %     'color',[0 0.4470 0.7410],'MarkerEdgeColor',[0 0.4470 0.7410],...
% %     'CapSize',0,'LineWidth',1)
% % hold on
% % errorbar(pars_cell_selected{5},[1:2]-1/3,errbar_cell_selected{5},...
% %     'horizontal','o','MarkerFaceColor',[0.8500    0.3250    0.0980],...
% %     'color',[0.8500    0.3250    0.0980],...
% %     'MarkerEdgeColor',[0.8500    0.3250    0.0980],'CapSize',0,'LineWidth',1)
% 
errorbar(pars_cell_selected{4},[1:2],errbar_cell_selected{4},...
    'horizontal','o','MarkerFaceColor',c_grayish,...
    'color',c_grayish,'MarkerEdgeColor',c_grayish,...
    'CapSize',0,'LineWidth',1)
hold on
errorbar(pars_cell_selected{5},[1:2]-1/3,errbar_cell_selected{5},...
    'horizontal','o','MarkerFaceColor','k',...
    'color','k',...
    'MarkerEdgeColor','k','CapSize',0,'LineWidth',1)
% 
tempXlim = [-20 20];
xlim(tempXlim)
ylim([-1.5 4.5])
% ylim([-2.5 3.5])
% yticks([-1 0 1 2 3 4])
yticks([1/2+1/3 3/2+1/3])
yticklabels({'v_{max}','k_{m,GLC}'}) % ,'k_{m,F6P}' taken out (1:4)

leg2 = legend('w_{HXK}','w_{HXK+G6P}');
legend('boxoff')
legend('orientation','horizontal')
% legend('location','north')
legend('location','southeast')

% title('\color{blue}GLT. Combining fluxomics and metabolomics data')
box on
% 
temp_xloc = sp2.XLim(1) + (sp2.XLim(2) - sp2.XLim(1)) * 0.1;
temp_xloc2 = sp2.XLim(1) + (sp2.XLim(2) - sp2.XLim(1)) * 0.8;
temp_yloc = sp2.YLim(1) + (sp2.YLim(2) - sp2.YLim(1)) * 0.9;
temp_txt = text(temp_xloc, temp_yloc, 'D','FontSize',12);
temp_txt2 = text(temp_xloc2, temp_yloc, 'HXT','FontSize',12);

line(tempXlim,[0+1/3 0+1/3],'LineStyle',':','color','k')
line(tempXlim,[1+1/3 1+1/3],'LineStyle',':','color','k')
line(tempXlim,[2+1/3 2+1/3],'LineStyle',':','color','k')
leg2.String(5) = [];
leg2.String(4) = [];
leg2.String(3) = [];
set(gca,'XTickLabel',[]);
hold off
% 
leg2.FontSize = 11;
sp2.FontSize = 11;
hold off


% %% 
sp3 = subplot(2,2,3);

% % errorbar(pars_cell_selected{6},[1:7],errbar_cell_selected{6},'horizontal','o')
% errorbar(pars_cell_selected{6}(4:7),[1:4],errbar_cell_selected{6}(4:7),'horizontal','o','MarkerFaceColor',[0 0.4470 0.7410],'color',[0 0.4470 0.7410],'MarkerEdgeColor',[0 0.4470 0.7410],'CapSize',0,'LineWidth',1)
% hold on
% % errorbar(pars_cell_selected{7},[1:7]-1/3,errbar_cell_selected{7},'horizontal','o')
% errorbar(pars_cell_selected{7}(4:7),[1:4]-1/3,errbar_cell_selected{7}(4:7),'horizontal','o','MarkerFaceColor',[0.8500    0.3250    0.0980],'color',[0.8500    0.3250    0.0980],'MarkerEdgeColor',[0.8500    0.3250    0.0980],'CapSize',0,'LineWidth',1)
errorbar(pars_cell_selected{6}(4:7),[1:4],errbar_cell_selected{6}(4:7),...
    'horizontal','o','MarkerFaceColor',c_grayish,...
    'color',c_grayish,'MarkerEdgeColor',c_grayish,...
    'CapSize',0,'LineWidth',1)
hold on
errorbar(pars_cell_selected{7}(4:7),[1:4]-1/3,errbar_cell_selected{7}(4:7),...
    'horizontal','o','MarkerFaceColor','k',...
    'color','k',...
    'MarkerEdgeColor','k','CapSize',0,'LineWidth',1)
% 
tempXlim = [-200 200];
xlim(tempXlim)
% % ylim([-0.25 5.75])
ylim([-1 6])

% yticks([0 1 2 3 4 5])
yticks([1/3+1/2 1/3+1/2+1 1/3+1/2+2 1/3+1/2+3])
yticklabels({'k_{m,ADP}','k_{m,ATP}','k_{eq}','k_{m,G6P}'}) % ,'k_{m,F6P}' taken out (1:4)

leg3 = legend('w_{SS}','w_{SS+GP}');
legend('boxoff')
legend('orientation','horizontal')
% legend('location','north')
legend('location','southeast')

% title('HXK. Combining steady state and glucose perturbation data')
box on
% 
temp_xloc = sp3.XLim(1) + (sp3.XLim(2) - sp3.XLim(1)) * 0.1;
temp_xloc2 = sp3.XLim(1) + (sp3.XLim(2) - sp3.XLim(1)) * 0.8;
temp_yloc = sp3.YLim(1) + (sp3.YLim(2) - sp3.YLim(1)) * 0.9;
temp_txt = text(temp_xloc, temp_yloc, 'E','FontSize',12);
temp_txt2 = text(temp_xloc2, temp_yloc, 'HXK','FontSize',12);
% 
line(tempXlim,[0+1/3 0+1/3],'LineStyle',':','color','k')
line(tempXlim,[1+1/3 1+1/3],'LineStyle',':','color','k')
line(tempXlim,[2+1/3 2+1/3],'LineStyle',':','color','k')
line(tempXlim,[3+1/3 3+1/3],'LineStyle',':','color','k')
line(tempXlim,[4+1/3 4+1/3],'LineStyle',':','color','k')
leg3.String(7) = [];
leg3.String(6) = [];
leg3.String(5) = [];
leg3.String(4) = [];
leg3.String(3) = [];
set(gca,'XTickLabel',[]);
hold off
% 
leg3.FontSize = 11;
sp3.FontSize = 11;
hold off

%%
% % % % print(gcf,'-dpdf','E:\OneDrive - TU Eindhoven\Documents\ch4_Y3M1\results\figures\manuscript\F4CDE');

% % saving the plots
% savefig(1001,'results\F4_proteomics.fig')
% savefig(6002,'results\F4_estimation.fig')






% % % % clc, clear, close all
% % % % set_paths
% % % % dbstop if error
% % % % colorpalette;
% % % % savefig_loc = 'E:\OneDrive - TU Eindhoven\Documents\ch4_Y3M1\results\figures\manuscript\';
% % % % setup.save = 0;
% % % % % savefig_loc = 'E:OneDrive - TU Eindhoven\Documents\ch4_Y3M1\results\figures\manuscript';
% % % % 
% % % % % inital simulations
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
% % % % % load('x108a.mat','x_1st','x_2nd'); x_2nd(35) = x_2nd(38);
% % % % % x = x_2nd; 
% % % % setup.GLT_km_same_extra = 1;
% % % % setup.VmGLT_mu_dependent = 1;
% % % % load('x108c.mat','x_3rd');
% % % % x = x_3rd;



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

% plot simulations
setup.plotResultsMode = 30;
namesHits = cell(1,2); namesHits{1} = 'reference'; namesHits{2} = 'PvHoek off';
simRes = [refSim, PvHoek_off_Sim];
%
% [~] = plotAll_Y3M1(simRes,legenda,canelas_SS,data,setup,[x_2nd; x_2nd]);
[~] = plotAll_Y3M1(simRes,legenda,canelas_SS,data,setup,[x_3rd; x_3rd]);
setup.plotResultsMode = 0;

%% plot manuscript
for plotting_proteomics = 1
    % checkup
    if exist('f2_h','var')
        clf(f2_h)
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

    % figure
    f1001_h = figure(1001);
    f1001_h.Position = [1978 275 1227 710];

    clf(1001)

    % simple visualization
    sp1 = subplot(2,3,1);
    pt1 = plot(mu,HXK_sim_protOFF,'o:');
    pt1.MarkerSize = 6.5;
    pt1.Color = 'k';
    pt1.MarkerFaceColor = 'r';
    hold on
    % 
    pt2 = plot(mu,HXK_sim,'o:');
    pt2.MarkerSize = 6.5;
    pt2.Color = 'k';
    pt2.MarkerFaceColor = 'b';
    % 
    pt3 = plot(mu,HXK_exp,'k.-','LineWidth',1,'MarkerSize',10);
    % 
    xlabel(Y3M1_labels.mu)
    ylabel(Y3M1_labels.q_mMs)
    title('simple visualization')
    legend('simulation, no proteomics','simulation proteomics',...
        'experimental data','location','northwest')
    box on
    hold off

    % simulated vs experimental
    sp2 = subplot(2,3,2);
    sc1 = scatter(HXK_exp,HXK_sim_protOFF);
    sc1.LineWidth = 0.6;
    sc1.MarkerEdgeColor = 'k';
    sc1.MarkerFaceColor = 'r';
    hold on
    % 
    sc2 = scatter(HXK_exp,HXK_sim);
    sc2.LineWidth = 0.6;
    sc2.MarkerEdgeColor = 'k';
    sc2.MarkerFaceColor = 'b';
    % 
    line([0 2],[0 2],'color','k','LineStyle','--')
    % 
    xlabel('experimental reaction rate (mmol L^{-1} s^{-1})')
    ylabel('simulated reaction rate (mmol L^{-1} s^{-1})')
    title('experimental vs simulation plot')
    legend('simulation, no proteomics','simulation proteomics',...
        'location','southeast')
    box on

    % error plots 1
    c = categorical({'error','error norm.'});
    sp3 = subplot(2,3,4);
    bp1 = bar(c,[[HXK_sim_error, HXK_sim_protOFF_error];[sum(abs(HXK_sim_error_full)), sum(abs(HXK_sim_protOFF_error_full))]]','FaceColor','flat');
    bp1(1).CData(1,:) = [0 0 .7];
    bp1(1).CData(2,:) = [1 0 0];
    bp1(2).CData(1,:) = [0 0 .7];
    bp1(2).CData(2,:) = [1 0 0];
    % legend('error sim NO prot','error sim','location','northwest')
    ylabel('error')
    title('error plots')

    % error plots 2
    % c2 = categorical({'error sim NO prot','error sim'});
    c2 = categorical({'sim + Prot','sim - prot'});
    sp4 = subplot(2,3,5);
    bp2 = bar(c2,[sum(abs(HXK_sim_error_full)), sum(abs(HXK_sim_protOFF_error_full))],'FaceColor','flat');
    bp2.CData(1,:) = [0 0 .7];
    bp2.CData(2,:) = [1 0 0];
    % legend('error sim NO prot','error sim','location','northwest')
    ylabel('error')
    title('error plots, only part')

    % simulated vs experimental + error plots
    a1 = subplot(2,3,3);
    copyobj(get(sp2,'children'),a1);
    % 
    xlabel('experimental reaction rate (mmol L^{-1} s^{-1})')
    ylabel('simulated reaction rate (mmol L^{-1} s^{-1})')
    title('proposed final plot','color','red')
    legend('simulation, no proteomics','simulation proteomics',...
        'location','southeast')
    box on
    % addition
    handaxes1 = axes('Position', [0.695 0.8125 a1.Position(3)/3.2 a1.Position(4)/3.2]);
    copyobj(get(sp4,'children'),handaxes1);
    handaxes1.XLim = [0 3];
    handaxes1.YLim = [0 3];
    set(gca,'XTickLabel',[],'YTickLabel',[]);
    box on
    txt1 = text(0.2,2.5,'error','FontSize',13);

    % 
    set(f1001_h,'color','w')

    % 
    if setup.save == 1
        saveName = [savefig_loc, 'fig7a_proteomics.fig'];
        savefig(f1001_h, saveName)
    %     % 
    %     set(gcf,'Units','inches');
    %     screenposition = get(3,'Position');
    %     set(gcf,...
    %         'PaperPosition',[0 0 screenposition(3:4)],...
    %         'PaperSize',[screenposition(3:4)]);
    %     print -dpdf -painters Y3M1_robustness_rates
    %     print -dpng -painters Y3M1_robustness_rates    
    end
end


%% (2) Comparison parameter fits SS AND/OR GP: fit parameters + Error + CI.
% recall setup
% select parameters (pars HXK)
% For (2)
% 3 cost functions (SS AND/OR GP) + control (where nothing is optimized for)
    % first onlt fit exp SS rates
    % for the other exp GP rates as well
    % (if it would always fit good, then add noise at start, or back to literature values, zeros)
% run in PC uni, fake parallel

% For (3)
    % fitting to rates HXK in one and the other
    % fitting also to tre data
    % (if it would always fit good, then add noise at start, or back to literature values, zeros)

% For (4)
    % within GP
    % First concentrations
    % Then rates
    % (if it would always fit good, then add noise at start, or back to literature values, zeros)


% blank and constant setup
setup.plotResultsMode = 0;
blankWeight = zeros(1,140);
lambdalist = 0;
setup.parEst.lambda = lambdalist(1); lam = setup.parEst.lambda;
options = optimoptions('lsqnonlin','Display','iter',...
    'FunctionTolerance', 1e-4, 'OptimalityTolerance', 1e-4);
setup.parEst.costfun = 10015;
% 
ntests = 6;

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


% % % % % loop estimation
% % % % cluster = parcluster('local');
% % % % pool = parpool(cluster,ntests);
% % % % parfor j = 1:ntests
% % % % % % % % for j = 1
% % % % 
% % % %     % select lambda
% % % %     lam = 0;
% % % % 
% % % %     % select specifics each round
% % % %     selPars = parComb_cell{j};
% % % %     w = warray_cell{j};
% % % %     saveName = names_cell{j};
% % % %     % adjustment after selection
% % % %     plength = length(selPars);
% % % %     lb = -3*ones(1,plength);
% % % %     ub = 3*ones(1,plength);
% % % %     x_temp = x(selPars);
% % % % 
% % % % % % % %     % test run
% % % % % % % %     [testRun_error]=costfunSystemY3M1_temp_separateParsGPSS3(x_temp,canelas_SS,setup,x,data,dataset,w,lam,selPars);
% % % % 
% % % %     % estimation
% % % %     sprintf('enzyme = %d, lambda = %d', 1, j)
% % % %     tic
% % % %     [xres,resnorm,residual,exitflag,output,tempLam,jacobian] = lsqnonlin(@costfunSystemY3M1_temp_separateParsGPSS3,x_temp,lb,ub,options,canelas_SS,setup,x,data,dataset,w,lam,selPars);
% % % % % % % %     xres = 1; resnorm = 1; residual = 1; exitflag = 1;
% % % %     t = toc; 
% % % %     disp(xres);
% % % %     
% % % %     % simulation
% % % %     xres_full = x;
% % % %     xres_full(selPars) = xres;
% % % %     [simRes_temp] = simulateY3M1_separateGPSS(xres_full, canelas_SS, data, dataset, setup);
% % % % 
% % % %     % error simulation
% % % %     x2 = x; x2(selPars) = xres;
% % % %     lam2 = 0;
% % % %     [error]=costfunSystemY3M1_temp_separateParsGPSS3(xres,canelas_SS,setup,x2,data,dataset,w,lam2,selPars);
% % % %     errorSim = sum(abs(error));
% % % %     % error parameters
% % % %     errorLam = sum(abs(xres));
% % % % 
% % % %     % saving result
% % % % % %         parsave_Y3M2_cluster_3(saveName,xres,resnorm,residual,exitflag,t,w,lam,simRes_temp)
% % % % %     parsave_Y3M2_cluster_4(saveName,xres,resnorm,residual,exitflag,t,w,lam,simRes_temp, error, errorSim, errorLam)
% % % %     parsave_Y3M2_cluster_5(saveName,xres,resnorm,residual,exitflag,t,w,lam,simRes_temp, error, errorSim, errorLam,output,tempLam,jacobian)
% % % % 
% % % % end
% % % % delete(pool)


%% DATA RECALL (not the multistart cases yet)
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


%% visualization 3: actualy simulations
xDummy = [x; ...
            x; x; x; ...
            x; x; x; ...
            ];
% simRes{1} = simRes{1}{1};
% % % % for i = 1:nEnzs
% % % %     simRes{i,1}{1} = simRes{i,1};
% % % % end
% % % % From here on it will probably change


%% 110A
% idxsSel = 1; % just reference simulation
% idxsSel = [1:17]; %
idxsSel_test1 = [1,2:4];
idxsSel_test2 = [1,5:7];

%% TPS1
close all
idxsSel_temp = idxsSel_test1;
setup.plotResultsMode = 30;
setup.namesHits = namesHits(idxsSel_temp);
[~] = plotAll_Y3M1(simRes(idxsSel_temp),legenda,canelas_SS,data,setup,xDummy(idxsSel_temp,:));
% looks as expected

%% GLT
close all
idxsSel_temp = idxsSel_test2;
setup.plotResultsMode = 30;
setup.namesHits = namesHits(idxsSel_temp);
[~] = plotAll_Y3M1(simRes(idxsSel_temp),legenda,canelas_SS,data,setup,xDummy(idxsSel_temp,:));
% looks as expected

%% confidence intervals
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
ci_cell = cell(1,nS);
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
    ci_cell{i} = nlparci(parameters,residual,'Jacobian',Jacobian);
end

%% nice trend for TPS1
abs([ci_cell{1}(:,1) - ci_cell{1}(:,2),...
    ci_cell{2}(:,1) - ci_cell{2}(:,2),...
    ci_cell{3}(:,1) - ci_cell{3}(:,2)])

%% can be only narrated form rate to rate+met
abs([ci_cell{4}(:,1) - ci_cell{4}(:,2),...
    ci_cell{5}(:,1) - ci_cell{5}(:,2),...
    ci_cell{6}(:,1) - ci_cell{6}(:,2)])


% % Narrative: can use the first story with HXK pretty well. Now for the
% others it does not look too good. Test in a second run with TPS1 the
% second case and GLT the third.

% % from
% warray_cell{2}(w2ss) = ones;
% % to
% warray_cell{4}([w2ss, w2gp]) = ones;
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



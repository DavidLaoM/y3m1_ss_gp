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
% % % % 
% % % % % (2022 01 21) Latest adjustments for individual GLT
% % % % setup.sameGLTHXK4all = 1;
% % % % load('x108a.mat','x_1st','x_2nd'), x_2nd_prev = x_2nd;
% % % % setup.GLT_km_same_extra = 1;
% % % % load('x108b.mat','x_1st','x_2nd')
% % % % setup.VmGLT_mu_dependent = 1;
% % % % % % % % x = x_2nd;
% % % % x = x_2nd_prev;
% % % % load('x108c.mat'),
% % % % x = x_3rd;
% % % % 
% % % % 
% % % % % %% Simulations
% % % % % Y3M1manus_labels;
% % % % % close all % start cleaning up
% % % % % % 
% % % % % setup.GLT_km_same = 1;
% % % % % setup.plotResultsMode = 30;
% % % % % [simRes_reference] = simulateY3M1_separateGPSS(x, canelas_SS, data, dataset, setup);
% % % % % setup.simGPdata = simRes_reference{1}.gs;
% % % % % setup.simSSdata = simRes_reference{1}.ss;
% % % % % setup.litParams = zeros(size(x));
% % % % % setup.plotResultsMode = 0;




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




%% Test simulation
setup.runSScanelas = 1;
setup.runGSvanHeerden = 1;
setup.plotBackColor = [1 1 1];
setup.plotResultsMode = 0;
simRes = simulateY3M1_separateGPSS(x, canelas_SS, data, dataset, setup);
setup.plotResultsMode = 0;


%% visualization
% clf(5001)
% clf(5002)
% clf(5003)
% clf(5004)
% clf(5005)

% % SS concentrations
% select metabolites_SS to plot
setup.metabolitesSS_2plot = [...
    5 4 3 ... % g6p f6p fbp
    19 18 17 14 ... % gly g3p gap dhap
    11 10 12 13 20 ... % p3g p2g pep pyr etoh
    9 15 16]; % atp adp amp
% select the subplot location of metabolites ss in figure
setup.metabolitesSS_location = [...
    10 11 12 ... %3 7 11 ... % g6p f6p fbp
    23 22 13 4 ... % 13 14 15 16 ... % gly g3p gap dhap
    5 6 7 8 9 ... %20 24 28 32 36 ... % p3g p2g pep pyr etoh
    25 26 27]; %22 26 30]; % atp adp amp

% % SS rates
% select metabolites_SS to plot
setup.ratesSS_2plot = [...
    1 2 3 4 5 7 6 ... % glt hxk pgi pfk ald tpi gpd
    8 9 10 11 12 13 41 ... % tdh pgk pgm eno pyk pdc adh
    34 35 36 37 38 39 40]; % sinkG6P, f6p, gap, p3g, pep, pyr, acer
% select the subplot location of metabolites ss in figure
setup.ratesSS_location = 1:21; %[...
%     1 2 3 4 5 6 7 ... % glt hxk pgi pfk ald tpi gpd
%     8 9 10 11 12 13 14 ... % tdh pgk pgm eno pyk pdc adh
%     % 34 35 36 37 38 39 40]; % sinkG6P, f6p, gap, p3g, pep, pyr, acer

% % gp concentrations
% select metabolites_SS to plot
setup.metabolitesGP_2plot = [...
    5 4 3 14 11 12 ... % g6p f6p fbp gap p2g pep
    21 26 25 9 15 16 ... % g1p t6p tre atp adp amp
    6 27]; % glci pi
% select the subplot location of metabolites ss in figure
setup.metabolitesGP_location = 1:14;

% % gp concentrations
% select metabolites_SS to plot
setup.ratesGP_2plot = [...
    2 3 4 5 7 ... % glk pgi pfk ald tdh
    17 21 19 20]; % pgm1 tps1 tps2 nth1
% select the subplot location of metabolites ss in figure
setup.ratesGP_location = 1:19;

setup.plotResultsMode = 30;
[~] = plotAll_Y3M1_manuscriptFig5(simRes,legenda,canelas_SS,data,setup,x);
setup.plotResultsMode = 0;

%% after treatment
fh_5001 = figure(5001);
fh_5002 = figure(5002);
fh_5004 = figure(5004);
fh_5005 = figure(5005);
fh_5001.Children(9).Children(1).String = 'GAP,';
fh_5001.Children(13).Children(1).String = 'FBP,';
fh_5001.Children(11).Children(1).String = 'G3P,';
fh_5001.Children(12).Children(1).String = 'GLYC,';
% remove arrows manually (in 5001)
% set names location correctly (in fig 5004)


%% saving
% fh_5001
set(fh_5001,'Units','Inches');
pos = get(fh_5001,'Position');
set(fh_5001,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(fh_5001,'E:\OneDrive - TU Eindhoven\Documents\ch4_Y3M1\results\figures\manuscript\SM_simulations_SS_conc','-dpdf','-r0')
% fh_5002
set(fh_5002,'Units','Inches');
pos = get(fh_5002,'Position');
set(fh_5002,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(fh_5002,'E:\OneDrive - TU Eindhoven\Documents\ch4_Y3M1\results\figures\manuscript\SM_simulations_SS_rate','-dpdf','-r0')
% fh_5004
set(fh_5004,'Units','Inches');
pos = get(fh_5004,'Position');
set(fh_5004,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(fh_5004,'E:\OneDrive - TU Eindhoven\Documents\ch4_Y3M1\results\figures\manuscript\SM_simulations_GP_conc','-dpdf','-r0')
% fh_5005
set(fh_5005,'Units','Inches');
pos = get(fh_5005,'Position');
set(fh_5005,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(fh_5005,'E:\OneDrive - TU Eindhoven\Documents\ch4_Y3M1\results\figures\manuscript\SM_simulations_GP_rate','-dpdf','-r0')





%%
function [numFigs] = plotAll_Y3M1_manuscriptFig5(simResults,legenda,canelas_SS,data,setup,xAll)

% recall data
[npSets,~] = size(xAll);
colorSet = cool(npSets);
recallSimsData_Y3M1;
% default reference parameter array
refParams = exist('setup.refParams');
if refParams == 0 % if 0, first array gets plotted in bold. If-1, none does
    setup.refParams = 1;
end
% background prot color (if exists. If not, default)
if isfield(setup,'plotBackColor')
    backColor = setup.plotBackColor;
else
    backColor = [0.94 0.94 0.94];
end


% plots profiles
if((setup.plotResultsMode == 1)||(setup.plotResultsMode == 2)||(setup.plotResultsMode == 10)||(setup.plotResultsMode == 30))

    
    % fig1. Steady state metabolite profile. All
    for simSSconc = 1
    if setup.runSScanelas == 1
%         figure('name','Steady state metabolite profile. All','units','normalized','outerposition',[0 0 1 1])
        fh_5001 = figure(5001);
        fh_5001.Position = [1969 402 1099 551];
        
% % % %         for k = 1:length(legenda.metabolites)
        for k_pre = 1:length(setup.metabolitesSS_2plot)
            k = setup.metabolitesSS_2plot(k_pre);
            k_loc = setup.metabolitesSS_location(k_pre);
% % % %             subplot(8,4,k)
% % % %             subplot(9,4,k_loc)
            sp = subplot(3,9,k_loc);
            % plot simulations
            for j = 1:npSets
                tempVal = zeros(1,length(dprofile{j}));
                for o = 1:length(dprofile{j})
                    tempVal(o) = Yss{j}.ss2can{o}(end,k);
                end
                if j == setup.refParams
                    % interp the pchip
                    darr = dprofile{1}(1):0.005:dprofile{1}(8);
                    parr = pchip(dprofile{j},tempVal,darr);
                    % plot the value
                    plot(darr, parr, 'k-','LineWidth',1)
                    hold on
                    % plot the points
                    plot(dprofile{j}, tempVal, 'k.','MarkerSize',5)
                else
                    plot(dprofile{j}, tempVal, 'color', colorSet(j,:))
                end
                hold on
            end
            % plot experimental data
            if canelas_SS.ordered.metabolites.onoff(k) == 1
                errorbar(canelas_SS.ordered.dprof{k},canelas_SS.ordered.metabolites.data{k},canelas_SS.ordered.metabolites.std{k},...
                    's','MarkerEdgeColor','k','MarkerFaceColor','k',...
                    'MarkerSize',3,'color','k')%,...
                    %'CapSize ',5)
% % % %                 errorbar(canelas_SS.ordered.dprof{k},canelas_SS.ordered.metabolites.data{k},canelas_SS.ordered.metabolites.std{k},'ks')
            end
            % arrows issues
            if k_loc == 9
                tempText = text(0.03,-sp.YLim(2)/10,'\uparrow');
                tempText.FontSize = 20;
                tempText = text(0.03+0.05,-sp.YLim(2)/6,'1');
                tempText.FontSize = 10;
            elseif k_loc == 25
                tempText = text(0.03,2.35,'\uparrow');
                tempText.FontSize = 20;
                tempText = text(0.03+0.05,2.25,'1');
                tempText.FontSize = 10;
            elseif k_loc == 26
                tempText = text(0.03,-sp.YLim(2)/10,'\uparrow');
                tempText.FontSize = 20;
                tempText = text(0.03+0.05,-sp.YLim(2)/6,'1');
                tempText.FontSize = 10;
            elseif k_loc == 27
                tempText = text(0.03,-sp.YLim(2)/10,'\uparrow');
                tempText.FontSize = 20;
                tempText = text(0.03+0.05,-sp.YLim(2)/6,'1');
                tempText.FontSize = 10;
            elseif k_loc == 10
                tempText = text(0.18,-sp.YLim(2)/10,'\uparrow');
                tempText.FontSize = 20;
                tempText = text(0.23,-sp.YLim(2)/6,'2');
                tempText.FontSize = 10;
            elseif k_loc == 11
                tempText = text(0.18,-sp.YLim(2)/10,'\uparrow');
                tempText.FontSize = 20;
                tempText = text(0.23,-sp.YLim(2)/6,'2');
                tempText.FontSize = 10;
            end
            
            set(gca,'XTickLabel',[]);
%             set(gca,'YTickLabel',[]);
            titleName = erase(legenda.metabolites{k},", [mM]");
            titleName = titleName(1:end-7);
            text(sp.XLim(2)*0.1,sp.YLim(2)*0.9,titleName,'FontSize',8)
            yticks(sp.YLim)
%             temp = title(titleName);%fontsize{16}
%             ylabel(titleName)%fontsize{16}
            if k == 25
%                 xlabel('dilution rate [h^{-1}]')
%                 ylabel('concentration [mM]')
            end
        end
        
%         % added plot NAD:NADH ratio
%         sp = subplot(3,9,19);
%         % recall data
%         tempVal = zeros(1,length(dprofile{j}));
%         for o = 1:length(dprofile{1})
%             tempVal(o) = Yss{1}.ss2can{o}(end,7)./Yss{1}.ss2can{o}(end,8);
%         end
%         % interp the pchip
%         darr = dprofile{1}(1):0.005:dprofile{1}(8);
%         parr = pchip(dprofile{j},tempVal,darr);
%         % plot the value
%         plot(darr, parr, 'k-','LineWidth',1)
%         hold on
%         % plot the points
%         plot(dprofile{j}, tempVal, 'k.','MarkerSize',5)
%         % 
%         plot(canelas_SS.mtlD.D,canelas_SS.mtlD.NAD_NADHratio,...
%             'ks','MarkerEdgeColor','k','MarkerFaceColor','k',...
%             'MarkerSize',3)%,...

        
%         suptitle('fig1. Steady state metabolite profile. All')
        % set(gcf,'color','w');
% % % %         set(gcf,'color',backColor);
    end
    end
    

    % fig2. Steady state flux profile. All
    for simSSrates = 1
    if setup.runSScanelas == 1
%         figure('name','Steady state flux profile. All','units','normalized','outerposition',[0 0 1 1])
%         figure(5002)
        fh_5002 = figure(5002);
        fh_5002.Position = [1969 404 1114 550];
        
% % % %         for k = 1:length(legenda.fluxes)
        for k_pre = 1:length(setup.ratesSS_2plot)
            k = setup.ratesSS_2plot(k_pre);
            k_loc = setup.ratesSS_location(k_pre);
            
            sp = subplot(3,7,k_loc);
            
% % % %             subplot(5,9,k)
            % plot simulations
            for j = 1:npSets
                tempVal = zeros(1,length(dprofile{j}));
                for o = 1:length(dprofile{j})
                    tempVal(o) = Vss{j}.ss2can{o}(end,k);
                end
                if j == setup.refParams
                    % interp the pchip
                    darr = dprofile{1}(1):0.005:dprofile{1}(8);
                    parr = pchip(dprofile{j},tempVal,darr);
                    % plot the value
                    plot(darr, parr, 'k-','LineWidth',1)
                    hold on
                    % plot the points
                    plot(dprofile{j}, tempVal, 'k.','MarkerSize',5)
%                     plot(dprofile{j}, tempVal, 'k','LineWidth',1.2)
                else
                    plot(dprofile{j}, tempVal, 'color', colorSet(j,:))
                end
                hold on
            end
            % plot experimental data
            if canelas_SS.ordered.fluxes.onoff(k) == 1
                plot(canelas_SS.ordered.dprof{end},canelas_SS.ordered.fluxes.data{k},...
                    'kd',...
                    'MarkerFaceColor','k','MarkerSize',3)
%                 errorbar(canelas_SS.ordered.dprof{k},canelas_SS.ordered.metabolites.data{k},canelas_SS.ordered.metabolites.std{k},...
%                     's','MarkerEdgeColor','k','MarkerFaceColor','k',...
%                     'MarkerSize',3,'color','k')%,...
            end
            % arrows issues
            if k_loc == 1
                tempText = text(0.03,-sp.YLim(2)/10,'\uparrow');
                tempText.FontSize = 20;
                tempText = text(0.03+0.05,-sp.YLim(2)/6,'1');
                tempText.FontSize = 10;
            elseif k_loc == 14
                tempText = text(0.03,-sp.YLim(2)/10,'\uparrow');
                tempText.FontSize = 20;
                tempText = text(0.03+0.05,-sp.YLim(2)/6,'1');
                tempText.FontSize = 10;
            end
                        
%             set(gca,'XTickLabel',[]);
%             set(gca,'YTickLabel',[]);
            titleName = erase(legenda.fluxes{k},", [mM]");
            if((k_loc == 6)||(k_loc == 8))
                titleName = titleName(1:8);
            elseif k_loc == 7
                titleName = titleName(1:9);
            elseif k_loc <= 14
                titleName = titleName(1:7);
            else
                titleName = titleName(1:11);
            end
            text(sp.XLim(2)*0.1,sp.YLim(2)*0.9,titleName,'FontSize',8)
            yticks(sp.YLim)
            
        end
       
        suptitle('fig2. Steady state flux profile. All')
        % set(gcf,'color','w');
% % % %         set(gcf,'color',backColor);
    end
    end


    % fig21. Physiology. Steady state profile. All
    for simSSphysioVal = 1
    if setup.runSScanelas == 1
%         figure('name','Physiology. Steady state profile. All','units','normalized','outerposition',[0 0 1 1])
        figure(5003)
%         fh_5003 = figure(5003);
%         fh_5003.Position = [1969 404 1114 550];
        % % Inputs
        % qS
        sp = subplot(3,4,1); %3,6,2)
        for j = 1:npSets
            if j == setup.refParams
                % interp the pchip
                darr = dprofile{1}(1):0.005:dprofile{1}(8);
                parr = pchip(dprofile{j},exp_v_GLT{j},darr);
                % plot the value
                plot(darr, parr, 'k-','LineWidth',1)
                hold on
                % plot the points
                plot(dprofile{j}(1:q(1)), exp_v_GLT{j}, 'k.','MarkerSize',5)
            else
                plot(dprofile{j}(1:q(1)), exp_v_GLT{j}, 'color', colorSet(j,:))
            end
        hold on
        %
        plot(canelas_SS.mtlD.D, canelas_SS.mtlD.v_GLT,...
            'kd',...
            'MarkerFaceColor','k','MarkerSize',3)
%         plot(canelas_SS.mtlD.D, canelas_SS.mtlD.v_GLT, 'r+')
        end
        % 
        set(gca,'XTickLabel',[]);
        titleName = 'q_{S}';
        text(sp.XLim(2)*0.1,sp.YLim(2)*0.9,titleName,'FontSize',8)
        yticks(sp.YLim)
        
        % % Exchange rates
        % qEtoh
        sp = subplot(3,4,5);%3,6,8)
        for j = 1:npSets
            if j == setup.refParams
                % interp the pchip
                darr = dprofile{1}(1):0.005:dprofile{1}(8);
                parr = pchip(dprofile{j},exp_v_ET_tr{j},darr);
                % plot the value
                plot(darr, parr, 'k-','LineWidth',1)
                hold on
                % plot the points
                plot(dprofile{j}(1:q(1)), exp_v_ET_tr{j}, 'k.','MarkerSize',5)
%                 plot(dprofile{j}(1:q(1)), exp_v_ET_tr{j}, 'k','LineWidth',1.2)
            else
                plot(dprofile{j}(1:q(1)), exp_v_ET_tr{j}, 'color', colorSet(j,:))
            end
        hold on
%         plot(canelas_SS.mtlD.D, canelas_SS.mtlD.v_ADH, 'r+')
        plot(canelas_SS.mtlD.D, canelas_SS.mtlD.v_ADH,...
            'kd',...
            'MarkerFaceColor','k','MarkerSize',3)
        end
%         title('q_{Etoh}')
        set(gca,'XTickLabel',[]);
        titleName = 'q_{Etoh}';
        text(sp.XLim(2)*0.1,sp.YLim(2)*0.9,titleName,'FontSize',8)
        yticks(sp.YLim)
            % arrows issues
            tempText = text(0.03,-sp.YLim(2)/10,'\uparrow');
            tempText.FontSize = 20;
            tempText = text(0.03+0.05,-sp.YLim(2)/6,'1');
            tempText.FontSize = 10;
        
        % qGlycerol
        sp = subplot(3,4,9);%3,6,9)
        for j = 1:npSets
            if j == setup.refParams
                darr = dprofile{1}(1):0.005:dprofile{1}(8);
                parr = pchip(dprofile{j},exp_v_GH_tr{j},darr);
                % plot the value
                plot(darr, parr, 'k-','LineWidth',1)
                hold on
                % plot the points
                plot(dprofile{j}(1:q(1)), exp_v_GH_tr{j}, 'k.','MarkerSize',5)
%                 plot(dprofile{j}(1:q(1)), exp_v_GH_tr{j}, 'k','LineWidth',1.2)
            else
                plot(dprofile{j}(1:q(1)), exp_v_GH_tr{j}, 'color', colorSet(j,:))
            end
        hold on
%         plot(canelas_SS.mtlD.D, canelas_SS.mtlD.v_G3PDH, 'r+')
        plot(canelas_SS.mtlD.D, canelas_SS.mtlD.v_G3PDH,...
            'kd',...
            'MarkerFaceColor','k','MarkerSize',3)
        end
%         title('q_{Glycerol}')
%         set(gca,'XTickLabel',[]);
        titleName = 'q_{Glycerol}';
        text(sp.XLim(2)*0.1,sp.YLim(2)*0.9,titleName,'FontSize',8)
        yticks(sp.YLim)
        
        % qCO2
        sp = subplot(3,4,6); %3,6,11)
        for j = 1:npSets
            if j == setup.refParams
                darr = dprofile{1}(1):0.005:dprofile{1}(8);
                parr = pchip(dprofile{j},(exp_v_sinkPYR{j} * 3 + 1 * exp_v_PDC{j})*1.7/2,darr);
                % plot the value
                plot(darr, parr, 'k-','LineWidth',1)
                hold on
                % plot the points
                plot(dprofile{j}(1:q(1)), (exp_v_sinkPYR{j} * 3 + 1 * exp_v_PDC{j})*1.7/2, 'k.','MarkerSize',5)

% %                 plot(dprofile{j}(1:q(1)), exp_v_sinkPYR{j} * 2 + 1 * exp_v_PDC{j}, 'k','LineWidth',1.2)
% %                 plot(dprofile{j}(1:q(1)), exp_v_sinkPYR{j} * 3 + 1 * exp_v_PDC{j}, 'k','LineWidth',1.2)
%                 plot(dprofile{j}(1:q(1)), (exp_v_sinkPYR{j} * 3 + 1 * exp_v_PDC{j})*1.7/2, 'k','LineWidth',1.2)
            else
%                 plot(dprofile{j}(1:q(1)), exp_v_sinkPYR{j} * 2 + 1 * exp_v_PDC{j}, 'color', colorSet(j,:))
%                 plot(dprofile{j}(1:q(1)), exp_v_sinkPYR{j} * 3 + 1 * exp_v_PDC{j}, 'color', colorSet(j,:))
                plot(dprofile{j}(1:q(1)), (exp_v_sinkPYR{j} * 3 + 1 * exp_v_PDC{j})*1.7/2, 'color', colorSet(j,:))
            end
        hold on
%         plot(canelas_SS.mtlD.D, canelas_SS.mtlD.qCO2, 'b+')
        plot(canelas_SS.mtlD.D, canelas_SS.mtlD.qCO2,...
            'kd',...
            'MarkerFaceColor','k','MarkerSize',3)
        end
%         title('q_{CO2}')
%         set(gca,'XTickLabel',[]);
        titleName = 'q_{CO2}';
        text(sp.XLim(2)*0.1,sp.YLim(2)*0.9,titleName,'FontSize',8)
        yticks(sp.YLim)
        
        % qO2
        sp = subplot(3,4,2); %3,6,12)
        for j = 1:npSets
            if j == setup.refParams
                darr = dprofile{1}(1):0.005:dprofile{1}(8);
                parr = pchip(dprofile{j}, (exp_v_sinkPYR{j} * 3 + 1 * exp_v_PDC{j})*1.7/2./[1.05 1.04 1.04 1.05 1.48 1.66 2.80 4.31]',darr);
                % plot the value
                plot(darr, parr, 'k-','LineWidth',1)
                hold on
                % plot the points
                plot(dprofile{j}(1:q(1)), (exp_v_sinkPYR{j} * 3 + 1 * exp_v_PDC{j})*1.7/2./[1.05 1.04 1.04 1.05 1.48 1.66 2.80 4.31]', 'k.','MarkerSize',5)

% %                 plot(dprofile{j}(1:q(1)), exp_v_mito{j}/1.1, 'k','LineWidth',1.2)
%                 plot(dprofile{j}(1:q(1)), (exp_v_sinkPYR{j} * 3 + 1 * exp_v_PDC{j})*1.7/2./[1.05 1.04 1.04 1.05 1.48 1.66 2.80 4.31]', 'k','LineWidth',1.2)
            else
%                 plot(dprofile{j}(1:q(1)), exp_v_mito{j}/1.1, 'color', colorSet(j,:))
                plot(dprofile{j}(1:q(1)), (exp_v_sinkPYR{j} * 3 + 1 * exp_v_PDC{j})*1.7/2./[1.05 1.04 1.04 1.05 1.48 1.66 2.80 4.31]', 'color', colorSet(j,:))
%                 plot(dprofile{j}(1:q(1)), (exp_v_sinkPYR{j} * 3 + 1 * exp_v_PDC{j})*1.7/2./[1.05 1.04 1.04 1.05 1.48 1.66 2.80 4.31], 'k','LineWidth',1.2)
            end
        hold on
%         plot(canelas_SS.mtlD.D,canelas_SS.mtlD.qO2, 'b+')
        plot(canelas_SS.mtlD.D, canelas_SS.mtlD.qO2,...
            'kd',...
            'MarkerFaceColor','k','MarkerSize',3)
        end
%         title('q_{O2}')
%         set(gca,'XTickLabel',[]);
        titleName = 'q_{O2}';
        text(sp.XLim(2)*0.1,sp.YLim(2)*0.9,titleName,'FontSize',8)
        yticks(sp.YLim)
        
        % calc qmito and qATPase
        for k = 27
            qATPase = zeros(1,length(dprofile{1}));
            qmito = zeros(1,length(dprofile{1}));
            for o = 1:length(dprofile{1})
                qATPase(o) = Vss{1}.ss2can{o}(end,k);
                qmito(o) = Vss{1}.ss2can{o}(end,k+1);
            end
        end
        % qmito
        sp = subplot(3,4,3); %3,6,5)
        for j = 1:npSets
            % calc qmito and qATPase
            qmito = zeros(1,length(dprofile{1}));
            for o = 1:length(dprofile{1})
                qmito(o) = Vss{j}.ss2can{o}(end,k+1);
            end            
            % sims            
            if j == setup.refParams
                darr = dprofile{1}(1):0.005:dprofile{1}(8);
                parr = pchip(dprofile{j},qmito,darr);
                % plot the value
                plot(darr, parr, 'k-','LineWidth',1)
                hold on
                % plot the points
                plot(dprofile{j}(1:q(1)), qmito, 'k.','MarkerSize',5)

%                 plot(dprofile{j}(1:q(1)), qmito, 'k','LineWidth',1.2)
            else
                plot(dprofile{j}(1:q(1)), qmito, 'color', colorSet(j,:))
            end
        hold on
%         plot(canelas_SS.mtlD.D,(exp_v_sinkPYR{j} * 3 + 1 * exp_v_PDC{j})*1.7/2./[1.05 1.04 1.04 1.05 1.48 1.66 2.80 4.31]'/1, 'b+')
%         plot(canelas_SS.mtlD.D,(canelas_SS.mtlD.sinkPYR * 3 + 1 * canelas_SS.mtlD.v_PDC)*1.7/2./[1.05 1.04 1.04 1.05 1.48 1.66 2.80 4.31]/1, 'b+')
        plot(canelas_SS.mtlD.D, (canelas_SS.mtlD.sinkPYR * 3 + 1 * canelas_SS.mtlD.v_PDC)*1.7/2./[1.05 1.04 1.04 1.05 1.48 1.66 2.80 4.31]/1,...
            'kd',...
            'MarkerFaceColor','k','MarkerSize',3)
        end
%         title('q_{mito}')
%         set(gca,'XTickLabel',[]);
        titleName = 'q_{mito}';
        text(sp.XLim(2)*0.1,sp.YLim(2)*0.9,titleName,'FontSize',8)
        yticks(sp.YLim)
        
        % qATPase
        sp = subplot(3,4,7); %3,6,6)
        for j = 1:npSets
% % % %                
            qATPase = zeros(1,length(dprofile{1}));
            for o = 1:length(dprofile{1})
                qATPase(o) = Vss{j}.ss2can{o}(end,k);
            end            
% % % %             
            if j == setup.refParams
                darr = dprofile{1}(1):0.005:dprofile{1}(8);
                parr = pchip(dprofile{j},qATPase,darr);
                % plot the value
                plot(darr, parr, 'k-','LineWidth',1)
                hold on
                % plot the points
                plot(dprofile{j}(1:q(1)), qATPase, 'k.','MarkerSize',5)

%                 plot(dprofile{j}(1:q(1)), qATPase, 'k','LineWidth',1.2)
            else
                plot(dprofile{j}(1:q(1)), qATPase, 'color', colorSet(j,:))
            end
            hold on
%             plot(canelas_SS.mtlD.D,(60*canelas_SS.mtlD.D+0.7)/1.7/3.6, 'b+')
%             plot(canelas_SS.mtlD.D,(60*canelas_SS.mtlD.D+0.7)/2.0/3.6, 'b+')
            plot(canelas_SS.mtlD.D, (60*canelas_SS.mtlD.D+0.7)/2.0/3.6,...
                'kd',...
                'MarkerFaceColor','k','MarkerSize',3)
            hold on
%             plot(canelas_SS.mtlD.D,(40*canelas_SS.mtlD.D+0.7)/2.0/3.6, 'g+')
            plot(canelas_SS.mtlD.D, (40*canelas_SS.mtlD.D+0.7)/2.0/3.6,...
                'kd',...
                'MarkerFaceColor','k','MarkerSize',3)
%             plot(canelas_SS.mtlD.D, canelas_SS.mtlD.v_ADH,...
%                 'kd',...
%                 'MarkerFaceColor','k','MarkerSize',3)
        end
%         title('q_{ATPase}')
%         set(gca,'XTickLabel',[]);
        titleName = 'q_{ATPase}';
        text(sp.XLim(2)*0.1,sp.YLim(2)*0.9,titleName,'FontSize',8)
        text(0.3, 3.75, 'gam60', 'FontSize', 8)
        text(0.3, 1.25, 'gam40', 'FontSize', 8)
        yticks(sp.YLim)
        
% % % %         
        % Yxs
        sp = subplot(3,4,4); %3,6,13)
        for j = 1:npSets
            if j == setup.refParams
                darr = dprofile{1}(1):0.005:dprofile{1}(8);
                parr = pchip(dprofile{j},dprofile{j}./exp_v_GLT{j}',darr);
                % plot the value
                plot(darr, parr, 'k-','LineWidth',1)
                hold on
                % plot the points
                plot(dprofile{j}(1:q(1)), dprofile{j}./exp_v_GLT{j}', 'k.','MarkerSize',5)

%                 plot(dprofile{j}, dprofile{j}./exp_v_GLT{j}', 'k','LineWidth',1.2)
            else
                plot(dprofile{j}, dprofile{j}./exp_v_GLT{j}', 'color', colorSet(j,:))
            end
        hold on
%         plot(canelas_SS.mtlD.D,(canelas_SS.mtlD.D./canelas_SS.mtlD.v_GLT), 'r+')
        plot(canelas_SS.mtlD.D, (canelas_SS.mtlD.D./canelas_SS.mtlD.v_GLT),...
            'kd',...
            'MarkerFaceColor','k','MarkerSize',3)
        end        
%         title('Y_{XS}')
        set(gca,'XTickLabel',[]);
        titleName = 'Y_{XS}';
        text(sp.XLim(2)*0.1,sp.YLim(2)*0.9,titleName,'FontSize',8)
        yticks(sp.YLim)
            % arrows issues
            tempText = text(0.03,0.15,'\uparrow');
            tempText.FontSize = 20;
            tempText = text(0.03+0.05,0.15,'1');
            tempText.FontSize = 10;
        
        % Yps
        sp = subplot(3,4,8); %3,6,14)
        for j = 1:npSets
            if j == setup.refParams
                darr = dprofile{1}(1):0.005:dprofile{1}(8);
                parr = pchip(dprofile{j}, exp_v_ET_tr{j}./exp_v_GLT{j}, darr);
                % plot the value
                plot(darr, parr, 'k-','LineWidth',1)
                hold on
                % plot the points
                plot(dprofile{j}(1:q(1)), exp_v_ET_tr{j}./exp_v_GLT{j}, 'k.','MarkerSize',5)

%                 plot(dprofile{j}, exp_v_ET_tr{j}./exp_v_GLT{j}, 'k','LineWidth',1.2)
            else
                plot(dprofile{j}, exp_v_ET_tr{j}./exp_v_GLT{j}, 'color', colorSet(j,:))
            end
        hold on
%         plot(canelas_SS.mtlD.D,(canelas_SS.mtlD.v_ADH./canelas_SS.mtlD.v_GLT), 'r+')
        plot(canelas_SS.mtlD.D, (canelas_SS.mtlD.v_ADH./canelas_SS.mtlD.v_GLT),...
            'kd',...
            'MarkerFaceColor','k','MarkerSize',3)
        end    
%         title('Y_{PS}')
%         set(gca,'XTickLabel',[]);
        titleName = 'Y_{PS}';
        text(sp.XLim(2)*0.1,sp.YLim(2)*0.9,titleName,'FontSize',8)
        yticks(sp.YLim)
            % arrows issues
            tempText = text(0.03,-sp.YLim(2)/10,'\uparrow');
            tempText.FontSize = 20;
            tempText = text(0.03+0.05,-sp.YLim(2)/6,'1');
            tempText.FontSize = 10;
        
        sp = subplot(3,4,12); %3,6,18)
        for j = 1:npSets
            sumAXP = exp_ATP{j} + exp_ADP{j} + exp_AMP{j};
            if j == setup.refParams
                darr = dprofile{1}(1):0.005:dprofile{1}(8);
                parr = pchip(dprofile{j},(exp_ATP{j} + exp_ADP{j}/2)./sumAXP,darr);
                % plot the value
                plot(darr, parr, 'k-','LineWidth',1)
                hold on
                % plot the points
                plot(dprofile{j}(1:q(1)), (exp_ATP{j} + exp_ADP{j}/2)./sumAXP, 'k.','MarkerSize',5)

%                 plot(dprofile{j},(exp_ATP{j} + exp_ADP{j}/2)./sumAXP, 'k','LineWidth',1.2)
            else
                plot(dprofile{j},(exp_ATP{j} + exp_ADP{j}/2)./sumAXP, 'color', colorSet(j,:))
            end
        hold on
%         plot(canelas_SS.mtlD.D, (canelas_SS.mtlD.ATP + canelas_SS.mtlD.ADP/2) ./ (canelas_SS.mtlD.ATP + canelas_SS.mtlD.ADP + canelas_SS.mtlD.AMP), 'r+')
        plot(canelas_SS.mtlD.D, (canelas_SS.mtlD.ATP + canelas_SS.mtlD.ADP/2) ./ (canelas_SS.mtlD.ATP + canelas_SS.mtlD.ADP + canelas_SS.mtlD.AMP),...
            'kd',...
            'MarkerFaceColor','k','MarkerSize',3)
        end
%         title('Energy charge: (ATP + ADP/2) / sumAXP')
%         set(gca,'XTickLabel',[]);
        titleName = 'Energy charge';
        text(sp.XLim(2)*0.1,0.93,titleName,'FontSize',8)
        yticks(sp.YLim)
            % arrows issues
            tempText = text(0.03,0.78,'\uparrow');
            tempText.FontSize = 20;
            tempText = text(0.03+0.05,0.78,'1');
            tempText.FontSize = 10;
        
        % set(gcf,'color','w'); suptitle('fig21. Physiology. Steady state profile. All')
% % % %         set(gcf,'color',backColor);
    end
    end

    
    %fig3. Dynamic glucose pulse metabolite profile. All
    for simGPconc = 1
    if setup.runGSvanHeerden == 1
%         figure('name','Dynamic glucose pulse metabolite profile. All','units','normalized','outerposition',[0 0 1 1])
%         figure(5004)
        fh_5004 = figure(5004);
        fh_5004.Position = [1990 400 1100 524];
        
%         for k = 1:length(legenda.metabolites)
%             subplot(4,8,k)
        for k_pre = 1:length(setup.metabolitesGP_2plot)
            k = setup.metabolitesGP_2plot(k_pre);
            k_loc = setup.metabolitesGP_location(k_pre);
% % % %             subplot(8,4,k)
% % % %             subplot(9,4,k_loc)
            sp = subplot(3,6,k_loc);
            
            % plot simulations
            for j = 1:npSets
                if j == setup.refParams
                    plot(Tgs{j}, Ygs{j}(:,k), 'k-','LineWidth',1)
                else
                    plot(Tgs{j}, Ygs{j}(:,k), 'color', colorSet(j,:))
                end
                hold on
            end
            % plot experimental data
            if data.ordered.metabolites.onoff(k) == 1
%                 plot(data.ordered.metabolites.time{k},data.ordered.metabolites.data{k},'r+')
                plot(data.ordered.metabolites.time{k},data.ordered.metabolites.data{k},...
                    'ks',...
                    'MarkerFaceColor','k','MarkerSize',3)
            end
%             titleName = erase(legenda.metabolites{k},", [mM]");
%             title(titleName)%fontsize{16}
            
            if k_loc == 14
                ylim([0 30])
            elseif k_loc == 4
                ylim([0 0.1])
            elseif k_loc == 10
                ylim([0 3])
            end
%             set(gca,'XTickLabel',[]);
%             set(gca,'YTickLabel',[]);
            titleName = erase(legenda.metabolites{k},", [mM]");
%             titleName = titleName(1:end-7);
            if k_loc == 9
                text(-80,68,titleName,'FontSize',8)
            else
                text(-80,sp.YLim(2)*0.9,titleName,'FontSize',8)
            end
            yticks(sp.YLim)
            xlim([-100 340])
            if k_loc == 5
                ylim([0 4])
            end
            if((k_loc == 10)||(k_loc == 14))
                % arrows issues
                tempText = text(0,-sp.YLim(2)/10,'\uparrow');
                tempText.FontSize = 20;
                tempText = text(30,-sp.YLim(2)/7,'3');
                tempText.FontSize = 10;
            end            
            
        end
        
        sp = subplot(3,6,16); % NAD:NADH ratio
        plot(Tgs{1}, Ygs{1}(:,7)./Ygs{1}(:,8), 'k-','LineWidth',1)
        titleName = 'NAD:NADH ratio';
        text(-80,sp.YLim(2)*0.9,titleName,'FontSize',8)
        yticks(sp.YLim)
        xlim([-100 340])
        ylim([0 4000])
                
        sp = subplot(3,6,17); % sumAXP
        plot(Tgs{1}, Ygs{1}(:,9) + Ygs{1}(:,15) + Ygs{1}(:,16), 'k-','LineWidth',1)
        ylim([0 4])
        titleName = 'sum AXP';
        text(-80,sp.YLim(2)*0.9,titleName,'FontSize',8)
        yticks(sp.YLim)
        xlim([-100 340])
        
        sp = subplot(3,6,18); % sumIXP
        plot(Tgs{1}, Ygs{1}(:,28) + Ygs{1}(:,29) + Ygs{1}(:,30), 'k-','LineWidth',1)
        ylim([0 3])
        titleName = 'sum IXP';
        text(-80,sp.YLim(2)*0.9,titleName,'FontSize',8)
        yticks(sp.YLim)
        xlim([-100 340])

%         suptitle('fig3. Dynamic glucose pulse metabolite profile. All')
        % set(gcf,'color','w');
% % % %         set(gcf,'color',backColor);
    end
    end
    
    
    %fig4. Dynamic glucose pulse flux profile. All
    for simGPrates = 1
    if setup.runGSvanHeerden == 1
%         figure('name','Dynamic glucose pulse flux profile. All','units','normalized','outerposition',[0 0 1 1])
%         for k = 1:length(legenda.fluxes)
%             subplot(5,9,k)
%         figure(5005)
        fh_5005 = figure(5005);
        fh_5005.Position = [1995 676 771 303];
        
        for k_pre = 1:length(setup.ratesGP_2plot)
            k = setup.ratesGP_2plot(k_pre);
            k_loc = setup.ratesGP_location(k_pre);
            sp = subplot(2,5,k_loc);
            
            % plot simulations
            for j = 1:npSets
                if j == setup.refParams
                    plot(Tgs{j}, Vgs{j}(:,k), 'k','LineWidth',1)
                else
                    plot(Tgs{j}, Vgs{j}(:,k), 'color', colorSet(j,:))
                end
                hold on
            end
            % plot experimental data
            if data.ordered.fluxes.onoff(k) == 1
%                 plot(data.ordered.fluxes.time{k},data.ordered.fluxes.data{k},'r+')
                plot(data.ordered.fluxes.time{k},data.ordered.fluxes.data{k},...
                    'k^',...
                    'MarkerFaceColor','k','MarkerSize',3)
            end
%             titleName = erase(legenda.fluxes{k},", [mM s^{-1}]");
%             title(titleName)%fontsize{16}

%             set(gca,'XTickLabel',[]);
%             set(gca,'YTickLabel',[]);
            titleName = erase(legenda.fluxes{k},", [mM s^{-1}]");
%             if k_loc == 
%             titleName = titleName(1:end-7);
            if k_loc <= 5
                titleName = titleName(1:end-7);
            else
                titleName = titleName(1:end-8);
            end
%                 text(-80,68,titleName,'FontSize',8)
%             else
                text(-80,sp.YLim(2)*0.9,titleName,'FontSize',8)
%             end
            yticks(sp.YLim)
            xlim([-100 340])
            if k_loc == 8
                ylim([0 0.05])
            elseif k_loc == 9
                ylim([0 0.025])
            elseif((k_loc == 4)||(k_loc == 5))
                ylim([0 0.31])
            end
            
        end 
        
%         suptitle('fig4. Dynamic glucose pulse flux profile. All')
        % set(gcf,'color','w');
% % % %         set(gcf,'color',backColor);
    end
    end
    
    
end

% cleaning up things
numFigs = setup.plotResultsMode; % add here the resulting validation I would say
end


% 
% % saving the plots
% savefig(5001,'results\manuscript_supplementary_materials\supF_SSmets.fig')
% savefig(5002,'results\manuscript_supplementary_materials\supF_SSrates.fig')
% savefig(5003,'results\manuscript_supplementary_materials\supF_SSfysio.fig')
% savefig(5004,'results\manuscript_supplementary_materials\supF_GPmets.fig')
% savefig(5005,'results\manuscript_supplementary_materials\supF_GPrates.fig')



% % % % %% MF5 HXK MPSA
% % % % clc, clear, close all
% % % % set_paths
% % % % path(path, strcat(folder,'/tempResults'))
% % % % dbstop if error
% % % % colorpalette;
% % % % 
% % % % %% Recall starting setup
% % % % for starting_setup = 1
% % % %     %
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
% % % %     % 
% % % % end
% % % % % latest setups options
% % % % setup.GLT_km_same = 1;
% % % % setup.sameGLTHXK4all = 1;
% % % % setup.GLT_km_same_extra = 1;
% % % % load('x108a.mat','x_1st','x_2nd'), x_2nd(35) = x_2nd(38); %x_2nd_prev = x_2nd;
% % % % % load('x108b.mat','x_1st','x_2nd')
% % % % % load('x108c.mat','x_1st','x_2nd'),



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
%% Starting simulation
% simulation
setup.plotResultsMode = 0;
[simRes_ref] = simulateY3M1_separateGPSS(x_2nd, canelas_SS, data, dataset, setup);
setup.simGPdata = simRes_ref{1}.gs;
setup.simSSdata = simRes_ref{1}.ss;
setup.litParams = zeros(size(x));
setup.plotResultsMode = 0;
% % plotting
% simRes = [simRes_ref,];
% xRes = [x_2nd;];
% setup.plotResultsMode = 30;
% [~] = plotAll_Y3M1(simRes,legenda,canelas_SS,data,setup,xRes);
% setup.plotResultsMode = 0;


% % %% Adjusted simulation: Manual adjustment parameters
% % setup.VmGLT_mu_dependent = 1;
% % % 
% % setup.plotResultsMode = 0;
% % [simAdj] = simulateY3M1_separateGPSS(x_2nd, canelas_SS, data, dataset, setup);
% % setup.plotResultsMode = 0;
% % % 
% % close all
% % simRes = [simRes_ref, simAdj];
% % xRes = [x_2nd; x_2nd];
% % setup.plotResultsMode = 30;
% % [~] = plotAll_Y3M1(simRes,legenda,canelas_SS,data,setup,xRes);
% % setup.plotResultsMode = 0;
% % 
% % 
% % %% plotting Exp vs Sim HXT
% % fh_2 = figure(2);
% % D = fh_2.Children(43).Children(1).XData;
% % HXK_exp = fh_2.Children(43).Children(1).YData;
% % % HXT_exp = fh_2.Children(43).Children(1).YData;
% % HXK_sim_single_pset = fh_2.Children(43).Children(3).YData;
% % HXK_sim_mu_dependent = fh_2.Children(43).Children(2).YData;
% % 
% % %% plot exp vs sim
% % % close(5000)
% % fh_5000 = figure(5000); 
% % fh_5000.Position = [2108 525 480 420]; 
% % % 
% % maxVal = 1.1 * 1.5; 
% % line([0 maxVal],[0 maxVal],'Color','r','LineStyle','--')
% % hold on
% % % 
% % plot(HXK_sim_single_pset,HXK_exp,...
% %     'o','MarkerSize',5,...
% %     'MarkerEdgeColor','k','MarkerFaceColor',[.7 .7 .7])
% % % 
% % plot(HXK_sim_mu_dependent,HXK_exp,...
% %     'o','MarkerSize',5,...
% %     'MarkerEdgeColor','k','MarkerFaceColor','k')
% % % 
% % xlim([0 maxVal])
% % ylim([0 maxVal])
% % % title('PDC')
% % % text(sp10.XLim(2) * 0.075, sp10.YLim(2) * 0.925,'PDC')
% % box on
% % tick_points = [0 .2 .5 1 1.5];
% % xticks(tick_points), xlim([0 maxVal])
% % yticks(tick_points), ylim([0 maxVal])
% % % 
% % xlabel('simulated reaction rate (mmol L^{-1} s^{-1})')
% % ylabel('experimental reaction rate (mmol L^{-1} s^{-1})')
% % % %% annotations
% % for i = 1:length(HXK_exp)
% %     temp_x = [HXK_sim_single_pset(i) HXK_sim_mu_dependent(i)];
% %     temp_y = [HXK_exp(i) HXK_exp(i)];
% %     temp_lin = line(temp_x,temp_y,'LineStyle','--','Color','k');
% % end
% % % %%
% % % embeding focused figure southeast
% % newh_fh_5000 = copyobj(fh_5000.Children,5000);
% % % 
% % newh_fh_5000.Position = [0.55 0.17 0.32 0.32];
% % newh_fh_5000.FontSize = 9;
% % newh_fh_5000.XLabel = [];
% % newh_fh_5000.YLabel = [];
% % newh_fh_5000.XTick = .2;
% % newh_fh_5000.YTick = .2;
% % % xticks([0.2])
% % % yticks([0.2])
% % % set(newh_fh_5000,'xtick',[])
% % % set(newh_fh_5000,'ytick',[])
% % newh_fh_5000.XLim = [0 0.22];
% % newh_fh_5000.YLim = [0 0.22];
% % % xtickpos = get(gca, 'xtick');
% % 
% % 
% % %% plot data vs dilRate
% % close(5001)
% % dilRate = canelas_SS.mtlD.D;
% % Cs = canelas_SS.mtlD.Cs;
% % fh_5001 = figure(5001); 
% % % fh_5001.Position = [2616 527 480 420];
% % fh_5001.Position = [2616 527 588 420]; 
% % % 
% % % maxVal = 1.1 * 1.5; 
% % % line([0 maxVal],[0 maxVal],'Color','r','LineStyle','--')
% % plot(dilRate,HXK_exp,'r.--')
% % hold on
% % % 
% % plot(dilRate,HXK_sim_single_pset,...
% %     'o','MarkerSize',5,...
% %     'MarkerEdgeColor','k','MarkerFaceColor',[.7 .7 .7])
% % % 
% % plot(dilRate,HXK_sim_mu_dependent,...%HXK_exp,...
% %     'o','MarkerSize',5,...
% %     'MarkerEdgeColor','k','MarkerFaceColor','k')
% % % %%  
% % % title('PDC')
% % % text(sp10.XLim(2) * 0.075, sp10.YLim(2) * 0.925,'PDC')
% % box on
% % maxVal_x = 1.5 * max(dilRate); 
% % maxVal_y = 1.2 * max(HXK_exp); 
% % xticks([0 0.1 0.2 0.3 0.4]), xlim([0 maxVal_x])
% % yticks([0 .2 .5 1.0 1.5]), ylim([0 maxVal_y])
% % % %% 
% % xlabel(Y3M1_labels.dilRate)
% % ylabel(Y3M1_labels.qT_mMs)
% % % 
% % txt1 = text(0.06,1.515,'A');
% % txt1.FontSize = 15;
% % % % %%
% % % % % %% annotations
% % % % for i = 1:length(HXK_exp)
% % % %     temp_x = [HXK_sim_single_pset(i) HXK_sim_mu_dependent(i)];
% % % %     temp_y = [HXK_exp(i) HXK_exp(i)];
% % % %     temp_lin = line(temp_x,temp_y,'LineStyle','--','Color','k');
% % % % end
% % % % %%
% % % % % embeding focused figure southeast
% % % % newh_fh_5001 = copyobj(fh_5001.Children,5001);
% % % % % 
% % % % % newh_fh_5001.Position = [0.6 0.17 0.27 0.20];
% % % % newh_fh_5001.Position = [.175 .32 0.24 0.24];
% % % % newh_fh_5001.FontSize = 9;
% % % % newh_fh_5001.XLabel = [];
% % % % newh_fh_5001.YLabel = [];
% % % % newh_fh_5001.XTick = .1;
% % % % newh_fh_5001.YTick = .2;
% % % % % xticks([0.2])
% % % % % yticks([0.2])
% % % % % set(newh_fh_5001,'xtick',[])
% % % % % set(newh_fh_5001,'ytick',[])
% % % % newh_fh_5001.XLim = [0 0.125];
% % % % newh_fh_5001.YLim = [0 0.3];
% % % % % xtickpos = get(gca, 'xtick');
% % 
% % % %% parameter bar plot
% % bplot_vals = [-1.2 -0.5 -0.3 0,...
% %     1 1 1 1];
% % bplot_c = categorical({'0.02','0.05','0.10','0.20',...
% %                 '0.30','0.325','0.35','0.375'});
% % 
% % % ax3 = axes('Position',[.25 .675 .24 .2]);
% % ax3 = axes('Position',[.2 .5 .24 .2]);
% % box on
% % cp3 = bar(bplot_c,bplot_vals,'FaceColor','k');
% % xtickangle(75)
% % % 
% % txt2 = text(1,1.25,'B');
% % txt2.FontSize = 15;
% % 
% % ax4 = axes('Position',[.625 .175 .24 .25]);
% % box on
% % plot(dilRate,Cs,'r.--')
% % % ax4.XTick = [0 2];
% % ax4.YTick = [0 2];
% % txt4 = text(ax4.XLim(2)*0.1, ax4.YLim(2)*0.8, 'C');
% % txt4.FontSize = 15;
% % 
% % % dilRate = canelas_SS.mtlD.D;
% % % Cs = canelas_SS.mtlD.Cs;
% % 
% % % %% ONLY IF SAVING
% % % savefig_loc = 'E:\OneDrive - TU Eindhoven\Documents\ch4_Y3M1\results\figures\manuscript\';
% % % saveName = [savefig_loc, 'fig5_hxt_d_dependent_linear']; savefig(fh_5000, saveName);
% % % 
% % % % savefig(1, 'tempfig21.fig')
% % % % specs printing (method 3)
% % % set(gcf,'Units','inches');
% % % screenposition = get(gcf,'Position');
% % % set(gcf,...
% % %     'PaperPosition',[0 0 screenposition(3:4)],...
% % %     'PaperSize',[screenposition(3:4)]);
% % % print -dpdf -painters fig5_hxt_d_dependent_linear
% % % print -dpng -painters fig5_hxt_d_dependent_linear
% % % %% ONLY IF SAVING
% % % savefig_loc = 'E:\OneDrive - TU Eindhoven\Documents\ch4_Y3M1\results\figures\manuscript\';
% % % saveName = [savefig_loc, 'fig5_hxt_d_dependent_dilRate']; savefig(fh_5001, saveName);
% % % 
% % % % savefig(1, 'tempfig21.fig')
% % % % specs printing (method 3)
% % % set(gcf,'Units','inches');
% % % screenposition = get(gcf,'Position');
% % % set(gcf,...
% % %     'PaperPosition',[0 0 screenposition(3:4)],...
% % %     'PaperSize',[screenposition(3:4)]);
% % % print -dpdf -painters fig5_hxt_d_dependent_dilRate
% % % print -dpng -painters fig5_hxt_d_dependent_dilRate
% % 
% % %% Meshgrid sampling of parameters (Vm and Km)
% % disp('here')
% % %% Plotting meshgrid



%%
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% % % %PROPERLY DEVELOPMENT OF THE 2D PLOTS % % % % % % % % % % % % % % % %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 

GLT_exp = canelas_SS.mtlD.v_GLT;
% create range Vm and Km
range_x_Vm = [-1:0.2:3]';
range_x_Km = [-2:0.2:2]';
% range_x_Vm2 = [-1:0.05:3]';
% range_x_Km2 = [-2:0.05:2]';
[mesh_Vm, mesh_Km] = meshgrid(range_x_Vm,range_x_Km);
% real values
real_range_x_Vm = 10 .^ range_x_Vm .* 3.67 * 0.2;
real_range_x_Km = 10 .^ range_x_Km .* 1.1918;
real_mesh_Vm = 10 .^ mesh_Vm .* 3.67 * 0.2;
real_mesh_Km = 10 .^ mesh_Km .* 1.1918;
% 
[j_array, i_array] = size(mesh_Vm);

% % % empty simulation and errors matrixes
% % vglt_sim_Vm_Km_mesh = cell(size(mesh_Vm));
% % error_full = cell(size(mesh_Vm));
% % error_full_abs = cell(size(mesh_Vm));
% % error_full_abs_sum = zeros(size(mesh_Vm));
% % error_d1 = zeros(size(mesh_Vm));
% % error_d2 = zeros(size(mesh_Vm));
% % error_d3 = zeros(size(mesh_Vm));
% % error_d4 = zeros(size(mesh_Vm));
% % error_d5 = zeros(size(mesh_Vm));
% % error_d6 = zeros(size(mesh_Vm));
% % error_d7 = zeros(size(mesh_Vm));
% % error_d8 = zeros(size(mesh_Vm));
% % % 
% % [j_array, i_array] = size(mesh_Vm);
% % 
% % % simulate
% % [j_array, i_array] = size(mesh_Vm);
% % cluster = parcluster('local');
% % pool = parpool(cluster,4);
% % parfor i = 1:i_array
% %     for j = 1:j_array
% %         tempName_pre = sprintf('simRes_mF5B_HXT_j%d_i%d.mat',j,i);
% %         if exist(tempName_pre,'file') == 2
% %         else
% %             %
% %             sprintf('simulation, row=%d, column=%d.',j,i)
% %             % 
% %             xtemp = x_2nd;
% %             xtemp(36) = mesh_Vm(j,i);
% %             xtemp(38) = mesh_Km(j,i);
% %             xtemp(35) = mesh_Km(j,i);
% %             [simRes_temp] = simulateY3M1_separateGPSS(xtemp, canelas_SS, data, dataset, setup);
% %             % 
% %             temp = zeros(1,8);
% %             for k = 1:8
% %                 temp(k) = simRes_temp{1}.ss.V.ss2can{k}(end,1);
% %             end
% %     %         vglt_sim_Vm_Km_mesh{j,i} = temp;
% %             save_MF5_array(j,i,temp)
% %         end
% %     end
% % end
% % delete(pool)


% % %% 
% % % 
% % GLT_exp = canelas_SS.mtlD.v_GLT;
% % % create range Vm and Km
% % % range_x_Vm = [-1:0.2:3]';
% % % range_x_Km = [-2:0.2:2]';
% % range_x_Vm = [-1:0.05:3]';
% % range_x_Km = [-2:0.05:2]';
% % [mesh_Vm, mesh_Km] = meshgrid(range_x_Vm,range_x_Km);
% % % real values
% % real_range_x_Vm = 10 .^ range_x_Vm .* 3.67 * 0.2;
% % real_range_x_Km = 10 .^ range_x_Km .* 1.1918;
% % real_mesh_Vm = 10 .^ mesh_Vm .* 3.67 * 0.2;
% % real_mesh_Km = 10 .^ mesh_Km .* 1.1918;
% % 
% % % empty simulation and errors matrixes
% % vglt_sim_Vm_Km_mesh = cell(size(mesh_Vm));
% % error_full = cell(size(mesh_Vm));
% % error_full_abs = cell(size(mesh_Vm));
% % error_full_abs_sum = zeros(size(mesh_Vm));
% % error_d1 = zeros(size(mesh_Vm));
% % error_d2 = zeros(size(mesh_Vm));
% % error_d3 = zeros(size(mesh_Vm));
% % error_d4 = zeros(size(mesh_Vm));
% % error_d5 = zeros(size(mesh_Vm));
% % error_d6 = zeros(size(mesh_Vm));
% % error_d7 = zeros(size(mesh_Vm));
% % error_d8 = zeros(size(mesh_Vm));
% % % 
% % [j_array, i_array] = size(mesh_Vm);
% % 
% % % simulate
% % [j_array, i_array] = size(mesh_Vm);
% % cluster = parcluster('local');
% % pool = parpool(cluster,16);
% % parfor i = 1:i_array
% %     for j = 1:j_array
% % %         tempName_pre = sprintf('simRes_mF5B_HXT_j%d_i%d.mat',j,i);
% %         tempName_pre = sprintf('simRes_mF5B2_HXT_j%d_i%d.mat',j,i);
% %         if exist(tempName_pre,'file') == 2
% %         else
% %             %
% %             sprintf('simulation, row=%d, column=%d.',j,i)
% %             % 
% %             xtemp = x_2nd;
% %             xtemp(36) = mesh_Vm(j,i);
% %             xtemp(38) = mesh_Km(j,i);
% %             xtemp(35) = mesh_Km(j,i);
% %             [simRes_temp] = simulateY3M1_separateGPSS(xtemp, canelas_SS, data, dataset, setup);
% %             % 
% %             temp = zeros(1,8);
% %             for k = 1:8
% %                 temp(k) = simRes_temp{1}.ss.V.ss2can{k}(end,1);
% %             end
% %     %         vglt_sim_Vm_Km_mesh{j,i} = temp;
% %             save_MF5_array(j,i,temp)
% %         end
% %     end
% % end
% % delete(pool)


% %% recall error data (run 1, about 21 dtps)
% for i = 1:i_array
%     for j = 1:j_array
%         % 
% %         loadName = sprintf('simRes_mF5B2_HXT_j%d_i%d.mat',j,i);
%         loadName = sprintf('simRes_mF5B_HXT_j%d_i%d.mat',j,i);
%         load(loadName);
%         vglt_sim_Vm_Km_mesh{j,i} = temp;
%         % 
%         tempError = zeros(1,8);
%         for k = 1:8
% %             tempError(k)= (vglt_sim_Vm_Km_mesh{j,i}(k) - vglt_exp(k))./vglt_exp(k);
%             tempError(k)= (vglt_sim_Vm_Km_mesh{j,i}(k) - GLT_exp(k))./GLT_exp(k);
%         end
%         error_full{j,i} = tempError;
%         error_full_abs{j,i} = abs(error_full{j,i});
%         error_full_abs_sum(j,i) = sum(error_full_abs{j,i});
%         % 
%         error_d1(j,i) = error_full_abs{j,i}(1);
%         error_d2(j,i) = error_full_abs{j,i}(2);
%         error_d3(j,i) = error_full_abs{j,i}(3);
%         error_d4(j,i) = error_full_abs{j,i}(4);
%         error_d5(j,i) = error_full_abs{j,i}(5);
%         error_d6(j,i) = error_full_abs{j,i}(6);
%         error_d7(j,i) = error_full_abs{j,i}(7);
%         error_d8(j,i) = error_full_abs{j,i}(8);
%     end
% end
% error_d_layers = zeros([size(mesh_Vm), 8]);
% error_d_layers(:,:,1) = error_d1;
% error_d_layers(:,:,2) = error_d2;
% error_d_layers(:,:,3) = error_d3;
% error_d_layers(:,:,4) = error_d4;
% error_d_layers(:,:,5) = error_d5;
% error_d_layers(:,:,6) = error_d6;
% error_d_layers(:,:,7) = error_d7;
% error_d_layers(:,:,8) = error_d8;


%% %
% % % % save('simRes_mF5B_HXT_meshSim.mat',...
% % % %     'vglt_sim_Vm_Km_mesh',...
% % % %     'error_full','error_full_abs','error_full_abs_sum',...
% % % %     'error_d1','error_d2','error_d3','error_d4','error_d5','error_d6','error_d7','error_d8'),
% save('simRes_mF5_HXT_meshSim.mat',...
%     'vglt_sim_Vm_Km_mesh',...
%     'error_full','error_full_abs','error_full_abs_sum',...
%     'error_d1','error_d2','error_d3','error_d4','error_d5','error_d6','error_d7','error_d8'),
% %% 
load('simRes_mF5_HXT_meshSim.mat',...
    'vglt_sim_Vm_Km_mesh',...
    'error_full','error_full_abs','error_full_abs_sum',...
    'error_d1','error_d2','error_d3','error_d4','error_d5','error_d6','error_d7','error_d8'),
error_d_layers = zeros([size(mesh_Vm), 8]);
error_d_layers(:,:,1) = error_d1;
error_d_layers(:,:,2) = error_d2;
error_d_layers(:,:,3) = error_d3;
error_d_layers(:,:,4) = error_d4;
error_d_layers(:,:,5) = error_d5;
error_d_layers(:,:,6) = error_d6;
error_d_layers(:,:,7) = error_d7;
error_d_layers(:,:,8) = error_d8;


%% plotting
d = canelas_SS.mtlD.D;
% preparation
yvalues = cellstr(num2str(mesh_Km(:,end))); 
xvalues = cellstr(num2str(mesh_Vm(1,:)'));
% finding closest value
% % x_2nd(38) = -0.0728 -> for Km
% % x_2nd(36) = 1.0445 -> for Vmax
perc_Vm = (x_2nd(36) - min(range_x_Vm)) / (max(range_x_Vm) - min(range_x_Vm));
perc_Km = (x_2nd(38) - min(range_x_Km)) / (max(range_x_Km) - min(range_x_Km));
% adjustments
Vm_adjusted = zeros(1,8);
    Vm_adjusted(1) = -0.05;
    Vm_adjusted(2) = 0.8;
    Vm_adjusted(3) = x_2nd(36);
    Vm_adjusted(4) = x_2nd(36);
    Vm_adjusted(5) = 1.15;
    Vm_adjusted(6) = 1.5;
    Vm_adjusted(7) = 1.26;
    Vm_adjusted(8) = 1.2;
% % 
% x_2nd(36) - Vm_adjusted' =
%     1.0945
%     0.2445
%          0
%          0
%    -0.1055
%    -0.4555
%    -0.2155
%    -0.1555
%   
% for the tick labels
xticks_nums = [-1 0 1 2 3];
yticks_nums = flip([-2 -1 0 1 2]);
len_temp = length(xticks_nums);
temp_xticks = cell(1,len_temp);
temp_yticks = cell(1,len_temp);
for i = 1:len_temp
    temp_xticks{i} = [strrep(num2str(10 .^ xticks_nums(i)' .* 3.67 * 0.2,'%10.2e\n'),'e','^{'),'}'];
    temp_yticks{i} = [strrep(num2str(10 .^ yticks_nums(i)' .* 1.1918,'%10.2e\n'),'e','^{'),'}'];
end


% plotting
if exist('figh_1041','var')
    clf(figh_1041)
end
figh_1041 = figure(1041);
figh_1041.Position = [1921 257 1536 600];%[1995 471 1317 400]; 
% figure(1042)
for i = 1:8
    sp_h = subplot(2,4,i);
% for i = 1
%     sp_h = subplot(1,1,i);
    
    % heatmap
    hm_h = heatmap(xvalues, yvalues, error_d_layers(:,:,i),...
        'GridVisible','off',...
        'ColorbarVisible','off');%,...
%         'Colormap', gray);
%     str = ['d = ',num2str(d(i)),' h^{-1}'];
    str = ['d = ',num2str(d(i))];
    
%     str = 'y_1 = x^2 and y_2 = 2x^{2 + k}';
    title(str)
%     subtitle(str)
    Ax = gca;
    Ax.XDisplayLabels = nan(size(Ax.XDisplayData));
    Ax.YDisplayLabels = nan(size(Ax.YDisplayData));
    
    % plot adjusted + estimated value
    ax = axes;
    ax.Position = hm_h.Position;
    % 
    lin_h = line([min([x_2nd(36) Vm_adjusted(i)]) max([x_2nd(36) Vm_adjusted(i)])],...
        [x_2nd(38) x_2nd(38)],...
        'Color','black','LineStyle','--','LineWidth',1);
    hold on
    % 
    pt_h_adjusted = plot(Vm_adjusted(i), x_2nd(38), 'ko',...
        'MarkerFaceColor','k','linewidth',1);
    % 
    pt_h = plot(x_2nd(36), x_2nd(38), 'ko',...
        'MarkerFaceColor','w','linewidth',1);
    xlim([min(range_x_Vm) max(range_x_Vm)])
    ylim([min(range_x_Km) max(range_x_Km)])
    ax.Color = 'none';
    
    % labels
    xlabel('Vm')
    ylabel('Km')
    if(i<5)
        xlabel([])
    end
    if((i==2)||(i==3)||(i==4)||(i==6)||(i==7)||(i==8))
        ylabel([])
    end
    
    % ticks and ticklabels, + moving location big
    % 
    xticks(xticks_nums)
    xticklabels(temp_xticks)
%     xtickangle(30)
    if(i<5)
        xticklabels([])
    end
    
    % 
    yticks_nums = [-2 -1 0 1 2];
    yticks(yticks_nums)
    yticklabels(temp_yticks)
%     ytickangle(30)
    if((i==2)||(i==3)||(i==4)||(i==6)||(i==7)||(i==8))
        yticklabels([])
    end
    
end
% lines literature value (correction 0.2 factor in Vmax)
% option of a 3D plot




%% Adjusted simulation: Manual adjustment parameters
% % adjusted values in map
% Vm_adjusted = zeros(1,8);
%     Vm_adjusted(1) = -0.05;
%     Vm_adjusted(2) = 0.8;
%     Vm_adjusted(3) = x_2nd(36); % 1.0445
%     Vm_adjusted(4) = x_2nd(36); % 1.0445
%     Vm_adjusted(5) = 1.15;
%     Vm_adjusted(6) = 1.5;
%     Vm_adjusted(7) = 1.26;
%     Vm_adjusted(8) = 1.2;
% % % adjusted values in
% % x_2nd(36) - Vm_adjusted' =
% %     1.0945
% %     0.2445
% %          0
% %          0
% %    -0.1055
% %    -0.4555
% %    -0.2155
% %    -0.1555
% %   

setup.VmGLT_mu_dependent = 1;
% 
setup.plotResultsMode = 0;
[simAdj] = simulateY3M1_separateGPSS(x_2nd, canelas_SS, data, dataset, setup);
setup.plotResultsMode = 0;
%% 
% close all
simRes = [simRes_ref, simAdj];
xRes = [x_2nd; x_2nd];
setup.plotResultsMode = 30;
[~] = plotAll_Y3M1(simRes,legenda,canelas_SS,data,setup,xRes);
setup.plotResultsMode = 0;


%% plotting Exp vs Sim HXT
fh_2 = figure(2);
D = fh_2.Children(43).Children(1).XData;
HXK_exp = fh_2.Children(43).Children(1).YData;
% HXT_exp = fh_2.Children(43).Children(1).YData;
HXK_sim_single_pset = fh_2.Children(43).Children(3).YData;
HXK_sim_mu_dependent = fh_2.Children(43).Children(2).YData;

%% plot exp vs sim
% close(5000)
fh_5000 = figure(5000); 
fh_5000.Position = [2108 525 480 420]; 
% 
maxVal = 1.1 * 1.5; 
line([0 maxVal],[0 maxVal],'Color','r','LineStyle','--')
hold on
% 
plot(HXK_sim_single_pset,HXK_exp,...
    'o','MarkerSize',5,...
    'MarkerEdgeColor','k','MarkerFaceColor',[.7 .7 .7])
% 
plot(HXK_sim_mu_dependent,HXK_exp,...
    'o','MarkerSize',5,...
    'MarkerEdgeColor','k','MarkerFaceColor','k')
% 
xlim([0 maxVal])
ylim([0 maxVal])
% title('PDC')
% text(sp10.XLim(2) * 0.075, sp10.YLim(2) * 0.925,'PDC')
box on
tick_points = [0 .2 .5 1 1.5];
xticks(tick_points), xlim([0 maxVal])
yticks(tick_points), ylim([0 maxVal])
% 
xlabel('simulated reaction rate (mmol L^{-1} s^{-1})')
ylabel('experimental reaction rate (mmol L^{-1} s^{-1})')
% %% annotations
for i = 1:length(HXK_exp)
    temp_x = [HXK_sim_single_pset(i) HXK_sim_mu_dependent(i)];
    temp_y = [HXK_exp(i) HXK_exp(i)];
    temp_lin = line(temp_x,temp_y,'LineStyle','--','Color','k');
end
% %%
% embeding focused figure southeast
newh_fh_5000 = copyobj(fh_5000.Children,5000);
% 
newh_fh_5000.Position = [0.55 0.17 0.32 0.32];
newh_fh_5000.FontSize = 9;
newh_fh_5000.XLabel = [];
newh_fh_5000.YLabel = [];
newh_fh_5000.XTick = .2;
newh_fh_5000.YTick = .2;
% xticks([0.2])
% yticks([0.2])
% set(newh_fh_5000,'xtick',[])
% set(newh_fh_5000,'ytick',[])
newh_fh_5000.XLim = [0 0.22];
newh_fh_5000.YLim = [0 0.22];
% xtickpos = get(gca, 'xtick');


%% plot data vs dilRate
% close(5001)
dilRate = canelas_SS.mtlD.D;
Cs = canelas_SS.mtlD.Cs;
fh_5001 = figure(5001); 
% fh_5001.Position = [2616 527 480 420];
fh_5001.Position = [2616 527 588 420]; 
% 
% maxVal = 1.1 * 1.5; 
% line([0 maxVal],[0 maxVal],'Color','r','LineStyle','--')
plot(dilRate,HXK_exp,'r.--')
hold on
% 
plot(dilRate,HXK_sim_single_pset,...
    'o','MarkerSize',5,...
    'MarkerEdgeColor','k','MarkerFaceColor',[.7 .7 .7])
% 
plot(dilRate,HXK_sim_mu_dependent,...%HXK_exp,...
    'o','MarkerSize',5,...
    'MarkerEdgeColor','k','MarkerFaceColor','k')
% %%  
% title('PDC')
% text(sp10.XLim(2) * 0.075, sp10.YLim(2) * 0.925,'PDC')
box on
maxVal_x = 1.5 * max(dilRate); 
maxVal_y = 1.2 * max(HXK_exp); 
xticks([0 0.1 0.2 0.3 0.4]), xlim([0 maxVal_x])
yticks([0 .2 .5 1.0 1.5]), ylim([0 maxVal_y])
% %% 
xlabel(Y3M1_labels.dilRate)
ylabel(Y3M1_labels.qT_mMs)
% 
txt1 = text(0.0625,1.515,'A');
txt1.FontSize = 15;
% % %%
% % % %% annotations
% % for i = 1:length(HXK_exp)
% %     temp_x = [HXK_sim_single_pset(i) HXK_sim_mu_dependent(i)];
% %     temp_y = [HXK_exp(i) HXK_exp(i)];
% %     temp_lin = line(temp_x,temp_y,'LineStyle','--','Color','k');
% % end
% % %%
% % % embeding focused figure southeast
% % newh_fh_5001 = copyobj(fh_5001.Children,5001);
% % % 
% % % newh_fh_5001.Position = [0.6 0.17 0.27 0.20];
% % newh_fh_5001.Position = [.175 .32 0.24 0.24];
% % newh_fh_5001.FontSize = 9;
% % newh_fh_5001.XLabel = [];
% % newh_fh_5001.YLabel = [];
% % newh_fh_5001.XTick = .1;
% % newh_fh_5001.YTick = .2;
% % % xticks([0.2])
% % % yticks([0.2])
% % % set(newh_fh_5001,'xtick',[])
% % % set(newh_fh_5001,'ytick',[])
% % newh_fh_5001.XLim = [0 0.125];
% % newh_fh_5001.YLim = [0 0.3];
% % % xtickpos = get(gca, 'xtick');

% %% parameter bar plot
safecopy_x36 = 1.0445;
bplot_vals = [-0.05-safecopy_x36 0.8-safecopy_x36 0 0,...
    1.15-safecopy_x36 1.5-safecopy_x36 1.26-safecopy_x36 1.2-safecopy_x36];
bplot_c = categorical({'0.02','0.05','0.10','0.20',...
                '0.30','0.325','0.35','0.375'});

% ax3 = axes('Position',[.25 .675 .24 .2]);
% ax3 = axes('Position',[.2 .5 .24 .275]);
ax3 = axes('Position',[.64 .225 .24 .275]);
box on
cp3 = bar(bplot_c,bplot_vals,'FaceColor','k');
xtickangle(75)
ax3.YLim = [-1.25 1.25];
% ax3.YTick = [-1 0 1];
% 
txt2 = text(1,0.9,'C');
txt2.FontSize = 15;

% ax4 = axes('Position',[.625 .175 .24 .25]);
ax4 = axes('Position',[.19 .46 .24 .29]);
box on
plot(dilRate,Cs,'r.--')
% ax4.XTick = [0 2];
ax4.YTick = [0 2];
txt4 = text(ax4.XLim(2)*0.1, ax4.YLim(2)*0.8, 'B');
txt4.FontSize = 15;

% dilRate = canelas_SS.mtlD.D;
% Cs = canelas_SS.mtlD.Cs;

% % % % %% ONLY IF SAVING: fig5000
% % % % savefig_loc = 'E:\OneDrive - TU Eindhoven\Documents\ch4_Y3M1\results\figures\manuscript\';
% % % % saveName = [savefig_loc, 'fig5_hxt_d_dependent_linear']; savefig(fh_5000, saveName);
% % % % figure(5000)
% % % % 
% % % % % savefig(1, 'tempfig21.fig')
% % % % % specs printing (method 3)
% % % % set(gcf,'Units','inches');
% % % % screenposition = get(gcf,'Position');
% % % % set(gcf,...
% % % %     'PaperPosition',[0 0 screenposition(3:4)],...
% % % %     'PaperSize',[screenposition(3:4)]);
% % % % print -dpdf -painters fig5_hxt_d_dependent_linear
% % % % print -dpng -painters fig5_hxt_d_dependent_linear
% % % % 
% % % % %% ONLY IF SAVING: fig5001
% % % % savefig_loc = 'E:\OneDrive - TU Eindhoven\Documents\ch4_Y3M1\results\figures\manuscript\';
% % % % saveName = [savefig_loc, 'fig5_hxt_d_dependent_dilRate']; savefig(fh_5001, saveName);
% % % % figure(5001)
% % % % 
% % % % % savefig(1, 'tempfig21.fig')
% % % % % specs printing (method 3)
% % % % set(gcf,'Units','inches');
% % % % screenposition = get(gcf,'Position');
% % % % set(gcf,...
% % % %     'PaperPosition',[0 0 screenposition(3:4)],...
% % % %     'PaperSize',[screenposition(3:4)]);
% % % % print -dpdf -painters fig5_hxt_d_dependent_dilRate
% % % % print -dpng -painters fig5_hxt_d_dependent_dilRate
% % % % 
% % % % %% ONLY IF SAVING: fig1041
% % % % savefig_loc = 'E:\OneDrive - TU Eindhoven\Documents\ch4_Y3M1\results\figures\manuscript\';
% % % % saveName = [savefig_loc, 'fig5sup_hxt_d_dependent_MPSA']; savefig(figh_1041, saveName);
% % % % figure(1041)
% % % % 
% % % % % savefig(1, 'tempfig21.fig')
% % % % % specs printing (method 3)
% % % % set(gcf,'Units','inches');
% % % % screenposition = get(gcf,'Position');
% % % % set(gcf,...
% % % %     'PaperPosition',[0 0 screenposition(3:4)],...
% % % %     'PaperSize',[screenposition(3:4)]);
% % % % print -dpdf -painters fig5sup_hxt_d_dependent_MPSA
% % % % print -dpng -painters fig5sup_hxt_d_dependent_MPSA




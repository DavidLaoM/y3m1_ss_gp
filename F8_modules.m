% % F7_HXT_DEPPEND.m
% This figure reproduces the growth rate dependency case for enzyme HXT.
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
setup.plotResultsMode = 30;
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
% [refSim] = simulateY3M1_separateGPSS(x_2nd, canelas_SS, data, dataset, setup);
[refSim] = simulateY3M1_separateGPSS(x, canelas_SS, data, dataset, setup);
setup.plotResultsMode = 0;

% %% plotting
setup.plotResultsMode = 30;
% namesHits = cell(1,2); namesHits{1} = 'start'; namesHits{2} = 'x108res';
namesHits = cell(1,1); namesHits{1} = 'x108res'; %namesHits{2} = 'x108res';
setup.namesHits = namesHits;
% simRes = [simRes_reference, refSim];
simRes = refSim;
% %%
% [~] = plotAll_Y3M1(simRes,legenda,canelas_SS,data,setup,x_2nd);
[~] = plotAll_Y3M1(simRes,legenda,canelas_SS,data,setup,x);
setup.plotResultsMode = 0;


%% glycolytic flux
% % (table/barplot quantification) Glycolytic flux: Start calculating ratios 
% % for the main flux and the sinks (1) Info missed if no sinks (account for 
% % TRE,PPP in SS), (2) info missed on respiration and fermentation.
% (1) Percentage of flux from uptake going to sinks vs fermentation(?). 
% Bar plots at the different growth rates
% (2) Percentage of flux from uptake to respiration and fermentitation. 
% Bar plots at the different growth rates

nS = 8;
% C-mol incoming
temp = zeros(1,nS); for i = 1:nS, temp(i) = simRes{1}.ss.V.ss2can{i}(end,1); end
Cmol_vGLT = 6 * temp;
% C-mol each sink + CO2_PDC
temp = zeros(1,nS); for i = 1:nS, temp(i) = simRes{1}.ss.V.ss2can{i}(end,34); end
Cmol_vsG6P = 6 * temp;
temp = zeros(1,nS); for i = 1:nS, temp(i) = simRes{1}.ss.V.ss2can{i}(end,35); end
Cmol_vsF6P = 6 * temp;
temp = zeros(1,nS); for i = 1:nS, temp(i) = simRes{1}.ss.V.ss2can{i}(end,36); end
Cmol_vsGAP = 3 * temp;
temp = zeros(1,nS); for i = 1:nS, temp(i) = simRes{1}.ss.V.ss2can{i}(end,37); end
Cmol_vsP3G = 3 * temp;
temp = zeros(1,nS); for i = 1:nS, temp(i) = simRes{1}.ss.V.ss2can{i}(end,38); end
Cmol_vsPEP = 3 * temp;
temp = zeros(1,nS); for i = 1:nS, temp(i) = simRes{1}.ss.V.ss2can{i}(end,39); end
Cmol_vsPYR = 3 * temp;
temp = zeros(1,nS); for i = 1:nS, temp(i) = simRes{1}.ss.V.ss2can{i}(end,40); end
Cmol_vsACE = 2 * temp;
temp = zeros(1,nS); for i = 1:nS, temp(i) = simRes{1}.ss.V.ss2can{i}(end,13); end
Cmol_CO2_PDC = 1 * temp;
% C-mol total sink
Cmol_vs_total = Cmol_vsG6P - Cmol_vsF6P - Cmol_vsGAP + Cmol_vsP3G + Cmol_vsPEP + Cmol_vsPYR + Cmol_vsACE;
% C-mol ethanol branch, glycerol branch
temp = zeros(1,nS); for i = 1:nS, temp(i) = simRes{1}.ss.V.ss2can{i}(end,23); end
Cmol_etohT = 2 * temp;
temp = zeros(1,nS); for i = 1:nS, temp(i) = simRes{1}.ss.V.ss2can{i}(end,24); end
Cmol_glycT = 3 * temp;
% C-mol fermentation (ethanol + glycerol?) + respiration (sink pyruvate)
Cmol_ferm = Cmol_etohT + Cmol_glycT + Cmol_CO2_PDC;
Cmol_resp = Cmol_vsPYR;

% %% check: values in the order of 10-7
% Cmol_vGLT - Cmol_vs_total - Cmol_CO2_PDC - Cmol_etohT - Cmol_glycT
% % % % %%
% percent_rest = 100 - percent_sink;
percent_ferm = (Cmol_etohT + Cmol_glycT)./Cmol_vGLT * 100;
percent_CO2 = Cmol_CO2_PDC./Cmol_vGLT * 100;
% 
Cmol_vs_total_pos = Cmol_vsG6P + Cmol_vsP3G + Cmol_vsPEP + Cmol_vsPYR + Cmol_vsACE;
percent_sinks_nopyr_overSinksPos = (Cmol_vs_total_pos - Cmol_vsPYR)./Cmol_vs_total_pos * 100;
percent_sinks_pyr_overSinksPos = Cmol_vsPYR./Cmol_vs_total_pos * 100;
percent_sinks_g6p_overSinksPos = Cmol_vsG6P./Cmol_vs_total_pos * 100;
percent_sinks_p3g_overSinksPos = Cmol_vsP3G./Cmol_vs_total_pos * 100;
percent_sinks_pep_overSinksPos = Cmol_vsPEP./Cmol_vs_total_pos * 100;
percent_sinks_ace_overSinksPos = Cmol_vsACE./Cmol_vs_total_pos * 100;
% 
percent_sinks_all = Cmol_vs_total./Cmol_vGLT * 100;
% percent_sinks_nopyr = (Cmol_vs_total - Cmol_vsPYR)./Cmol_vGLT * 100;
% percent_sinks_pyr = Cmol_vsPYR./Cmol_vGLT * 100;
percent_sinks_nopyr = percent_sinks_all .* percent_sinks_pyr_overSinksPos /100; % option not considering the negative direction F6P, GAP.
percent_sinks_pyr = percent_sinks_all .* (1 - percent_sinks_pyr_overSinksPos /100); % option not considering the negative direction F6P, GAP.


% % % % %%
c = categorical({'0.02','0.05','0.10','0.20',...
                '0.30','0.325','0.35','0.375'});
% 
% clf(7001)
if exist('fh_7001','var')
    clf(fh_7001)
end
fh_7001 = figure(7001);
fh_7001.Position = [1986 415 1375 504];
% % 
% subplot(1,3,1)
% bp1 = bar(c,[percent_sink; percent_rest]','stacked','FaceColor','flat');
% legend('sink reactions','rest','location','SouthOutside','Orientation','Horizontal')
% ylabel('percentage of carbon flux')
% xlabel('dilution rate, h^{-1}')
% xtickangle(30)
% yticks([0 25 50 75 100])
% 
subplot(1,3,3)
% cp3 = bar(c,[percent_sinks_pyr; percent_sinks_nopyr; percent_CO2; percent_ferm]','stacked','FaceColor','flat');
% legend('sinks(pyr)','sinks (No pyr)','ferm. (CO_{2}, v_{PDC})','ferm. (etoh_T+glyc_T)','location','SouthOutside','Orientation','Horizontal')
cp3 = bar(c,[percent_sinks_pyr + percent_sinks_nopyr; percent_CO2 + percent_ferm]','stacked','FaceColor','flat');
legend('sinks','ferm','location','SouthOutside','Orientation','Horizontal')
ylabel('percentage of carbon flux')
xlabel('dilution rate, h^{-1}')
xtickangle(30)
ylim([0 100])
yticks([0 25 50 75 100])
for i = 1:nS
% % % %     cp3(1).CData(i,:) = [44,127,184]/255; % [1 1 1];
% % % %     cp3(2).CData(i,:) = [127,205,187]/255;
% % % %     cp3(3).CData(i,:) = [253,187,132]/255;
% % % %     cp3(4).CData(i,:) = [227,74,51]/255;
    cp3(1).CData(i,:) = [44,127,184]/255;
    cp3(2).CData(i,:) = [227,74,51]/255;
end

% 
subplot(1,3,1)
bp1 = bar(c,[percent_sinks_pyr_overSinksPos; percent_sinks_nopyr_overSinksPos]','stacked','FaceColor','flat');
legend('sink (pyr)','sinks (no pyr)','location','SouthOutside','Orientation','Horizontal')
ylabel('percentage of carbon flux')
xlabel('dilution rate, h^{-1}')
xtickangle(30)
yticks([0 25 50 75 100])
ylim([0 100])
for i = 1:nS
    bp1(1).CData(i,:) = [44,127,184]/255; % [1 1 1];
    bp1(2).CData(i,:) = [127,205,187]/255;
end

% 
subplot(1,3,2)
bp2 = bar(c,[percent_sinks_pyr_overSinksPos; percent_sinks_g6p_overSinksPos; percent_sinks_p3g_overSinksPos; percent_sinks_pep_overSinksPos; percent_sinks_ace_overSinksPos]','stacked','FaceColor','flat');
legend('sink (pyr)','sink (g6p)','sink (p3g)','sink (pep)','sink (ace)'...
    ,'location','SouthOutside','Orientation','Horizontal')
ylabel('percentage of carbon flux')
xlabel('dilution rate, h^{-1}')
xtickangle(30)
yticks([0 25 50 75 100])
ylim([0 100])
for i = 1:nS
    bp2(1).CData(i,:) = [44,127,184]/255; % [1 1 1];
    bp2(2).CData(i,:) = [127,205,187]/255;
end

% % 
% set(fh_7001,'color','w')

% % % % % 
% % % % if setup.save == 1
% % % %     saveName = [savefig_loc, 'fig8a_sinks.fig'];
% % % %     savefig(fh_7001, saveName)
% % % % %     % 
% % % % %     set(gcf,'Units','inches');
% % % % %     screenposition = get(3,'Position');
% % % % %     set(gcf,...
% % % % %         'PaperPosition',[0 0 screenposition(3:4)],...
% % % % %         'PaperSize',[screenposition(3:4)]);
% % % % %     print -dpdf -painters Y3M1_robustness_rates
% % % % %     print -dpng -painters Y3M1_robustness_rates    
% % % % end



%% respiration branch
% NOTE, TO CALCULATE THE ATP THINK, ALSO DISCOUNT THE ATP GENERATED BY THE
% RESPIRATION BY MULTIPLYING BY THE QCO2 (coming from vPDC CO2 released).


% % (table/barplot quantification/small diagram around) (Could be already 
% % covered by the previous) so cutting presence of oxygen -> PDC and sinkPYR 
% % reaction. I wonder if this and another implementations will be valid for 
% % the experts.
% (1) Simulate cutting oxygen (sinkPYR)
% (2) See both GP and SS simulations
% (3) Decide on: (table/barplot quantification/small diagram around)
close(1:6)

% simulation
setup.oxygen_OFF = 1;
% [sim_oxOFF] = simulateY3M1_separateGPSS(x_2nd, canelas_SS, data, dataset, setup);
[sim_oxOFF] = simulateY3M1_separateGPSS(x, canelas_SS, data, dataset, setup);
setup.oxygen_OFF = 0;

% plotting
setup.plotResultsMode = 30;
namesHits = cell(1,2); namesHits{1} = 'reference'; namesHits{2} = 'oxOFF';
simRes = [refSim, sim_oxOFF];
% [~] = plotAll_Y3M1(simRes,legenda,canelas_SS,data,setup,[x_2nd; x_2nd]);
[~] = plotAll_Y3M1(simRes,legenda,canelas_SS,data,setup,[x; x]);
setup.plotResultsMode = 0;

%% visualization
c = categorical({'0.02','0.05','0.10','0.20',...
                '0.30','0.325','0.35','0.375'});
% 
% clf(7002)
if exist('fh_7002','var')
    clf(fh_7002)
end
fh_7002 = figure(7002);
fh_7002.Position = [1986 415 1375 504];
% 
sp3 = subplot(1,3,3);
cp3 = bar(c,[percent_sinks_pyr; percent_sinks_nopyr; percent_CO2; percent_ferm]','stacked','FaceColor','flat');
% legend('resp(sinks(pyr))','sinks (No pyr)','ferm. (CO_{2}, v_{PDC})','ferm. (etoh_T+glyc_T)','location','SouthOutside','Orientation','Horizontal')
legend('resp(sinks(pyr))','sinks (No pyr)','(CO_{2}, v_{PDC})','ferm. (etoh_T+glyc_T)','location','SouthOutside','Orientation','Horizontal')
ylabel('percentage of carbon flux')
xlabel('dilution rate, h^{-1}')
xtickangle(30)
ylim([0 100])
yticks([0 25 50 75 100])
for i = 1:nS
    cp3(1).CData(i,:) = [44,127,184]/255; % [1 1 1];
    cp3(2).CData(i,:) = [254,232,200]/255;
    cp3(3).CData(i,:) = [253,187,132]/255;
    cp3(4).CData(i,:) = [227,74,51]/255;
end

% node_pyr
node_pyr_resp = zeros(1,nS); for i = 1:nS, node_pyr_resp(i) = simRes{1}.ss.V.ss2can{i}(end,39); end
node_pyr_ferm = zeros(1,nS); for i = 1:nS, node_pyr_ferm(i) = simRes{1}.ss.V.ss2can{i}(end,13); end
node_pyr_tot = node_pyr_resp + node_pyr_ferm;

% 
sp2 = subplot(1,3,2);
bp2 = bar(c,[100*node_pyr_resp./node_pyr_tot; 100*node_pyr_ferm./node_pyr_tot]','stacked','FaceColor','flat');
% legend('resp,sink_{pyr}','ferm,v_{PDC}','location','SouthOutside','Orientation','Horizontal')
legend('resp,sink_{pyr}','v_{PDC}','location','SouthOutside','Orientation','Horizontal')
ylabel('percentage of carbon flux')
xlabel('dilution rate, h^{-1}')
xtickangle(30)
yticks([0 25 50 75 100])
ylim([0 100])
for i = 1:nS
    bp2(1).CData(i,:) = [44,127,184]/255; % [1 1 1];
    bp2(2).CData(i,:) = [253,187,132]/255;
end

% what if resp (sinkPYR) off?
fh_2 = figure(2);
vPDCexp = fh_2.Children(31).Children(1).YData;
vPDCsim_respON = fh_2.Children(31).Children(3).YData;
vPDCsim_respOFF = fh_2.Children(31).Children(2).YData;

% 
figure(7002)
sp1 = subplot(1,3,1);
scatter(vPDCexp,vPDCsim_respON,40,'MarkerEdgeColor','k',...
              'MarkerFaceColor',[0 0 0],...
              'LineWidth',1.5)
hold on
scatter(vPDCexp,vPDCsim_respOFF,40,'MarkerEdgeColor','k',...
              'MarkerFaceColor',[1 1 1],...
              'LineWidth',1.5)
line([0 2.5],[0 2.5],'linestyle','--','color','k')
xlabel('v_{PDC,exp}')
ylabel('v_{PDC,sim}')
legend('resp_{ON}','resp_{OFF}','location','northwest')
box on

% %%
figure(7002), 
sp1.Position = [0.1300    0.3226    0.2134    0.6024];
fh_7002.Children(5).Position = [0.5011 0.0502 0.4036 0.0526];
fh_7002.Children(6).Position = [0.6916 0.3587 0.2134 0.5663]; % rerun

% % 
% set(fh_7002,'color','w')

% % % % % 
% % % % if setup.save == 1
% % % %     saveName = [savefig_loc, 'fig8b_resp.fig'];
% % % %     savefig(fh_7002, saveName)
% % % % %     % 
% % % % %     set(gcf,'Units','inches');
% % % % %     screenposition = get(3,'Position');
% % % % %     set(gcf,...
% % % % %         'PaperPosition',[0 0 screenposition(3:4)],...
% % % % %         'PaperSize',[screenposition(3:4)]);
% % % % %     print -dpdf -painters Y3M1_robustness_rates
% % % % %     print -dpng -painters Y3M1_robustness_rates    
% % % % end


%% trehalose cycle
% % (simulation/sampling?) Trehalose cycle knock-out.
close(1:6)

% (1) Classic assay
% x_treOFF = x_2nd; x_treOFF([86 126 120 123]) = -10 * ones;
x_treOFF = x; x_treOFF([86 126 120 123]) = -10 * ones;

% simulation
[sim_treOFF] = simulateY3M1_separateGPSS(x_treOFF, canelas_SS, data, dataset, setup);

% plotting
setup.plotResultsMode = 30;
namesHits = cell(1,2); namesHits{1} = 'reference'; namesHits{2} = 'treOFF';
simRes = [refSim, sim_treOFF];
% [~] = plotAll_Y3M1(simRes,legenda,canelas_SS,data,setup,[x_2nd; x_treOFF]);
[~] = plotAll_Y3M1(simRes,legenda,canelas_SS,data,setup,[x; x_treOFF]);
setup.plotResultsMode = 0;

%% visualization: recall data
fh_3 = figure(3);
fh_4 = figure(4);

% pt1 balanced metabolites
f6p_sim_conc = fh_3.Children(30).Children(3).YData;
f6p_sim_time = fh_3.Children(30).Children(3).XData;
atp_sim_conc = fh_3.Children(25).Children(3).YData;
atp_sim_time = fh_3.Children(25).Children(3).XData;
t6p_sim_conc = fh_3.Children(8).Children(3).YData;
t6p_sim_time = fh_3.Children(8).Children(3).XData;
% 
f6p_exp_conc = fh_3.Children(30).Children(1).YData;
f6p_exp_time = fh_3.Children(30).Children(1).XData;
atp_exp_conc = fh_3.Children(25).Children(1).YData;
atp_exp_time = fh_3.Children(25).Children(1).XData;
t6p_exp_conc = fh_3.Children(8).Children(1).YData;
t6p_exp_time = fh_3.Children(8).Children(1).XData;

% pt2 imbalanced metabolites
f6p_sim_conc_imb = fh_3.Children(30).Children(2).YData;
f6p_sim_time_imb = fh_3.Children(30).Children(2).XData;
atp_sim_conc_imb = fh_3.Children(25).Children(2).YData;
atp_sim_time_imb = fh_3.Children(25).Children(2).XData;
t6p_sim_conc_imb = fh_3.Children(8).Children(2).YData;
t6p_sim_time_imb = fh_3.Children(8).Children(2).XData;

% pt3 balanced rates
vupp_hxk_sim_time = fh_4.Children(43).Children(3).XData;
vupp_hxk_sim_rate = fh_4.Children(43).Children(3).YData;
vupp_hxk_exp_time = fh_4.Children(43).Children(1).XData;
vupp_hxk_exp_rate = fh_4.Children(43).Children(1).YData;
vupp_gapdh_sim_time = fh_4.Children(37).Children(2).XData;
vupp_gapdh_sim_rate = fh_4.Children(37).Children(2).YData;
vupp_g3pdh_sim_time = fh_4.Children(39).Children(2).XData;
vupp_g3pdh_sim_rate = fh_4.Children(39).Children(2).YData;

% pt4 imbalanced rates
vupp_hxk_sim_time_imb = fh_4.Children(43).Children(2).XData;
vupp_hxk_sim_rate_imb = fh_4.Children(43).Children(2).YData;
vupp_gapdh_sim_time_imb = fh_4.Children(37).Children(1).XData;
vupp_gapdh_sim_rate_imb = fh_4.Children(37).Children(1).YData;
vupp_g3pdh_sim_time_imb = fh_4.Children(39).Children(1).XData;
vupp_g3pdh_sim_rate_imb = fh_4.Children(39).Children(1).YData;


%% visualization: plotting
% pt1, pt2, remember to normalize
% clf(7003)
if exist('fh_7003','var')
    clf(fh_7003)
end
fh_7003 = figure(7003);
fh_7003.Position = [1986 280 793 707];

% pt1 balanced metabolites
sp1 = subplot(2,2,1);
% atp
% % % % plot(atp_sim_time_imb,atp_sim_conc_imb,'b--','LineWidth',1.5)
hold on
plot(atp_sim_time,atp_sim_conc,'b-','LineWidth',2)
% % % % plot(atp_exp_time,atp_exp_conc,'bo','MarkerFaceColor','b','MarkerSize',4)
% f6p
% % % % plot(f6p_sim_time_imb,f6p_sim_conc_imb,'r--','LineWidth',1.5)
hold on
plot(f6p_sim_time,f6p_sim_conc,'r-','LineWidth',2)
% % % % plot(f6p_exp_time,f6p_exp_conc,'ro','MarkerFaceColor','r','MarkerSize',4)
% 
xlim([-100 340])
ylim([0 3.5])
yticks([0 1 2 3])
xlabel(Y3M1_labels.time)
ylabel(Y3M1_labels.conc_mMs)
txt2 = text(-90,3.3,'ATP','FontSize',13,'Color','b');
txt3 = text(-90,3.0,'F6P','FontSize',13,'Color','r');
box on
title('balanced concentrations')

% addition
% handaxes1 = axes('Position', [0.695 0.8125 a1.Position(3)/3.2 a1.Position(4)/3.2]);
handaxes1 = axes('Position', [0.35 0.81 0.1 0.1]);
% t6p
% % % % plot(t6p_sim_time_imb,t6p_sim_conc_imb,'k--','LineWidth',2)
hold on
plot(t6p_sim_time,t6p_sim_conc,'k-','LineWidth',2)
% % % % plot(t6p_exp_time,t6p_exp_conc,'ko','MarkerFaceColor','k','MarkerSize',4)
handaxes1.XLim = [-100 340];
handaxes1.YLim = [0 10];
set(gca,'XTickLabel',[],'YTickLabel',[]);
txt1 = text(-90,8.5,'TRE','FontSize',13);
box on
    

% pt2 balanced metabolites
sp2 = subplot(2,2,2);
% atp
plot(atp_sim_time_imb,atp_sim_conc_imb,'b','LineWidth',1.5)
hold on
% % % % plot(atp_sim_time,atp_sim_conc,'b-','LineWidth',2)
% % % % plot(atp_exp_time,atp_exp_conc,'bo','MarkerFaceColor','b','MarkerSize',4)
% f6p
plot(f6p_sim_time_imb,f6p_sim_conc_imb,'r','LineWidth',1.5)
hold on
% % % % plot(f6p_sim_time,f6p_sim_conc,'r-','LineWidth',2)
% % % % plot(f6p_exp_time,f6p_exp_conc,'ro','MarkerFaceColor','r','MarkerSize',4)
% 
xlim([-100 340])
ylim([0 3.5])
yticks([0 1 2 3])
xlabel(Y3M1_labels.time)
ylabel(Y3M1_labels.conc_mMs)
txt2 = text(-90,3.3,'ATP','FontSize',13,'Color','b');
txt3 = text(-90,3.0,'F6P','FontSize',13,'Color','r');
box on
title('imbalanced concentrations')

% addition
% handaxes1 = axes('Position', [0.695 0.8125 a1.Position(3)/3.2 a1.Position(4)/3.2]);
handaxes1 = axes('Position', [0.79 0.81 0.1 0.1]);
% t6p
plot(t6p_sim_time_imb,t6p_sim_conc_imb,'k','LineWidth',2)
hold on
% % % % plot(t6p_sim_time,t6p_sim_conc,'k-','LineWidth',2)
% % % % plot(t6p_exp_time,t6p_exp_conc,'ko','MarkerFaceColor','k','MarkerSize',4)
handaxes1.XLim = [-100 340];
handaxes1.YLim = [0 10];
set(gca,'XTickLabel',[],'YTickLabel',[]);
txt1 = text(-90,8.5,'TRE','FontSize',13);
box on

% pt1 balanced metabolites
sp3 = subplot(2,2,3);
% % atp
% hold on
% plot(vupp_hxk_sim_time,vupp_hxk_sim_rate,'b--','LineWidth',2)
% % f6p
% hold on
% plot(vupp_gapdh_sim_time,vupp_gapdh_sim_rate,'r--','LineWidth',2)
% vHXK - vGAPDH
% plot(vupp_hxk_sim_time,vupp_hxk_sim_rate - vupp_gapdh_sim_rate - vupp_g3pdh_sim_rate,'k--','LineWidth',2.5)
% % % % plot(vupp_hxk_sim_time,vupp_hxk_sim_rate - vupp_gapdh_sim_rate,'k--','LineWidth',2.5)
plot(vupp_hxk_sim_time,vupp_hxk_sim_rate - vupp_gapdh_sim_rate,'k--','LineWidth',2.5)
hold on
plot(vupp_hxk_sim_time_imb, vupp_hxk_sim_rate_imb - vupp_gapdh_sim_rate_imb,'--','LineWidth',2.5,'color',[.6 .6 .6])
line([-100 340],[0 0],'LineWidth',0.5,'Linestyle','-','color','k')
% 
xlim([-100 340])
ylim([-1 2])
yticks([-1 0 1 2])
xlabel(Y3M1_labels.time)
ylabel(Y3M1_labels.q_mumolgDWs)
% txt2 = text(-90,0.5,'HXK','FontSize',13,'Color','b');
% txt3 = text(-90,0.25,'GAPDH','FontSize',13,'Color','r');
% % % % txt3 = text(-90,0.25,'GAPDH','FontSize',13,'Color','r');
box on
% title({'balanced rates:';'\color{red}addition? ATP consumed vs produced'}) 
% % % % title('balanced rates:') 
% % % % text(-90,1.75,'addition? ATP consumed vs produced','FontSize',10,'Color','r');

% pt1 balanced metabolites
sp4 = subplot(2,2,4);
% vHXK - vGAPDH
% plot(vupp_hxk_sim_time_imb, vupp_hxk_sim_rate_imb - vupp_gapdh_sim_rate_imb - vupp_g3pdh_sim_rate_imb,'k--','LineWidth',2.5)
plot(vupp_hxk_sim_time_imb, vupp_hxk_sim_rate_imb - vupp_gapdh_sim_rate_imb,'k--','LineWidth',2.5)
hold on
line([-100 340],[0 0],'LineWidth',0.5,'Linestyle','-','color','k')
% 
xlim([-100 340])
ylim([-1 2])
yticks([-1 0 1 2])
xlabel(Y3M1_labels.time)
ylabel(Y3M1_labels.q_mumolgDWs)
% txt2 = text(-90,0.5,'HXK','FontSize',13,'Color','b');
% txt3 = text(-90,0.25,'GAPDH','FontSize',13,'Color','r');
% % % % txt3 = text(-90,0.25,'GAPDH','FontSize',13,'Color','r');
box on
% title({'imbalanced rates;';' \color{red}negative due to sinks (fixed to WT)';'\color{blue}but close overlap plot.'})     
title('imbalanced rates;')     
text(-90,1.75,'negative due to sinks (fixed to WT)','FontSize',10,'Color','r');
text(-90,1.5,'but close overlap plot.','FontSize',10,'Color','b');

%
box on

% % 
% set(fh_7003,'color','w')

% % % % % 
% % % % if setup.save == 1
% % % %     saveName = [savefig_loc, 'fig8c_tre.fig'];
% % % %     savefig(fh_7003, saveName)
% % % % %     % 
% % % % %     set(gcf,'Units','inches');
% % % % %     screenposition = get(3,'Position');
% % % % %     set(gcf,...
% % % % %         'PaperPosition',[0 0 screenposition(3:4)],...
% % % % %         'PaperSize',[screenposition(3:4)]);
% % % % %     print -dpdf -painters Y3M1_robustness_rates
% % % % %     print -dpng -painters Y3M1_robustness_rates    
% % % % end



%% inosine salvage pathway
close(1:6)
% % Inosine salvage pathway knock-out? At least trials would still run.
% (1) Knock out to show only the mis fit in AXP balance then
% x_inoOFF = x_2nd; x_inoOFF([112 115]) = -10 * ones;
% x_inoOFF = x_2nd; x_inoOFF(112:118) = -10 * ones;
x_inoOFF = x; x_inoOFF(112:118) = -10 * ones;

% simulation
[sim_inoOFF] = simulateY3M1_separateGPSS(x_inoOFF, canelas_SS, data, dataset, setup);

% plotting
setup.plotResultsMode = 30;
namesHits = cell(1,2); namesHits{1} = 'reference'; namesHits{2} = 'inoOFF';
simRes = [refSim, sim_inoOFF];
% [~] = plotAll_Y3M1(simRes,legenda,canelas_SS,data,setup,[x_2nd; x_inoOFF]);
[~] = plotAll_Y3M1(simRes,legenda,canelas_SS,data,setup,[x; x_inoOFF]);
setup.plotResultsMode = 0;

%%
fh_3 = figure(3);
% pt1 balanced metabolites
atp_sim_conc = fh_3.Children(25).Children(3).YData;
atp_sim_time = fh_3.Children(25).Children(3).XData;
%  exp data
atp_exp_conc = fh_3.Children(25).Children(1).YData;
atp_exp_time = fh_3.Children(25).Children(1).XData;
% pt2 imbalanced metabolites
atp_sim_conc_imb = fh_3.Children(25).Children(2).YData;
atp_sim_time_imb = fh_3.Children(25).Children(2).XData;


%% visualization: plotting
% pt1, pt2, remember to normalize
% clf(7003)
if exist('fh_7003b','var')
    clf(fh_7003b)
end
fh_7003b = figure(70030);
fh_7003b.Position = [2056 580 322 264];

% pt1 balanced metabolites
% sp1 = subplot(2,2,1);
% 
plot(atp_sim_time_imb,atp_sim_conc_imb,'-','LineWidth',2,'color',[.7 .7 .7])
hold on
plot(atp_sim_time,atp_sim_conc,'k-','LineWidth',2)
plot(atp_exp_time,atp_exp_conc,'ko','MarkerFaceColor','k','MarkerSize',4)
% %% 
xlim([-100 340])
ylim([0 2.5])
yticks([0 1 2 3])
xlabel(Y3M1_labels.time)
ylabel(Y3M1_labels.conc_mMs)
box on
% title('balanced concentrations')
%%
% % % % % 
% % % % if setup.save == 1
% % % %     saveName = [savefig_loc, 'fig8c2_ixp.fig'];
% % % %     savefig(fh_7003b, saveName)
% % % % %     % 
% % % % %     set(gcf,'Units','inches');
% % % % %     screenposition = get(3,'Position');
% % % % %     set(gcf,...
% % % % %         'PaperPosition',[0 0 screenposition(3:4)],...
% % % % %         'PaperSize',[screenposition(3:4)]);
% % % % %     print -dpdf -painters Y3M1_robustness_rates
% % % % %     print -dpng -painters Y3M1_robustness_rates    
% % % % end



%% growth rate dependency
close(1:6)
% % (some barplot threshold) Growth rate dependency in cofactors. For 
% % instance. NADH regeneration capacity.
% (1) Simulate the lack of mitochondrial increase (x132:136 = 0)
% (2) Then decide on a plot (simulation or threshold?) or 
% also (table/barplot quantification/small diagram around)
% x_mitoOFF = x_2nd; x_mitoOFF(132:136) = -10 * ones;
x_mitoOFF = x; x_mitoOFF(132:136) = -10 * ones;

% simulation
[sim_mitoOFF] = simulateY3M1_separateGPSS(x_mitoOFF, canelas_SS, data, dataset, setup);

% plotting
setup.plotResultsMode = 30;
namesHits = cell(1,2); namesHits{1} = 'reference'; namesHits{2} = 'mitoOFF';
simRes = [refSim, sim_mitoOFF];
% [~] = plotAll_Y3M1(simRes,legenda,canelas_SS,data,setup,[x_2nd; x_mitoOFF]);
[~] = plotAll_Y3M1(simRes,legenda,canelas_SS,data,setup,[x; x_mitoOFF]);
setup.plotResultsMode = 0;

% visualization

%% only NADH
% also (table/barplot quantification/small diagram around)
% x_NADmitoOFF = x_2nd; x_NADmitoOFF(132) = -10 * ones;
x_NADmitoOFF = x; x_NADmitoOFF(132) = -10 * ones;

% simulation
[sim_NADmitoOFF] = simulateY3M1_separateGPSS(x_NADmitoOFF, canelas_SS, data, dataset, setup);

% plotting
setup.plotResultsMode = 30;
namesHits = cell(1,2); namesHits{1} = 'reference'; namesHits{2} = 'NADmitoOFF';
simRes = [refSim, sim_NADmitoOFF];
% [~] = plotAll_Y3M1(simRes,legenda,canelas_SS,data,setup,[x_2nd; x_NADmitoOFF]);
[~] = plotAll_Y3M1(simRes,legenda,canelas_SS,data,setup,[x; x_NADmitoOFF]);
setup.plotResultsMode = 0;

close(3)
close(4)
close(6)

%%
% 
fh_2 = figure(2);
% 
vGLK_sim = fh_2.Children(42).Children(3).YData;
vGLK_sim_NADmitoOFF = fh_2.Children(42).Children(2).YData;
vGLK_exp = fh_2.Children(42).Children(1).YData;
% 
vGPD_sim = fh_2.Children(38).Children(3).YData;
vGPD_sim_NADmitoOFF = fh_2.Children(38).Children(2).YData;
vGPD_exp = fh_2.Children(38).Children(1).YData;
% v(7)= + v_G3PDH - v_GAPDH + v_ADH + v_mitoNADH;%NAD
% v(8)= - v_G3PDH + v_GAPDH - v_ADH - v_mitoNADH;%NADH
vGAPDH_sim = fh_2.Children(36).Children(3).YData;
vADH_sim = fh_2.Children(3).Children(3).YData;
vmitoNADH_sim = fh_2.Children(18).Children(2).YData;
% - vGAPDH_sim + vADH_sim + vmitoNADH_sim + vGPD_sim

%% visualization
% clf(7004)
if exist('fh_7004','var')
    clf(fh_7004)
end
fh_7004 = figure(7004);
fh_7004.Position = [1986 436 1375 483];
% 
sp1 = subplot(1,3,1);
scatter(vGLK_exp,vGLK_sim,40,'MarkerEdgeColor','k',...
              'MarkerFaceColor',[0 0 0],...
              'LineWidth',1.5)
hold on
scatter(vGLK_exp,vGLK_sim_NADmitoOFF,40,'MarkerEdgeColor','k',...
              'MarkerFaceColor',[1 1 1],...
              'LineWidth',1.5)
line([0 2],[0 2],'linestyle','--','color','k')
xlabel('v_{GLK,exp}')
ylabel('v_{GLK,sim}')
legend('v_{GLK,sim}','v_{GLK,sim,mitoOFF}','location','northwest')
box on
% 
sp2 = subplot(1,3,2);
scatter(vGPD_exp,vGPD_sim,40,'MarkerEdgeColor','k',...
              'MarkerFaceColor',[0 0 0],...
              'LineWidth',1.5)
hold on
scatter(vGPD_exp,vGPD_sim_NADmitoOFF,40,'MarkerEdgeColor','k',...
              'MarkerFaceColor',[1 1 1],...
              'LineWidth',1.5)
line([0 0.2],[0 0.2],'linestyle','--','color','k')
xlabel('v_{GPD,exp}')
ylabel('v_{GPD,sim}')
legend('v_{GPD,sim}','v_{GPD,sim,mitoOFF}','location','southeast')
% xlim([0 0.1])
% ylim([0 0.1])
box on
% 
c_sp3 = categorical({'0.02','0.05','0.10','0.20',...
                '0.30','0.325','0.35','0.375'});
sp3 = subplot(1,3,3); % - vGAPDH_sim + vADH_sim + vmitoNADH_sim + vGPD_sim
cp3 = bar(c_sp3,100*[vGPD_sim./vGAPDH_sim; vADH_sim./vGAPDH_sim; vmitoNADH_sim./vGAPDH_sim]','stacked','FaceColor','flat');
legend('v_{GPD}','v_{ADH}','v_{mitoNADH}','location','SouthOutside','Orientation','Horizontal')
ylabel('percentage of NADH produced')
xlabel('dilution rate, h^{-1}')
title('NADH recycle')
xtickangle(30)
ylim([0 100])
yticks([0 25 50 75 100])
for i = 1:nS
    cp3(1).CData(i,:) = [44,127,184]/255; % [1 1 1];
    cp3(2).CData(i,:) = [254,232,200]/255;
    cp3(3).CData(i,:) = [253,187,132]/255;
end
box on

figure(7004),
sp1.Position = [0.1300    0.3346    0.2134    0.5904];
sp2.Position = [0.4108    0.3346    0.2134    0.5904];

% % 
% set(fh_7004,'color','w')

% % % % % 
% % % % if setup.save == 1
% % % %     saveName = [savefig_loc, 'fig8d_NADHrecycle.fig'];
% % % %     savefig(fh_7004, saveName)
% % % % %     % 
% % % % %     set(gcf,'Units','inches');
% % % % %     screenposition = get(3,'Position');
% % % % %     set(gcf,...
% % % % %         'PaperPosition',[0 0 screenposition(3:4)],...
% % % % %         'PaperSize',[screenposition(3:4)]);
% % % % %     print -dpdf -painters Y3M1_robustness_rates
% % % % %     print -dpng -painters Y3M1_robustness_rates    
% % % % end


% % % % %% main plotting
% % % % close all
% % % % f1_h = openfig('fig8a_sinks.fig');
% % % % f2_h = openfig('fig8b_resp.fig');
% % % % f3_h = openfig('fig8c_tre.fig');
% % % % f3b_h = openfig('fig8c2_ixp.fig');
% % % % f4_h = openfig('fig8d_NADHrecycle.fig');
% % % % 
% % % % % 
% % % % f8001_h = figure(8001);
% % % % % f8001_h.Position = [962 42 958 954];
% % % % f8001_h.Position = [1258 42 662 954];
% % % % % % % % f8001_h.Position = [1 41 1920 963];
% % % % 
% % % % % plot all of them + locate in the right location
% % % % % prepare plots for IXP
% % % % % correct text
% % % % 
% % % % % %%
% % % % new_f1_h = copyobj(f1_h.Children(6),8001);
% % % % % new_f1_h.Position = [0.1300    0.7810    0.3347    0.1440];
% % % % % % % % new_f1_h.Position = [0.1300    0.5838    0.1566    0.3412];
% % % % new_f1_h.Position = [0.1300    0.7093    0.3347    0.2157];
% % % % 
% % % % 
% % % % 
% % % % % plotting added dilrate
% % % % for temp_plot = 1
% % % % % % % %     catSave = new_f1_h.Children(4).XData;
% % % % % % % % new_f1_h.Children(1).XData
% % % % % % % % 
% % % % % % % % ans = 
% % % % % % % % 
% % % % % % % %   1×8 categorical array
% % % % % % % % 
% % % % % % % %   Columns 1 through 5
% % % % % % % % 
% % % % % % % %      0.02      0.05      0.10      0.20      0.30 
% % % % % % % % 
% % % % % % % %   Columns 6 through 8
% % % % % % % % 
% % % % % % % %      0.325      0.35      0.375 
% % % %     catSave = new_f1_h.Children(1).XData;
% % % %     loadData_Y3M1;
% % % %     hold on
% % % %     yyaxis right
% % % % % % % % %     plot(catSave,canelas_SS.mtlD.v_GLT,'ko-','MarkerFaceColor','k')
% % % % % % % % simrate
% % % % % % % % 
% % % % % % % % simrate =
% % % % % % % % 
% % % % % % % %   Columns 1 through 5
% % % % % % % % 
% % % % % % % %     0.0448    0.0888    0.1594    0.1514    0.6554
% % % % % % % % 
% % % % % % % %   Columns 6 through 8
% % % % % % % % 
% % % % % % % %     0.7489    1.0272    1.3587
% % % %     simrate = f4_h.Children(6).Children(3).YData;
% % % %     plot(catSave,simrate,'ko-','MarkerFaceColor','k')
% % % %     ylabel('Transport Rate (mmol L-1 s-1)')
% % % % end
% % % % % set(new_f1_h,'defaultAxesColorOrder',[[0 0 0]; [0 0 0]]);
% % % % 
% % % % % % % % % %%
% % % % % % % % new_f1_h2 = copyobj(f1_h.Children(2),8001);
% % % % % % % % % new_f1_h2.Position = [0.5703    0.7810    0.3347    0.1440];
% % % % % % % % % % % % new_f1_h2.Position = [0.3361    0.5838    0.1566    0.3412];
% % % % % % % % new_f1_h2.Position = [0.5703    0.7093    0.3347    0.2157];
% % % % 
% % % % % % plotting added dilrate
% % % % % for temp_plot = 1
% % % % % %     catSave = new_f1_h.Children(4).XData;
% % % % % %     loadData_Y3M1;
% % % % %     hold on
% % % % %     yyaxis right
% % % % %     plot(catSave,canelas_SS.mtlD.v_GLT,'ko-','MarkerFaceColor','k')
% % % % %     ylabel('transport rate (mmol L^{-1} s^{-1})')
% % % % % end
% % % % 
% % % % 
% % % % 
% % % % % %%
% % % % new_f1_h3 = copyobj(f2_h.Children(2),8001);
% % % % % new_f1_h3.Position = [0.1300    0.5619    0.3347    0.1440];
% % % % % % % % new_f1_h3.Position = [0.5422    0.5838    0.1566    0.3412];
% % % % new_f1_h3.Position = [0.1300    0.4096    0.3347    0.2157];
% % % % % % % % % %% 
% % % % % % % % new_f1_h3b = copyobj(f2_h.Children(4),8001);
% % % % % % % % % new_f1_h3b.Position = [0.5703    0.5619    0.3347    0.1440];
% % % % % % % % new_f1_h3b.Position = [0.6171    0.6015    0.0756    0.1296];
% % % % % % % % % %% 
% % % % % % % % new_f1_h3b.XLabel.Visible = 'off';
% % % % % % % % new_f1_h3b.XTickLabel = [];
% % % % % % % % new_f1_h3b.YLabel.Visible = 'off';
% % % % % % % % new_f1_h3b.YTickLabel = [];
% % % % % %%
% % % % new_f1_h4 = copyobj(f4_h.Children(6),8001);
% % % % % % % % new_f1_h4.Position = [0.7484    0.5838    0.1566    0.3412];
% % % % new_f1_h4.Position = [0.5703    0.4096    0.3347    0.2157];
% % % % % % % % % 
% % % % % % % % new_f1_h4b = copyobj(f4_h.Children(2),8001);
% % % % % % % % new_f1_h4b.Position = [0.7540 0.7747 0.0742 0.1345];
% % % % % % % % new_f1_h4b.XLabel.Visible = 'off';
% % % % % % % % new_f1_h4b.XTickLabel = [];
% % % % % % % % new_f1_h4b.YLabel.Visible = 'off';
% % % % % % % % new_f1_h4b.YTickLabel = [];
% % % % % % % % new_f1_h4b.Title.Visible = 'off';
% % % % % %%
% % % % new_f1_h5 = copyobj(f3_h.Children(2),8001);
% % % % % % % % new_f1_h5.Position = [0.1300    0.1100    0.1566    0.3412];
% % % % new_f1_h5.Position = [0.1300    0.1100    0.3347    0.2157];
% % % % % 
% % % % new_f1_h6 = copyobj(f3b_h.Children,8001);
% % % % % % % % new_f1_h6.Position = [0.3361    0.1100    0.1566    0.3412];
% % % % new_f1_h6.Position = [0.5703    0.1100    0.3347    0.2157];
% % % % 
% % % % % 
% % % % f8001_h = figure(8001);
% % % % 
% % % % 
% % % % %%
% % % % % % % % print(gcf,'-dpdf','E:\OneDrive - TU Eindhoven\Documents\ch4_Y3M1\results\figures\manuscript\F8AE');



% % saving the plots
% savefig(7001,'results\F8_sinks.fig')
% savefig(7002,'results\F8_respiration.fig')
% savefig(7003,'results\F8_reductive.fig')
% savefig(7004,'results\F8_trehalose.fig')



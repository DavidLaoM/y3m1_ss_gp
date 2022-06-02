% % F2_PHYSIO.m
% This code runs and plots the steady state (ss) and glucose perturbation
% (gp) simulations and reproduces the physiological varaibles plots in the
% manuscript.
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

[simRes_reference] = simulateY3M1_separateGPSS(x, canelas_SS, data, dataset, setup);
setup.simGPdata = simRes_reference{1}.gs;
setup.simSSdata = simRes_reference{1}.ss;
setup.litParams = zeros(size(x));
setup.plotResultsMode = 0;

% recall data for plots 
fh_2 = figure(2);
fh_2_children = get(fh_2,'Children');
% 
fh_5 = figure(5);
fh_5_children = get(fh_5,'Children');
% 
fh_qO2 = fh_5_children(9);
fh_qO2_axes = findobj(fh_5_children(9),'type','axes'); 
fh_qCO2 = fh_5_children(10);


% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% plot gas exchange % % % % % % % % % % % % % % % % % % % % % % % % % % % 
if exist('fh_1000','var')
    clf(1000)
end
fh_1000 = figure(1000); 
fh_1000.Position = [1968 351 766 566];
pos_sp1 = [0.1300    0.5838    0.3347    0.3412];
pos_sp2 = [0.5703    0.5838    0.3347    0.3412];
% qO2
copy_h = copyobj(fh_5_children(9),1000); % qO2
set(copy_h,'Position',pos_sp1)
% qCO2
copy_h = copyobj(fh_5_children(10),1000); % qCO2
set(copy_h,'Position',pos_sp2)
% recalling data
qCO2_sim = fh_1000.Children(1).Children(2).YData;
qCO2_exp = fh_1000.Children(1).Children(1).YData;
qO2_sim = fh_1000.Children(2).Children(2).YData;
qO2_exp = fh_1000.Children(2).Children(1).YData;
dil_rate = fh_1000.Children(1).Children(1).XData;
% splines of simulated data
spline_period = dil_rate(1):0.005:dil_rate(end);
qCO2_sim_spline = spline(dil_rate,qCO2_sim,spline_period);
qO2_sim_spline = spline(dil_rate,qO2_sim,spline_period); 
% new plot mMs
hold off
% important figure

% fh_gasEx_mMs = subplot(2,2,3);
fh_gasEx_mMs = figure(10001);
fh_gasEx_mMs.Position = [1934 625 542 355];
%
plot(dil_rate, qCO2_sim, ...
    'k-','LineWidth',1.5)
hold on
plot(dil_rate, qO2_sim, ...
    'k--','LineWidth',2.5)
plot(dil_rate, qCO2_exp, ...
    'ko','MarkerFaceColor','k')
plot(dil_rate, qO2_exp, ...
    'ko','MarkerFaceColor','w')
%
xlabel(Y3M1_labels.mu)
ylabel(Y3M1_labels.qE_mMs)
% 
ylim([0 4])
yticks([0 1 2 3 4])
xlim([0 0.4])
%
line([0.2 0.2],[0 4],'color','black')
temp1 = text(0.4*0.25,4*0.925,'respiratory','Rotation',0,...
    'HorizontalAlignment','center');
temp2 = text(0.4*0.75,4*0.925,'respirofermentative','Rotation',0,...
    'HorizontalAlignment','center');
% 
temp_legNames = cell(1,4);
temp_legNames{1} = 'q_{CO2.sim}'; 
temp_legNames{2} = 'q_{O2.sim}';
temp_legNames{3} = 'q_{CO2.exp}';
temp_legNames{4} = 'q_{O2.exp}';
legend(temp_legNames,'Location','WestOutside');
% 
fh_gasEx_mMs.Children(1).FontSize = 10;
fh_gasEx_mMs.Children(2).FontSize = 10;
hold off


% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% Plot energetics % % % % % % % % % % % % % % % % % % % % % % % % % % % %  
if exist('fh_1001','var')
    clf(1001)
end
fh_1001 = figure(1001); % here important figures extracted
fh_1001.Position = [2660 280 768 700]; % at home
% 
pos1001_sp1 = [0.1300    0.7673    0.1566    0.1577];
pos1001_sp2 = [0.3361    0.7673    0.1566    0.1577];
pos1001_sp3 = [0.5422    0.7673    0.1566    0.1577];
pos1001_sp4 = [0.7484    0.7673    0.1566    0.1577];
% 
pos1001_sp5 = [0.1300    0.5482    0.1566    0.1577];
pos1001_sp6 = [0.3361    0.5482    0.1566    0.1577];
pos1001_sp7 = [0.5422    0.5482    0.1566    0.1577];
pos1001_sp8 = [0.7484    0.5482    0.1566    0.1577];
% 
% 
pos1001_sp914 = [0.1300    0.1100    0.3628    0.3768];
pos1001_sp1116 = [0.5422    0.1100    0.3628    0.3768];
% selecting relevant
% recall important plots
fh_eCharge = fh_5_children(1);
fh_mATP = fh_5_children(7);
fh_qATPresp = fh_5_children(8);
fh_qATPferm_vGLK = fh_2_children(42);
fh_qATPferm_vPFK = fh_2_children(40);
fh_qATPferm_vPGK = fh_2_children(35);
fh_qATPferm_vPYK = fh_2_children(32);
% replacing plots 
% fh_qATPresp
copy_h = copyobj(fh_qATPresp,1001); 
set(copy_h,'Position',pos1001_sp1)
% fh_mATP
copy_h = copyobj(fh_mATP,1001); 
set(copy_h,'Position',pos1001_sp5)
% fh_eCharge
copy_h = copyobj(fh_eCharge,1001); 
set(copy_h,'Position',pos1001_sp6)
% fh_qATPferm
copy_h = copyobj(fh_qATPferm_vGLK,1001); 
set(copy_h,'Position',pos1001_sp3)
copy_h = copyobj(fh_qATPferm_vPFK,1001);
set(copy_h,'Position',pos1001_sp4) 
copy_h = copyobj(fh_qATPferm_vPGK,1001); 
set(copy_h,'Position',pos1001_sp7)
copy_h = copyobj(fh_qATPferm_vPYK,1001); 
set(copy_h,'Position',pos1001_sp8)
% replot qATPferm
qATP_ferm_sim = - fh_qATPferm_vGLK.Children(2).YData - fh_qATPferm_vPFK.Children(2).YData + fh_qATPferm_vPGK.Children(2).YData + fh_qATPferm_vPYK.Children(2).YData; 
qATP_ferm_exp = - fh_qATPferm_vGLK.Children(1).YData - fh_qATPferm_vPFK.Children(1).YData + fh_qATPferm_vPGK.Children(1).YData + fh_qATPferm_vPYK.Children(1).YData; 
% 
subplot(4,4,2)
plot(dil_rate,qATP_ferm_sim,'k-','Linewidth',1.5)
hold on
plot(dil_rate,qATP_ferm_exp,'r+')
title('qATPferm')
hold off
% recall data
    qATP_ferm_sim_spline = spline(dil_rate,qATP_ferm_sim,spline_period);
qATP_resp_sim = fh_qATPresp.Children(2).YData; % from plot
    qATP_resp_sim_spline = spline(dil_rate,qATP_resp_sim,spline_period);
qATP_resp_exp = fh_qATPresp.Children(1).YData; % from plot
mATP_sim = fh_mATP.Children(3).YData; % from plot
    mATP_sim_spline = spline(dil_rate,mATP_sim,spline_period);
mATP_exp60 = fh_mATP.Children(2).YData; % from plot
mATP_exp40 = fh_mATP.Children(1).YData; % from plot
eCharge_sim = fh_eCharge.Children(2).YData; % from plot
    eCharge_sim_spline = spline(dil_rate,eCharge_sim,spline_period);
eCharge_exp = fh_eCharge.Children(1).YData; % from plot
% GAM % from theory, use spline_period
% NGAM % from theory
   
% %% new plot.main
% blue [0 0 1]
% red [1 0 0]
% white [1 1 1]
hold off

% important figure
if exist('fh_energetics_mMs','var')
    delete(fh_energetics_mMs)
end
% fh_energetics_mMs = subplot(4,4,[9 10 13 14]);
fh_energetics_mMs = figure(10002);
fh_energetics_mMs.Position = [1912 441 564 355];
% 
hold on
% area plots qATP
% ar_h = area(spline_period, [qATP_resp_sim_spline; qATP_ferm_sim_spline]','LineStyle','none');
ar_h = area(dil_rate, [qATP_resp_sim; qATP_ferm_sim]','LineStyle','none');
ar_h(1).FaceColor = [.5 .5 1];
ar_h(2).FaceColor = [1 .5 .5];
% m40
plot(dil_rate, mATP_exp40, ...
    'k--')
% m60
plot(dil_rate, mATP_exp60, ...
    'k--')
%
xlabel(Y3M1_labels.mu) 
ylabel(Y3M1_labels.q_mMs)
% 
ylim([-0.02 4])
yticks([0 1 2 3 4])
xlim([0 0.4])
%
line([0.2 0.2],[0 3],'color','black')
if exist('temp1b','var')
    delete(temp1b)
    delete(temp2b)
end
line([0.2 0.2],[0 4],'color','black')
temp1b = text(0.4*0.25,4*0.925,'respiratory','Rotation',0,...
    'HorizontalAlignment','center');
temp2b = text(0.4*0.75,4*0.925,'respirofermentative','Rotation',0,...
    'HorizontalAlignment','center');
% 
leg_h = legend('q_{ATP,resp,sim}','q_{ATP,glyc,sim}','m_{ATP,exp}', ...
    'Orientation','Vertical','Location','WestOutside');
temp3b = text(0.3,1.6,'GAM40','Rotation',25,...
    'HorizontalAlignment','center');
temp4b = text(0.28,2.65,'GAM60','Rotation',41,...
    'HorizontalAlignment','center');
% legend boxoff  
box on
fh_energetics_mMs.Children(1).FontSize = 10;
fh_energetics_mMs.Children(2).FontSize = 10;
% %%
hold off


% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% Plot sugar metabolism % % % % % % % % % % % % % % % % % % % % % % % % % 
% % ref C-flux GLT and sinks
% C_vGLT = vGLT * 6; 
% C_vSinks = 6 * exp_v_sinkG6P{npSets} - 6 * exp_v_sinkF6P{npSets} - 3 * exp_v_sinkGAP{npSets} + 3 * exp_v_sinkP3G{npSets} + 3 * exp_v_sinkPEP{npSets} + 3 * exp_v_sinkPYR{npSets} + 2 * exp_v_sinkACE{npSets};
fh_vGLT = fh_2_children(43);
fh_vsinkG6P = fh_2_children(10);
fh_vsinkF6P = fh_2_children(9);
fh_vsinkGAP = fh_2_children(8);
fh_vsinkP3G = fh_2_children(7);
fh_vsinkPEP = fh_2_children(6);
fh_vsinkPYR = fh_2_children(5);
fh_vsinkACE = fh_2_children(4);
% 
C_vGLT_sim = fh_vGLT.Children(2).YData * 6; 
    C_vGLT_sim_spline = spline(dil_rate([1:5,7:8]),C_vGLT_sim([1:5,7:8]),spline_period);
%     C_vGLT_sim_spline = spline(dil_rate([1:8]),C_vGLT_sim([1:8]),spline_period);
C_vGLT_exp = fh_vGLT.Children(1).YData * 6; 
% 
C_vsinks_sim = 6 * fh_vsinkG6P.Children(2).YData - 6 * fh_vsinkF6P.Children(2).YData - 3 * fh_vsinkGAP.Children(2).YData + 3 * fh_vsinkP3G.Children(2).YData + 3 * fh_vsinkPEP.Children(2).YData + 3 * fh_vsinkPYR.Children(2).YData + 2 * fh_vsinkACE.Children(2).YData; 
    C_vsinks_sim_spline = spline(dil_rate,C_vsinks_sim,spline_period);
C_vsinks_exp = 6 * fh_vsinkG6P.Children(1).YData - 6 * fh_vsinkF6P.Children(1).YData - 3 * fh_vsinkGAP.Children(1).YData + 3 * fh_vsinkP3G.Children(1).YData + 3 * fh_vsinkPEP.Children(1).YData + 3 * fh_vsinkPYR.Children(1).YData + 2 * fh_vsinkACE.Children(1).YData; 
% main plot
if exist('fh_1002','var')
    clf(fh_1002)
end
fh_1002 = figure(1002);
fh_1002.Position = [1923 319 575 345]; % at home
% fh_1002.Position = [100 100 560*0.75 460*0.75]; % at tu/e
hold on
%
% th(1) = area(spline_period,C_vGLT_sim_spline,'FaceColor',[.75 .75 1],'LineStyle','none');
th(1) = area(dil_rate,C_vGLT_sim,'FaceColor',[.75 .75 1],'LineStyle','none');
% th(2) = area(spline_period,C_vsinks_sim_spline,'FaceColor',[1 .75 .75],'LineStyle','none');
th(2) = area(dil_rate,C_vsinks_sim,'FaceColor',[1 .75 .75],'LineStyle','none');
%
% plot(spline_period,C_vGLT_sim_spline,'b-')
plot(dil_rate,C_vGLT_sim,'b-')
% plot(spline_period,C_vsinks_sim_spline,'r-')
plot(dil_rate,C_vsinks_sim,'r-')
%
th(3) = plot(dil_rate,C_vGLT_exp,'ko','MarkerFaceColor','b');
plot(dil_rate,C_vsinks_exp,'ko','MarkerFaceColor','r')
% 
xlabel(Y3M1_labels.mu) 
ylabel(Y3M1_labels.qC_mMs)
%
ylim([-0.02 10])
yticks([0 2.5 5 7.5 10])
xlim([0 0.4])
xticks([0 0.1 0.2 0.3 0.4])
%
line([0.2 0.2],[0 10],'color','black')
if exist('temp1c','var')
    delete(temp1c)
    delete(temp2c)
end
line([0.2 0.2],[0 4],'color','black')
temp1c = text(0.4*0.25,10*0.925,'respiratory','Rotation',0,...
    'HorizontalAlignment','center');
temp2c = text(0.4*0.75,10*0.925,'respirofermentative','Rotation',0,...
    'HorizontalAlignment','center');
%
[leg_h,icons_h] = legend(th,'q_{C,GLT}','q_{C,Sinks}','Experimental', ...
    'Orientation','Vertical','location','WestOutside');%,'Location','SouthOutside');
% leg_h.Box = 'off';
icons_h(7).MarkerFaceColor = 'none';
% 
box on
hold off


% % saving the plots
% savefig(10002,'results\F2_energetics.fig')
% savefig(1002,'results\F2_c_flux.fig')
% savefig(10001,'results\F2_gas_exchange.fig')


%% (to be located outside in the end) Initial setup:
robust_setup.simulation_reference = refSim;
robust_setup.simulation_robustness = simRobust;
%     robust_setup.simulation_robustness = selResIntial;
%     for i = 1:39
%         robust_setup.simulation_robustness = [robust_setup.simulation_robustness, selResIntial];
%     end
robust_setup.forced_density = 0;
% robust_setup.simulation_robustness_idxs = 1:10; 
% robust_setup.simulation_robustness_idxs = 1:500; %1:100 %1:40
robust_setup.simulation_robustness_idxs = 1:1000; %1:100 %1:40
% robust_setup.simulation_robustness_idxs = 1:10000; %1:100 %1:40
% robust_setup.simulation_robustness_idxs = 1:3000; %1:100 %1:40
% % % % robust_setup.simulation_robustness_idxs = 1:100; %1:40
robust_setup.interpSpace = [0:340]';% = [0:400]';
% robust_setup.dpVal = 256;
robust_setup.dpVal = 340; %400;
% robust_setup.dpVal_yaxis = 40;
robust_setup.dpVal_yaxis = 75;
% robust_setup.dpVal_yaxis = 400;
robust_setup.maxRange_time = 340;% 400;
robust_setup.LineWidth = 2;
robust_setup.MarkerSize = 30;
robust_setup.contourLevels = 500; %32
robust_setup.Ylim_foldIncrease = 3;
robust_setup.shading_increase_ratio = 0; % 0.1

% GP: metabolite-specific
% % % % robust_setup.metabolite_id = 15;
% % % % robust_setup.maxRange_concentration = 5;
robust_setup.nMets = 32; % 38
robust_setup.metabolite_id_range = 1:32; %1:38; %1:15; % 1:38
robust_setup.maxRange_concentration_range = zeros(1,robust_setup.nMets);
% loop to get the ymax from experimental data and, if not, from the
% simulation
f3_h = figure(3);
f3_c_h = get(f3_h,'children');
% 
for i = 1:length(robust_setup.maxRange_concentration_range)
    % option based on simulated axes elsewhere
    robust_setup.maxRange_concentration_range(i) = f3_c_h(end+1-i).YLim(2);
%     Another option
%     if isfield(ExpData.metabolites{i}, 'conc')
%         robust_setup.maxRange_concentration_range(i) = max(ExpData.metabolites{i}.conc)*robust_setup.Ylim_foldIncrease;
%     else
%         robust_setup.maxRange_concentration_range(i) = max(selResIntial{1}.Y_FF01(:,i))*robust_setup.Ylim_foldIncrease;
%     end
end
% 
robust_setup.maxRange_concentration_range(1) = 0.5; % ACE
robust_setup.maxRange_concentration_range(2) = 0.003; % BPG
robust_setup.maxRange_concentration_range(3) = 12; % FBP
robust_setup.maxRange_concentration_range(4) = 2.5; % F6P
robust_setup.maxRange_concentration_range(5) = 15; % G6P
robust_setup.maxRange_concentration_range(6) = 30; % 
robust_setup.maxRange_concentration_range(7) = 1.6; % 
robust_setup.maxRange_concentration_range(8) = 0.025; % 
% 
robust_setup.maxRange_concentration_range(9) = 3; % 
robust_setup.maxRange_concentration_range(10) = 0.25; % 
robust_setup.maxRange_concentration_range(11) = 3; % 
robust_setup.maxRange_concentration_range(12) = 1.5; % 
robust_setup.maxRange_concentration_range(13) = 2.25; % 
robust_setup.maxRange_concentration_range(14) = 0.1; % 
robust_setup.maxRange_concentration_range(15) = 1.5; % 
robust_setup.maxRange_concentration_range(16) = 0.65; % 
% 
robust_setup.maxRange_concentration_range(17) = 6; % 
robust_setup.maxRange_concentration_range(18) = 0.45; % 
robust_setup.maxRange_concentration_range(19) = 0.35; % 
robust_setup.maxRange_concentration_range(20) = 12.5; % 
robust_setup.maxRange_concentration_range(21) = 0.5; % 
robust_setup.maxRange_concentration_range(22) = 2; % 
robust_setup.maxRange_concentration_range(23) = 1.5; % 
robust_setup.maxRange_concentration_range(24) = 0.1; % 
% 
robust_setup.maxRange_concentration_range(25) = 90; % 
robust_setup.maxRange_concentration_range(26) = 15; % 
robust_setup.maxRange_concentration_range(27) = 30; % 
robust_setup.maxRange_concentration_range(28) = 1; % 
robust_setup.maxRange_concentration_range(29) = 2.25; % 
robust_setup.maxRange_concentration_range(30) = 2; % 
robust_setup.maxRange_concentration_range(31) = 0.75; % 
robust_setup.maxRange_concentration_range(32) = 1; % 
% 
legenda.metabolites{3} = 'FBP, x_{3}, [mM]';
legenda.metabolites{14} = 'GAP, x_{14}, [mM]';
legenda.metabolites{18} = 'G3P, x_{18}, [mM]';
legenda.metabolites{19} = 'GLYC, x_{19}, [mM]';

% % % % % troubleshooting
% % % % robust_setup.maxRange_concentration_range = 2 * robust_setup.maxRange_concentration_range;

% % Specific changes, upon a given case
% robust_setup.maxRange_concentration_range(9) = 8;
% robust_setup.maxRange_concentration_range(15) = 2.5;
% robust_setup.maxRange_concentration_range(38) = 80;


% GP: reaction rate-specific
robust_setup.nRates = 42; %51;
robust_setup.rate_id_range = 1:42;% = 1:51;
robust_setup.maxRange_rate_range = zeros(1,robust_setup.nRates);
% loop to get the ymax from experimental data and, if not, from the
% simulation
f4_h = figure(4);
f4_c_h = get(f4_h,'children');
% 
for i = 1:length(robust_setup.maxRange_rate_range)
    % option based on simulated axes elsewhere
    robust_setup.maxRange_rate_range(i) = f4_c_h(end+1-i).YLim(2);
end
robust_setup.maxRange_rate_range(1) = 10;
robust_setup.maxRange_rate_range(2) = 3;
robust_setup.maxRange_rate_range(3) = 1;
robust_setup.maxRange_rate_range(4) = 1;
robust_setup.maxRange_rate_range(5) = 0.6;
robust_setup.maxRange_rate_range(6) = 0.06;
robust_setup.maxRange_rate_range(7) = 0.6;
robust_setup.maxRange_rate_range(8) = 1;
robust_setup.maxRange_rate_range(9) = 1;
% 
robust_setup.maxRange_rate_range(10) = 1;
robust_setup.maxRange_rate_range(11) = 1;
robust_setup.maxRange_rate_range(12) = 1;
robust_setup.maxRange_rate_range(13) = 1;
robust_setup.maxRange_rate_range(14) = 0.5;
robust_setup.maxRange_rate_range(15) = 0.05;
robust_setup.maxRange_rate_range(16) = 1;
robust_setup.maxRange_rate_range(17) = 0.25;
robust_setup.maxRange_rate_range(18) = 1;
% 
robust_setup.maxRange_rate_range(19) = 0.1;
robust_setup.maxRange_rate_range(20) = 0.05;
robust_setup.maxRange_rate_range(21) = 0.25;
robust_setup.maxRange_rate_range(22) = 0.1;
robust_setup.maxRange_rate_range(23) = 0.5;
robust_setup.maxRange_rate_range(24) = 0.05;
robust_setup.maxRange_rate_range(25) = 1;
robust_setup.maxRange_rate_range(26) = 0.5;
robust_setup.maxRange_rate_range(27) = 1.25;
% 
robust_setup.maxRange_rate_range(28) = 1.5;
robust_setup.maxRange_rate_range(29) = 1;
robust_setup.maxRange_rate_range(30) = 0.5;
robust_setup.maxRange_rate_range(31) = 0.25;
robust_setup.maxRange_rate_range(32) = 0.05;
robust_setup.maxRange_rate_range(33) = 0.2;
% robust_setup.maxRange_rate_range(34) = 1;
% robust_setup.maxRange_rate_range(35) = 1;
% robust_setup.maxRange_rate_range(36) = 1;
% % 
% robust_setup.maxRange_rate_range(37) = 1;
% robust_setup.maxRange_rate_range(38) = 1;
% robust_setup.maxRange_rate_range(39) = 1;
% robust_setup.maxRange_rate_range(40) = 1;
% robust_setup.maxRange_rate_range(41) = 1;
% robust_setup.maxRange_rate_range(42) = 1;



% %% Figure for metabolites
robust_setup.plotSS = 1;
if exist('fig_h_c_gs','var')
    clear fig_h_c_gs
end
fig_h_c_gs = figure(2003);
fig_h_c_gs.Position = [1921 257 1536 747];
for j = robust_setup.metabolite_id_range
    if j ~= 1
        toc
    end
    disp(j)
    tic
    % minimum setup each loop
    robust_setup.metabolite_id = j;
    robust_setup.maxRange_concentration = robust_setup.maxRange_concentration_range(j);
    % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
    % %% arrange data
    % heatmap data
    if exist('tData','var')
        clear tData yData
    end
    tData = [];
    yData = [];
% % % %     for i = 1:robust_setup.simulation_robustness_idxs
    for i = robust_setup.simulation_robustness_idxs
        tData = [tData; robust_setup.interpSpace];
% % % %         yData = [yData; interp1(robust_setup.simulation_robustness{i}.T_FF01, robust_setup.simulation_robustness{i}.Y_FF01(:,15), robust_setup.interpSpace, 'pchip')];
        yData = [yData; interp1(robust_setup.simulation_robustness{i}.gs.T, robust_setup.simulation_robustness{i}.gs.Y(:,j), robust_setup.interpSpace, 'pchip')];
    end
    tData = [0; robust_setup.maxRange_time; tData];
    yData = [0; robust_setup.maxRange_concentration; yData];

    % simulation data
    sim_time = robust_setup.simulation_reference{1}.gs.T(3002:end) * robust_setup.dpVal / max(tData); %[1 2 10 100];
% % % %     sim_conc = robust_setup.dpVal - robust_setup.simulation_reference{1}.Y_FF01(:,15) * robust_setup.dpVal / max(yData);
    sim_conc = robust_setup.dpVal - robust_setup.simulation_reference{1}.gs.Y(3002:end,j) * robust_setup.dpVal / max(yData);

    % experimental data
%     if isfield(ExpData.metabolites{j}, 'conc')
    if data.ordered.metabolites.onoff(j) == 1
%         exp_time = ExpData.metabolites{robust_setup.metabolite_id}.time * robust_setup.dpVal / max(tData);
        exp_time = data.ordered.metabolites.time{robust_setup.metabolite_id} * robust_setup.dpVal / max(tData);
%         exp_conc = robust_setup.dpVal - ExpData.metabolites{robust_setup.metabolite_id}.conc * robust_setup.dpVal / max(yData);
        exp_conc = robust_setup.dpVal - data.ordered.metabolites.data{robust_setup.metabolite_id} * robust_setup.dpVal / max(yData);
    end
    % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
    % %% Visualization (in loop)
    % % % % if exist('fig_h_c_gs','var')
    % % % %     clear fig_h_c_gs
    % % % % end
    % % % % fig_h_c_gs = figure(1004);
    figure(fig_h_c_gs);

%     if j == 26
%         disp('Stop here.')
%     end

    % make the heatmap plot
%     sp_h = subplot(4,8,robust_setup.metabolite_id);
    % late change
    if robust_setup.metabolite_id == 1, tempVal = 6;
    elseif robust_setup.metabolite_id == 2, tempVal = 7;
    elseif robust_setup.metabolite_id == 3, tempVal = 3;
    elseif robust_setup.metabolite_id == 4, tempVal = 2;
    elseif robust_setup.metabolite_id == 5, tempVal = 1;
    elseif robust_setup.metabolite_id == 6, tempVal = 8;
    elseif robust_setup.metabolite_id == 7, tempVal = 14;
    elseif robust_setup.metabolite_id == 8, tempVal = 15;
    % 
    elseif robust_setup.metabolite_id == 9, tempVal = 17;
    elseif robust_setup.metabolite_id == 10, tempVal = 10;
    elseif robust_setup.metabolite_id == 11, tempVal = 9;
    elseif robust_setup.metabolite_id == 12, tempVal = 11;
    elseif robust_setup.metabolite_id == 13, tempVal = 12;
    elseif robust_setup.metabolite_id == 14, tempVal = 4;
    elseif robust_setup.metabolite_id == 15, tempVal = 18;
    elseif robust_setup.metabolite_id == 16, tempVal = 19;
    % 
    elseif robust_setup.metabolite_id == 17, tempVal = 5;
    elseif robust_setup.metabolite_id == 18, tempVal = 16;
    elseif robust_setup.metabolite_id == 19, tempVal = 20;
    elseif robust_setup.metabolite_id == 20, tempVal = 13;
    elseif robust_setup.metabolite_id == 21, tempVal = 25;
    elseif robust_setup.metabolite_id == 22, tempVal = 21;
    elseif robust_setup.metabolite_id == 23, tempVal = 22;
    elseif robust_setup.metabolite_id == 24, tempVal = 23;
    % 
    elseif robust_setup.metabolite_id == 25, tempVal = 27;
    elseif robust_setup.metabolite_id == 26, tempVal = 26;
    elseif robust_setup.metabolite_id == 27, tempVal = 24;
    elseif robust_setup.metabolite_id == 28, tempVal = 28;
    elseif robust_setup.metabolite_id == 29, tempVal = 29;
    elseif robust_setup.metabolite_id == 30, tempVal = 30;
    elseif robust_setup.metabolite_id == 31, tempVal = 31;
    elseif robust_setup.metabolite_id == 32, tempVal = 32;
    end
    sp_h = subplot(4,8,tempVal);
    
    
%     [f_h] = DataDensityPlot(tData, yData, robust_setup.contourLevels, fig_h_c_gs, sp_h);
    [f_h, range90] = DataDensityPlot_1D_Y3M1(tData, yData, robust_setup.contourLevels, fig_h_c_gs, sp_h, robust_setup);
    hold on

    % simulation data
% % % %     plot(sim_time, sim_conc,'k-','LineWidth', robust_setup.LineWidth)
% % % %     plot(sim_time, sim_conc/(max(sim_conc)-min(sim_conc))*robust_setup.dpVal_yaxis,'k-','LineWidth', robust_setup.LineWidth)
% % % %     plot(sim_time, sim_conc/max(sim_conc)*robust_setup.dpVal_yaxis,'k-','LineWidth', robust_setup.LineWidth)
    sp_h.XLim(1) = 0;
    plot(sim_time, sim_conc/robust_setup.dpVal*robust_setup.dpVal_yaxis,'k-','LineWidth', robust_setup.LineWidth)

%     % 
%     figure,
%     subplot(131), plot(range90(1,:),'.-')
%     subplot(132), plot(range90(2,:),'.-')
%     subplot(133), plot(range90(3,:),'.-')
    temp_range90min = robust_setup.dpVal - range90(2,:) * robust_setup.dpVal / max(yData);
    temp_range90max = robust_setup.dpVal - range90(3,:) * robust_setup.dpVal / max(yData);
    plot(range90(1,:), temp_range90min/robust_setup.dpVal*robust_setup.dpVal_yaxis,'k--','LineWidth', 1)
    plot(range90(1,:), temp_range90max/robust_setup.dpVal*robust_setup.dpVal_yaxis,'k--','LineWidth', 1)
    
% % % %     sim_conc = robust_setup.dpVal - robust_setup.simulation_reference{1}.Y_FF01(:,j) * robust_setup.dpVal / max(yData);
    
    
    % experimental data
%     if isfield(ExpData.metabolites{j}, 'conc')
    if data.ordered.metabolites.onoff(j) == 1
% % % %         scatter(exp_time, exp_conc, ...
        scatter(exp_time, exp_conc/robust_setup.dpVal*robust_setup.dpVal_yaxis, ...
                    robust_setup.MarkerSize, 'MarkerEdgeColor', 'black', ...
                    'MarkerFaceColor', 'white', ...
                    'LineWidth', 1.5)
    end
    % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
    % add text label
    metName_pre = erase(legenda.metabolites{j}, ", [mM]");
    metName = erase(metName_pre(1:end-7),",");
    % exceptions
%     if j == 14
%         metName = 'GAP';
%     elseif j == 19
%         metName = 'GLYC';
%     end
    text(sp_h.XLim(2)*0.9,sp_h.YLim(2)*0.15,metName,'FontSize',10,...
        'HorizontalAlignment','right')

end       


% %% Figure for reaction rates
robust_setup.maxRange_concentration = 0; % for plots of rates, rathe than concentrations
if exist('fig_h_r_gs','var')
    clear fig_h_r_gs
end
% %%
fig_h_r_gs = figure(2004);
fig_h_r_gs.Position = [1921 257 1536 747];
for j = robust_setup.rate_id_range
    % minimum setup each loop
    robust_setup.rate_id = j;
    robust_setup.maxRange_rate = robust_setup.maxRange_rate_range(j);
    % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
    % %% arrange data
    % heatmap data
    if exist('tData','var')
        clear tData yData
    end
    tData = [];
    yData = [];
% % % %     for i = 1:robust_setup.simulation_robustness_idxs
    for i = robust_setup.simulation_robustness_idxs
        tData = [tData; robust_setup.interpSpace];
% % % %         yData = [yData; interp1(robust_setup.simulation_robustness{i}.T_FF01, robust_setup.simulation_robustness{i}.Y_FF01(:,15), robust_setup.interpSpace, 'pchip')];
        yData = [yData; interp1(robust_setup.simulation_robustness{i}.gs.T, robust_setup.simulation_robustness{i}.gs.V(:,j), robust_setup.interpSpace, 'pchip')];
    end
    tData = [0; robust_setup.maxRange_time; tData];
    yData = [0; robust_setup.maxRange_rate; yData];

    % simulation data
    sim_time = robust_setup.simulation_reference{1}.gs.T(3002:end) * robust_setup.dpVal / max(tData); %[1 2 10 100];
% % % %     sim_conc = robust_setup.dpVal - robust_setup.simulation_reference{1}.Y_FF01(:,15) * robust_setup.dpVal / max(yData);
    sim_conc = robust_setup.dpVal - robust_setup.simulation_reference{1}.gs.V(3002:end,j) * robust_setup.dpVal / max(yData);
    
    % experimental data
%     if isfield(ExpData.fluxes{j}, 'rate')
    if data.ordered.fluxes.onoff(j) == 1
        exp_time = data.ordered.fluxes.time{robust_setup.rate_id} * robust_setup.dpVal / max(tData); exp_time(1) = 0.5;
        exp_flux = robust_setup.dpVal - data.ordered.fluxes.data{robust_setup.rate_id} * robust_setup.dpVal / max(yData);
    end
    
    % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
    % %% Visualization (in loop)
    % % % % if exist('fig_h','var')
    % % % %     clear fig_h
    % % % % end
    % % % % fig_h = figure(1004);
    figure(fig_h_r_gs);

    % make the heatmap plot
    sp_h2 = subplot(5,9,robust_setup.rate_id);
%     [f_h] = DataDensityPlot(tData, yData, robust_setup.contourLevels, fig_h_r_gs, sp_h2);
    [f_h, range90] = DataDensityPlot_1D_Y3M1(tData, yData, robust_setup.contourLevels, fig_h_r_gs, sp_h2, robust_setup);
    hold on

    % simulation data
% % % %     plot(sim_time, sim_conc,'k-','LineWidth', robust_setup.LineWidth)
    plot(sim_time, sim_conc/robust_setup.dpVal*robust_setup.dpVal_yaxis,'k-','LineWidth', robust_setup.LineWidth)

    % 90% range
    temp_range90min = robust_setup.dpVal - range90(2,:) * robust_setup.dpVal / max(yData);
    temp_range90max = robust_setup.dpVal - range90(3,:) * robust_setup.dpVal / max(yData);
    plot(range90(1,:), temp_range90min/robust_setup.dpVal*robust_setup.dpVal_yaxis,'k--','LineWidth', 1)
    plot(range90(1,:), temp_range90max/robust_setup.dpVal*robust_setup.dpVal_yaxis,'k--','LineWidth', 1)
    
    % experimental data
%     if isfield(ExpData.fluxes{j}, 'rate')
    if data.ordered.fluxes.onoff(j) == 1
% % % %         scatter(exp_time, exp_conc, ...
        scatter(exp_time, exp_flux/robust_setup.dpVal*robust_setup.dpVal_yaxis, ...
                    robust_setup.MarkerSize, 'MarkerEdgeColor', 'black', ...
                    'MarkerFaceColor', 'white', ...
                    'LineWidth', 1.5)
    end
    % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
    % add text label
    rateName_pre = erase(legenda.fluxes{j}, ", [mM s^{-1}]");
    rateName = erase(rateName_pre(1:end-7),",");
%     text(sp_h.XLim(2)*0.7,sp_h.YLim(2)*0.25,rateName,'FontSize',10)
    % 
    text(sp_h.XLim(2)*0.9,sp_h.YLim(2)*0.15,rateName,'FontSize',10,...
        'HorizontalAlignment','right')
    % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
end       
robust_setup.plotSS = 0;


%% SS simulations

% 
robust_setup.interpSpace = [0.02:0.005:0.375]';% = [0:400]';
robust_setup.dpVal = length(robust_setup.interpSpace); %400;
robust_setup.maxRange_time = 0.375;
% % % % %
% % % % robust_setup.forcedSSdata = 1;
% % % % robust_setup.forcedSSdata_vals = [3    0.0200    10.0000    3.000   15.0000,...
% % % %     4.0000    3.000    0.300    6.0000    1,...
% % % %     10.0000    3.0000    10.0000    0.100    3.0000,...
% % % %     0.6000    2    0.8000    0.8000   90.0000,...
% % % %     0.5000    2.0000    2.0000    2.0000  150.0000,...
% % % %     0.2000   11.0000    2.0000    2.0000    3.0000,...
% % % %     25.0000    0.75000];

%% Figure for SS metabolites
% simulation
f1_h = figure(1);
f1_c_h = get(f1_h,'children');
% 
for i = 1:length(robust_setup.maxRange_concentration_range)
    % option based on simulated axes elsewhere
    robust_setup.maxRange_concentration_range(i) = f1_c_h(end+1-i).YLim(2);
end
% %
robust_setup.maxRange_concentration_range(1) = 23;
robust_setup.maxRange_concentration_range(2) = 0.03;
robust_setup.maxRange_concentration_range(3) = 8;
robust_setup.maxRange_concentration_range(4) = 2;
robust_setup.maxRange_concentration_range(5) = 14;
robust_setup.maxRange_concentration_range(6) = 3;
% robust_setup.maxRange_concentration_range(7) = 1.6;
% robust_setup.maxRange_concentration_range(8) = 1;
% %
robust_setup.maxRange_concentration_range(9) = 5;
robust_setup.maxRange_concentration_range(10) = 1;
robust_setup.maxRange_concentration_range(11) = 12.5;
robust_setup.maxRange_concentration_range(12) = 3;
robust_setup.maxRange_concentration_range(13) = 7;
robust_setup.maxRange_concentration_range(14) = 0.05;
robust_setup.maxRange_concentration_range(15) = 2.5;
robust_setup.maxRange_concentration_range(16) = 2;
% 
robust_setup.maxRange_concentration_range(17) = 1.25;
robust_setup.maxRange_concentration_range(18) = 2;
robust_setup.maxRange_concentration_range(19) = 1.25;
robust_setup.maxRange_concentration_range(20) = 70;
% robust_setup.maxRange_concentration_range(21) = 1;
% robust_setup.maxRange_concentration_range(22) = 1;
% robust_setup.maxRange_concentration_range(23) = 1;
% robust_setup.maxRange_concentration_range(24) = 1;
% % 
% robust_setup.maxRange_concentration_range(25) = 1;
% robust_setup.maxRange_concentration_range(26) = 1;
% robust_setup.maxRange_concentration_range(27) = 1;
% robust_setup.maxRange_concentration_range(28) = 1;
% robust_setup.maxRange_concentration_range(29) = 1;
% robust_setup.maxRange_concentration_range(30) = 1;
% robust_setup.maxRange_concentration_range(31) = 1;
% robust_setup.maxRange_concentration_range(32) = 1;

% plotting
if exist('fig_h_c_ss','var')
    clear fig_h_c_ss
end
fig_h_c_ss = figure(2001);
fig_h_c_ss.Position = [1921 257 1536 747];
for j = robust_setup.metabolite_id_range
% %     if j == 3
% %         robust_setup.specialcase = 13; % PFK_SS range correction
% %     else
% %         robust_setup.specialcase = 0;
% %     end
    % 
    if j ~= 1
        toc
    end
    disp(j)
    tic
    % minimum setup each loop
    robust_setup.metabolite_id = j;
    robust_setup.maxRange_concentration = robust_setup.maxRange_concentration_range(j);
% % % %     robust_setup.maxRange_concentration_forced = robust_setup.forcedSSdata_vals(j);
    % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
    % %% arrange data
    % heatmap data
    if exist('tData','var')
        clear tData yData
    end
    tData = [];
    yData = [];
% % % %     for i = 1:robust_setup.simulation_robustness_idxs
    for i = robust_setup.simulation_robustness_idxs
%         tData = [tData; robust_setup.interpSpace];
        tData = [tData; [0:71]'];
% % % %         yData = [yData; interp1(robust_setup.simulation_robustness{i}.T_FF01, robust_setup.simulation_robustness{i}.Y_FF01(:,15), robust_setup.interpSpace, 'pchip')];
        yData = [yData; interp1(canelas_SS.mtlD.D, robust_setup.simulation_robustness{i}.ss.Y(:,j)', robust_setup.interpSpace, 'pchip')];
%         % 
%         figure, plot(tData,yData)
    end
    tData = [0; robust_setup.maxRange_time; tData];
    yData = [0; robust_setup.maxRange_concentration; yData];
    % 
    if j == 3 % PFK_SS range correction
        temp_maxVal_y = 7;
        indices = find(abs(yData)>temp_maxVal_y);
        yData(indices) = temp_maxVal_y;
    elseif j == 14 % GAP_SS range correction
        temp_maxVal_y = 0.04;
        indices = find(abs(yData)>temp_maxVal_y);
        yData(indices) = temp_maxVal_y;
    elseif j == 17 % DHAP_SS range correction
        temp_maxVal_y = 1;
        indices = find(abs(yData)>temp_maxVal_y);
        yData(indices) = temp_maxVal_y;
    else
    end
    
    
    
    % simulation data
%     sim_time = robust_setup.simulation_reference{1}.gs.T(3002:end) * robust_setup.dpVal / max(tData); %[1 2 10 100];
    sim_time = canelas_SS.mtlD.D * robust_setup.dpVal / max(tData); %[1 2 10 100];
% % % %     sim_conc = robust_setup.dpVal - robust_setup.simulation_reference{1}.Y_FF01(:,15) * robust_setup.dpVal / max(yData);
%     sim_conc = robust_setup.dpVal - robust_setup.simulation_reference{1}.gs.Y(3002:end,j) * robust_setup.dpVal / max(yData);
    temp_val = zeros(8,1);
        temp_val(1) = robust_setup.simulation_reference{1}.ss.Y.ss2can{1}(end,j);
        temp_val(2) = robust_setup.simulation_reference{1}.ss.Y.ss2can{2}(end,j);
        temp_val(3) = robust_setup.simulation_reference{1}.ss.Y.ss2can{3}(end,j);
        temp_val(4) = robust_setup.simulation_reference{1}.ss.Y.ss2can{4}(end,j);
        temp_val(5) = robust_setup.simulation_reference{1}.ss.Y.ss2can{5}(end,j);
        temp_val(6) = robust_setup.simulation_reference{1}.ss.Y.ss2can{6}(end,j);
        temp_val(7) = robust_setup.simulation_reference{1}.ss.Y.ss2can{7}(end,j);
        temp_val(8) = robust_setup.simulation_reference{1}.ss.Y.ss2can{8}(end,j);
    sim_conc = robust_setup.dpVal - temp_val * robust_setup.dpVal / max(yData);

    % experimental data
%     if isfield(ExpData.metabolites{j}, 'conc')
    if canelas_SS.ordered.metabolites.onoff(j) == 1
%         exp_time = ExpData.metabolites{robust_setup.metabolite_id}.time * robust_setup.dpVal / max(tData);
        exp_time = canelas_SS.mtlD.D * robust_setup.dpVal / max(tData);
%         exp_conc = robust_setup.dpVal - ExpData.metabolites{robust_setup.metabolite_id}.conc * robust_setup.dpVal / max(yData);
        exp_conc = robust_setup.dpVal - canelas_SS.ordered.metabolites.data{robust_setup.metabolite_id} * robust_setup.dpVal / max(yData);
    end
    % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
    % %% Visualization (in loop)
    % % % % if exist('fig_h_c_gs','var')
    % % % %     clear fig_h_c_gs
    % % % % end
    % % % % fig_h_c_gs = figure(1004);
    figure(fig_h_c_ss);

%     if j == 26
%         disp('Stop here.')
%     end

    % make the heatmap plot
%     sp_h = subplot(4,8,robust_setup.metabolite_id);
    % late change
    if robust_setup.metabolite_id == 1, tempVal = 6;
    elseif robust_setup.metabolite_id == 2, tempVal = 7;
    elseif robust_setup.metabolite_id == 3, tempVal = 3;
    elseif robust_setup.metabolite_id == 4, tempVal = 2;
    elseif robust_setup.metabolite_id == 5, tempVal = 1;
    elseif robust_setup.metabolite_id == 6, tempVal = 8;
    elseif robust_setup.metabolite_id == 7, tempVal = 14;
    elseif robust_setup.metabolite_id == 8, tempVal = 15;
    % 
    elseif robust_setup.metabolite_id == 9, tempVal = 17;
    elseif robust_setup.metabolite_id == 10, tempVal = 10;
    elseif robust_setup.metabolite_id == 11, tempVal = 9;
    elseif robust_setup.metabolite_id == 12, tempVal = 11;
    elseif robust_setup.metabolite_id == 13, tempVal = 12;
    elseif robust_setup.metabolite_id == 14, tempVal = 4;
    elseif robust_setup.metabolite_id == 15, tempVal = 18;
    elseif robust_setup.metabolite_id == 16, tempVal = 19;
    % 
    elseif robust_setup.metabolite_id == 17, tempVal = 5;
    elseif robust_setup.metabolite_id == 18, tempVal = 16;
    elseif robust_setup.metabolite_id == 19, tempVal = 20;
    elseif robust_setup.metabolite_id == 20, tempVal = 13;
    elseif robust_setup.metabolite_id == 21, tempVal = 25;
    elseif robust_setup.metabolite_id == 22, tempVal = 21;
    elseif robust_setup.metabolite_id == 23, tempVal = 22;
    elseif robust_setup.metabolite_id == 24, tempVal = 23;
    % 
    elseif robust_setup.metabolite_id == 25, tempVal = 27;
    elseif robust_setup.metabolite_id == 26, tempVal = 26;
    elseif robust_setup.metabolite_id == 27, tempVal = 24;
    elseif robust_setup.metabolite_id == 28, tempVal = 28;
    elseif robust_setup.metabolite_id == 29, tempVal = 29;
    elseif robust_setup.metabolite_id == 30, tempVal = 30;
    elseif robust_setup.metabolite_id == 31, tempVal = 31;
    elseif robust_setup.metabolite_id == 32, tempVal = 32;
    end
    sp_h = subplot(4,8,tempVal);
    
    
    
    
%     [f_h] = DataDensityPlot(tData, yData, robust_setup.contourLevels, fig_h_c_gs, sp_h);
    [f_h, range90] = DataDensityPlot_1D_Y3M1(tData, yData, robust_setup.contourLevels, fig_h_c_ss, sp_h, robust_setup);
    hold on
    set(gca, 'XTickLabel', [0 0.4]);
    % simulation data
% % % %     plot(sim_time, sim_conc,'k-','LineWidth', robust_setup.LineWidth)
% % % %     plot(sim_time, sim_conc/(max(sim_conc)-min(sim_conc))*robust_setup.dpVal_yaxis,'k-','LineWidth', robust_setup.LineWidth)
% % % %     plot(sim_time, sim_conc/max(sim_conc)*robust_setup.dpVal_yaxis,'k-','LineWidth', robust_setup.LineWidth)
%     plot(sim_time*robust_setup.dpVal_yaxis/0.4, sim_conc/robust_setup.dpVal*robust_setup.dpVal_yaxis,'k-','LineWidth', robust_setup.LineWidth)

%     % 
%     figure,
%     subplot(131), plot(range90(1,:),'.-')
%     subplot(132), plot(range90(2,:),'.-')
%     subplot(133), plot(range90(3,:),'.-')
    temp_range90min = robust_setup.dpVal - range90(2,:) * robust_setup.dpVal / max(yData);
    temp_range90max = robust_setup.dpVal - range90(3,:) * robust_setup.dpVal / max(yData);
    plot(range90(1,:), temp_range90min/robust_setup.dpVal*robust_setup.dpVal_yaxis,'k--','LineWidth', 1)
    plot(range90(1,:), temp_range90max/robust_setup.dpVal*robust_setup.dpVal_yaxis,'k--','LineWidth', 1)
    
% % % %     sim_conc = robust_setup.dpVal - robust_setup.simulation_reference{1}.Y_FF01(:,j) * robust_setup.dpVal / max(yData);
    
    
    % experimental data
%     if isfield(ExpData.metabolites{j}, 'conc')
%     if data.ordered.metabolites.onoff(j) == 1
    if canelas_SS.ordered.metabolites.onoff(j) == 1
% % % % %         scatter(exp_time, exp_conc, ...
%         scatter(exp_time, exp_conc/robust_setup.dpVal*robust_setup.dpVal_yaxis, ...
%                     robust_setup.MarkerSize, 'MarkerEdgeColor', 'black', ...
%                     'MarkerFaceColor', 'white', ...
%                     'LineWidth', 1.5)
% %         scatter(exp_time*robust_setup.dpVal_yaxis/0.4, exp_conc/robust_setup.dpVal*robust_setup.dpVal_yaxis, ...
% %                     robust_setup.MarkerSize, 'MarkerEdgeColor', 'black', ...
% %                     'MarkerFaceColor', 'white', ...
% %                     'LineWidth', 1.5)
        scatter(exp_time*robust_setup.dpVal_yaxis/0.4, exp_conc/robust_setup.dpVal*robust_setup.dpVal_yaxis, ...
                    robust_setup.MarkerSize, 'MarkerEdgeColor', 'black', ...
                    'MarkerFaceColor', 'white', ...
                    'LineWidth', 1.5)
    end
    % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
    % add text label
    metName_pre = erase(legenda.metabolites{j}, ", [mM]");
    metName = erase(metName_pre(1:end-7),",");
    % exceptions
%     if j == 14
%         metName = 'GAP';
%     elseif j == 19
%         metName = 'GLYC';
%     end
    text(sp_h.XLim(2)*0.7,sp_h.YLim(2)*0.25,metName,'FontSize',10)

end       
% % % % robust_setup.forcedSSdata = 0;


%% Figure for SS metabolites
% simulation
f2_h = figure(2);
f2_c_h = get(f2_h,'children');
% 
for i = 1:length(robust_setup.maxRange_concentration_range)
    % option based on simulated axes elsewhere
    robust_setup.maxRange_concentration_range(i) = f2_c_h(end+1-i).YLim(2);
end
% 
robust_setup.maxRange_rate_range(1) = 2;
robust_setup.maxRange_rate_range(2) = 2; 
robust_setup.maxRange_rate_range(3) = 2; 
robust_setup.maxRange_rate_range(4) = 2; 
robust_setup.maxRange_rate_range(5) = 2; 
robust_setup.maxRange_rate_range(6) = 0.1; 
robust_setup.maxRange_rate_range(7) = 2; 
robust_setup.maxRange_rate_range(8) = 3; 
robust_setup.maxRange_rate_range(9) = 3;  
% 
robust_setup.maxRange_rate_range(10) = 3; 
robust_setup.maxRange_rate_range(11) = 3; 
robust_setup.maxRange_rate_range(12) = 3; 
robust_setup.maxRange_rate_range(13) = 3; 

% plotting
robust_setup.maxRange_concentration = 0; % for plots of rates, rathe than concentrations
if exist('fig_h_r_ss','var')
    clear fig_h_r_ss
end
fig_h_r_ss = figure(2002);
fig_h_r_ss.Position = [1921 257 1536 747];
% for j = robust_setup.metabolite_id_range
for j = robust_setup.rate_id_range
    if j ~= 1
        toc
    end
    disp(j)
    tic
    % minimum setup each loop
    robust_setup.rate_id = j;
    robust_setup.maxRange_rate = robust_setup.maxRange_rate_range(j);
% % % %     robust_setup.maxRange_concentration_forced = robust_setup.forcedSSdata_vals(j);
    % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
    % %% arrange data
    % heatmap data
    if exist('tData','var')
        clear tData yData
    end
    tData = [];
    yData = [];
% % % %     for i = 1:robust_setup.simulation_robustness_idxs
    for i = robust_setup.simulation_robustness_idxs
%         tData = [tData; robust_setup.interpSpace];
        tData = [tData; [0:71]'];
% % % %         yData = [yData; interp1(robust_setup.simulation_robustness{i}.T_FF01, robust_setup.simulation_robustness{i}.Y_FF01(:,15), robust_setup.interpSpace, 'pchip')];
        yData = [yData; interp1(canelas_SS.mtlD.D, robust_setup.simulation_robustness{i}.ss.V(:,j)', robust_setup.interpSpace, 'pchip')];
%         % 
%         figure, plot(tData,yData)
    end
    tData = [0; robust_setup.maxRange_time; tData];
    yData = [0; robust_setup.maxRange_rate; yData];
    % 
    if j == 1 % GLTs_SS range correction
%         temp_maxVal_y = 1.5;
        temp_maxVal_y = robust_setup.maxRange_rate_range(j);
        indices = find(abs(yData)>temp_maxVal_y);
        yData(indices) = temp_maxVal_y;
    else
    end

    % simulation data
%     sim_time = robust_setup.simulation_reference{1}.gs.T(3002:end) * robust_setup.dpVal / max(tData); %[1 2 10 100];
    sim_time = canelas_SS.mtlD.D * robust_setup.dpVal / max(tData); %[1 2 10 100];
% % % %     sim_conc = robust_setup.dpVal - robust_setup.simulation_reference{1}.Y_FF01(:,15) * robust_setup.dpVal / max(yData);
%     sim_conc = robust_setup.dpVal - robust_setup.simulation_reference{1}.gs.Y(3002:end,j) * robust_setup.dpVal / max(yData);
    temp_val = zeros(8,1);
        temp_val(1) = robust_setup.simulation_reference{1}.ss.V.ss2can{1}(end,j);
        temp_val(2) = robust_setup.simulation_reference{1}.ss.V.ss2can{2}(end,j);
        temp_val(3) = robust_setup.simulation_reference{1}.ss.V.ss2can{3}(end,j);
        temp_val(4) = robust_setup.simulation_reference{1}.ss.V.ss2can{4}(end,j);
        temp_val(5) = robust_setup.simulation_reference{1}.ss.V.ss2can{5}(end,j);
        temp_val(6) = robust_setup.simulation_reference{1}.ss.V.ss2can{6}(end,j);
        temp_val(7) = robust_setup.simulation_reference{1}.ss.V.ss2can{7}(end,j);
        temp_val(8) = robust_setup.simulation_reference{1}.ss.V.ss2can{8}(end,j);
    sim_conc = robust_setup.dpVal - temp_val * robust_setup.dpVal / max(yData);

    % experimental data
%     if isfield(ExpData.metabolites{j}, 'conc')
    if canelas_SS.ordered.fluxes.onoff(j) == 1
%         exp_time = ExpData.metabolites{robust_setup.metabolite_id}.time * robust_setup.dpVal / max(tData);
        exp_time = canelas_SS.mtlD.D * robust_setup.dpVal / max(tData);
%         exp_conc = robust_setup.dpVal - ExpData.metabolites{robust_setup.metabolite_id}.conc * robust_setup.dpVal / max(yData);
        exp_conc = robust_setup.dpVal - canelas_SS.ordered.fluxes.data{robust_setup.rate_id} * robust_setup.dpVal / max(yData);
    end
    % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
    % %% Visualization (in loop)
    % % % % if exist('fig_h_c_gs','var')
    % % % %     clear fig_h_c_gs
    % % % % end
    % % % % fig_h_c_gs = figure(1004);
    figure(fig_h_r_ss);

%     if j == 26
%         disp('Stop here.')
%     end

    % make the heatmap plot
    sp_h2 = subplot(5,9,robust_setup.rate_id);
%     [f_h] = DataDensityPlot(tData, yData, robust_setup.contourLevels, fig_h_c_gs, sp_h);
    [f_h, range90] = DataDensityPlot_1D_Y3M1(tData, yData, robust_setup.contourLevels, fig_h_r_ss, sp_h2, robust_setup);
    hold on
    set(gca, 'XTickLabel', [0 0.4]);

    % simulation data
% % % %     plot(sim_time, sim_conc,'k-','LineWidth', robust_setup.LineWidth)
% % % %     plot(sim_time, sim_conc/(max(sim_conc)-min(sim_conc))*robust_setup.dpVal_yaxis,'k-','LineWidth', robust_setup.LineWidth)
% % % %     plot(sim_time, sim_conc/max(sim_conc)*robust_setup.dpVal_yaxis,'k-','LineWidth', robust_setup.LineWidth)
    plot(sim_time, sim_conc/robust_setup.dpVal*robust_setup.dpVal_yaxis,'k-','LineWidth', robust_setup.LineWidth)

%     % 
%     figure,
%     subplot(131), plot(range90(1,:),'.-')
%     subplot(132), plot(range90(2,:),'.-')
%     subplot(133), plot(range90(3,:),'.-')
    temp_range90min = robust_setup.dpVal - range90(2,:) * robust_setup.dpVal / max(yData);
    temp_range90max = robust_setup.dpVal - range90(3,:) * robust_setup.dpVal / max(yData);
    plot(range90(1,:), temp_range90min/robust_setup.dpVal*robust_setup.dpVal_yaxis,'k--','LineWidth', 1)
    plot(range90(1,:), temp_range90max/robust_setup.dpVal*robust_setup.dpVal_yaxis,'k--','LineWidth', 1)
    
% % % %     sim_conc = robust_setup.dpVal - robust_setup.simulation_reference{1}.Y_FF01(:,j) * robust_setup.dpVal / max(yData);
    
    
    % experimental data
%     if isfield(ExpData.metabolites{j}, 'conc')
%     if data.ordered.metabolites.onoff(j) == 1
    if canelas_SS.ordered.fluxes.onoff(j) == 1
% % % % %         scatter(exp_time, exp_conc, ...
%         scatter(exp_time, exp_conc/robust_setup.dpVal*robust_setup.dpVal_yaxis, ...
%                     robust_setup.MarkerSize, 'MarkerEdgeColor', 'black', ...
%                     'MarkerFaceColor', 'white', ...
%                     'LineWidth', 1.5)
        scatter(exp_time*robust_setup.dpVal_yaxis/0.4, exp_conc/robust_setup.dpVal*robust_setup.dpVal_yaxis, ...
                    robust_setup.MarkerSize, 'MarkerEdgeColor', 'black', ...
                    'MarkerFaceColor', 'white', ...
                    'LineWidth', 1.5)
    end
    % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
    % add text label
    rateName_pre = erase(legenda.fluxes{j}, ", [mM s^{-1}]");
    rateName = erase(rateName_pre(1:end-7),",");
    text(sp_h.XLim(2)*0.1,sp_h.YLim(2)*0.25,rateName,...
        'FontSize',10,'HorizontalAlignment','Left')
    % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 

end       
% % % % robust_setup.forcedSSdata = 0;


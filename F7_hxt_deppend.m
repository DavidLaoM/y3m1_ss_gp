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
% % % % % % 
% % % % % x_2nd = x_3rd;


%% Starting simulation
% simulation
setup.plotResultsMode = 0;
[simRes_ref] = simulateY3M1_separateGPSS(x_2nd, canelas_SS, data, dataset, setup);
setup.simGPdata = simRes_ref{1}.gs;
setup.simSSdata = simRes_ref{1}.ss;
setup.litParams = zeros(size(x));
setup.plotResultsMode = 0;


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


%% plot data vs dilRate
% close(5001)
dilRate = canelas_SS.mtlD.D;
Cs = canelas_SS.mtlD.Cs;
fh_5001 = figure(5001); 
% fh_5001.Position = [2616 527 480 420];
% fh_5001.Position = [2616 527 588 420]; 
fh_5001.Position = [1657 185 247 799]; 
% 
% maxVal = 1.1 * 1.5; 
% line([0 maxVal],[0 maxVal],'Color','r','LineStyle','--')
subplot(3,1,1)
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
xticks([0 0.1 0.2 0.3 0.4]), xlim([0 0.4])%, xlim([0 maxVal_x])
yticks([0 .5 1.0 1.5]), ylim([0 1.5])%, ylim([0 maxVal_y])
% %% 
xlabel(Y3M1_labels.dilRate)
ylabel(Y3M1_labels.qT_mMs)
% 
% txt1 = text(0.0625,1.515,'A');
% txt1.FontSize = 15;

% %% parameter bar plot
safecopy_x36 = 1.0445;
bplot_vals = [-0.05-safecopy_x36 0.8-safecopy_x36 0 0,...
    1.15-safecopy_x36 1.5-safecopy_x36 1.26-safecopy_x36 1.2-safecopy_x36];
bplot_c = categorical({'0.02','0.05','0.10','0.20',...
                '0.30','0.325','0.35','0.375'});

% ax3 = axes('Position',[.25 .675 .24 .2]);
% ax3 = axes('Position',[.2 .5 .24 .275]);
% % % % ax3 = axes('Position',[.64 .225 .24 .275]);
subplot(3,1,2)
box on
cp3 = bar(bplot_c,bplot_vals,'FaceColor','k');
xtickangle(75)
% ax3.YLim = [-1.25 1.25];
% ax3.YTick = [-1 0 1];
% 
% txt2 = text(1,0.9,'C');
% txt2.FontSize = 15;
ylabel('Fold Change')
xlabel(Y3M1_labels.dilRate)
ylim([-1.5 1.5])

% ax4 = axes('Position',[.625 .175 .24 .25]);
% % % % ax4 = axes('Position',[.19 .46 .24 .29]);
subplot(3,1,3)
box on
plot(dilRate,Cs,'r.--')
% ax4.XTick = [0 2];
% ax4.YTick = [0 2];
% txt4 = text(ax4.XLim(2)*0.1, ax4.YLim(2)*0.8, 'B');
% txt4.FontSize = 15;
ylabel('Residual glucose (mM)')
xlabel(Y3M1_labels.dilRate)


% % saving the plots
% savefig(5001,'results\F7_glt_mu_dependency.fig')
% savefig(1041,'results\supF_glt_mu_dependency.fig')








%SINKS DEVELOPMENT
% 0- Points where inequilibrium exists
% 1- Obtaining the polynomials
% 2- Monod version
% 3- Hill version
clc, clear, close all
set_paths;
CanelasData2;


%% (0a) Points where inequilibrium exists
massBalance.Glci    = canelas_SS.mtlD.v_HK - canelas_SS.mtlD.v_GLT; %Glci, GLT_HXK
massBalance.G6P     = canelas_SS.mtlD.v_PGI - canelas_SS.mtlD.v_HK; %G6P, PGI_HXK
massBalance.F6P     = canelas_SS.mtlD.v_PFK - canelas_SS.mtlD.v_PGI; %F6P, PFK_PGI
massBalance.FBP     = canelas_SS.mtlD.v_FBA - canelas_SS.mtlD.v_PFK; %FBP, FBA_PFK
massBalance.DHAP    = canelas_SS.mtlD.v_G3PDH + canelas_SS.mtlD.v_TPI - canelas_SS.mtlD.v_FBA; %DHAP, G3PDH_TPI_FBA
massBalance.GAP     = canelas_SS.mtlD.v_GAPDH + canelas_SS.mtlD.v_G3PDH - 2*canelas_SS.mtlD.v_FBA; %UG_LG
massBalance.BPG     = canelas_SS.mtlD.v_PGK - canelas_SS.mtlD.v_GAPDH; %BPG, PGK_GAPDH
massBalance.P3G     = canelas_SS.mtlD.v_PGM - canelas_SS.mtlD.v_PGK; %P3G, PGM_PGK
massBalance.P2G     = canelas_SS.mtlD.v_ENO - canelas_SS.mtlD.v_PGM; %P2G, ENO_PGM
massBalance.PEP     = canelas_SS.mtlD.v_PYK - canelas_SS.mtlD.v_ENO; %PEP,PYK_ENO
massBalance.PYR     = canelas_SS.mtlD.v_PDC - canelas_SS.mtlD.v_PYK; %PYR, PDC_PYK
massBalance.ADH     = canelas_SS.mtlD.v_ADH - canelas_SS.mtlD.v_PDC; %ADH, ADH_PDC

% Q = struct2cell(massBalance);
% titles = {'Glci','G6P','F6P','FBP','DHAP','LG_UG','BPG','P3G','P2G','PEP','PYR','ADH'};
% figure(1)
% for i = 1:length(Q)
%     subplot(3,4,i)
% %     plot(canelas_SS.mtlD.D, Q{i}, 'k-.', 'MarkerSize', 2)
%     area([0 0.4], [0.001 0.001],'FaceColor',[1.0 0.8 0.4],'EdgeColor','none');
%     hold on
%     area([0 0.4], [-0.001 -0.001],'FaceColor',[1.0 0.8 0.4],'EdgeColor','none');
%     hold on
%     plot(canelas_SS.mtlD.D, Q{i}, 'k-o', 'MarkerSize', 2)
%     title(titles{i})
%     hold off
%     if i == 9
%         xlabel('Growth rate, /mu, [h^{-1}]')
%         ylabel('Reaction rate, v, [mM s^{-1}]')
%     end
% %     xlabel('Growth rate, mu, [h-1]')
% %     ylabel('Reaction rate, v, [mM s-1]')
% end
% suptitle('mass Balances for metabolites in the canelas dataset')


%% (0b) importance of the sinks
canelas = canelas_SS;
balance_legend  = {'GLCi';    'G6P';     'F6P';     'FBP';     'DHAP';     'GAP';     'BPG';     'P3G';     'P2G';     'PEP';     'PYR';     'ACE'};
balance_in      = [canelas.mtlD.v_GLT', canelas.mtlD.v_HK', canelas.mtlD.v_PGI', canelas.mtlD.v_PFK',     canelas.mtlD.v_FBA',     (canelas.mtlD.v_FBA + canelas.mtlD.v_TPI)',     canelas.mtlD.v_GAPDH',     canelas.mtlD.v_PGK',     canelas.mtlD.v_PGM',     canelas.mtlD.v_ENO',     canelas.mtlD.v_PYK',     canelas.mtlD.v_PDC'];
balance_out     = [canelas.mtlD.v_HK', canelas.mtlD.v_PGI', canelas.mtlD.v_PFK', canelas.mtlD.v_FBA',     (canelas.mtlD.v_TPI + canelas.mtlD.v_G3PDH)',     canelas.mtlD.v_GAPDH',     canelas.mtlD.v_PGK',     canelas.mtlD.v_PGM',     canelas.mtlD.v_ENO',   canelas.mtlD.v_PYK',     canelas.mtlD.v_PDC',     canelas.mtlD.v_ADH'];

% figure(2)
% for i = 1:12
%     subplot(3,4,i)
%     plot(canelas.mtlD.D, balance_in(:,i), 'o-', canelas.mtlD.D, balance_out(:,i), '*-')
%     title(balance_legend(i))
% end
% legend('produced', 'consumed')
% suptitle('Fluxes producing and consuming metabolite. Sinks not considered')
% ylabel('Reaction rate [mM h^{-1}]')
% xlabel('Growth rate [h^{-1}]')
% xlabh = get(gca,'XLabel');
% set(xlabh,'Position',get(xlabh,'Position') - [.8 .3 0])
% ylabh = get(gca,'YLabel');
% set(ylabh,'Position',get(ylabh,'Position') - [1.65 -3.5 0])


%% (1) Polynomial fit
d = canelas_SS.mtlD.D;
% polynomial fits
dataFit{1} = polyfit(d,massBalance.G6P,3);
dataFit{2} = polyfit(d,massBalance.F6P,6);
dataFit{3} = polyfit(d,massBalance.GAP,6);
dataFit{4} = polyfit(d,massBalance.P3G,2);
dataFit{5} = polyfit(d,massBalance.PEP,2);
dataFit{6} = polyfit(d,massBalance.PYR,6);
dataFit{7} = polyfit(d,massBalance.ADH,6);
% recalling experimental data
expData{1} = massBalance.G6P;
expData{2} = massBalance.F6P;
expData{3} = massBalance.GAP;
expData{4} = massBalance.P3G;
expData{5} = massBalance.PEP;
expData{6} = massBalance.PYR;
expData{7} = massBalance.ADH;

titles2 = {'G6P','F6P','GAP','P3G','PEP','PYR','ADH'};
%%
if exist('fh10','var')
    clf(fh10)
end
fh10 = figure(10);
fh10.Position = [1900 545 982 435];
for i = 1:length(dataFit)
    sp_h = subplot(2,4,i);
    plot(d, expData{i}, 'ko', 'MarkerSize', 4, 'MarkerFaceColor', 'k')
    hold on
    plot(d, polyval(dataFit{i}, d), 'b-', 'MarkerSize', 2, 'LineWidth', 1.5)
    %
%     xlim([0 0.5])
    ylim([sp_h.YLim(1) sp_h.YLim(2)*1.25])
    temp_X = sp_h.XLim(2) * 0.8;
    temp_Y = (sp_h.YLim(2) - sp_h.YLim(1)) * 0.8 + sp_h.YLim(1);
    text(temp_X,temp_Y,titles2{i},'HorizontalAlignment','Center')
    ylabel('Reaction Rate (mM h^{-1})')
    xlabel('Dilution Rate (h^{-1})')
end
Qleg = legend('Calculated','Fitted');
set(Qleg, 'Position', [0.8,0.125,0.075,0.10]);
% suptitle('Comparison difference between reactions vs growth rate')

%%
h10 = figure(10);
set(h10,'Units','Inches');
pos = get(h10,'Position');
set(h10,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(h10,'E:\OneDrive - TU Eindhoven\Documents\ch4_Y3M1\results\figures\manuscript\SM_sinks','-dpdf','-r0')

%%
if exist('fh11','var')
    clf(fh11)
end
fh11 = figure(11);
fh11.Position = [2093 692 305 256];
% for i = 1:length(dataFit)
%     sp_h = subplot(2,4,i);
    i = 6; 
    sp_h = plot(d, expData{i}, 'ko', 'MarkerSize', 4, 'MarkerFaceColor', 'k');
    hold on
    plot(d, polyval(dataFit{i}, d), 'b-', 'MarkerSize', 2, 'LineWidth', 1.5)
    %
%     xlim([0 0.5])
%     ylim([sp_h.YLim(1) sp_h.YLim(2)*1.25])
%     temp_X = sp_h.XLim(2) * 0.8;
%     temp_Y = (sp_h.YLim(2) - sp_h.YLim(1)) * 0.8 + sp_h.YLim(1);
    text(0.35,-0.05,titles2{i},'HorizontalAlignment','Center',...
        'FontSize',10)
    ylabel('Reaction Rate (mM h^{-1})')
    xlabel('Dilution Rate (h^{-1})')
% end
Qleg = legend('Calculated','Fitted');
% set(Qleg, 'Position', [0.8,0.125,0.075,0.10]);
set(Qleg, 'Location', 'SouthWest');
% suptitle('Comparison difference between reactions vs growth rate')
%%
% % % % print(gcf,'-dpdf','E:\OneDrive - TU Eindhoven\Documents\ch4_Y3M1\results\figures\manuscript\SM_sinkPYR_preAI');





% % saving the plots
% savefig(10,'results\manuscript_supplementary_materials\supF_sinks.fig')





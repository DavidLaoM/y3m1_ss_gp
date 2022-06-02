% figure parameter comparison
set_paths
clear, close all

% (1) Recall a code with plot in vitro vs in vivo
% (2) Rerun the part needed from it (using pE99b? + change in km) and split
% in parts. First just in suplots
% (3) Then adjust the size of the suplots
    % line 1:   2 + 7 + 4 + 14 + 5 = 32
        % 35 to 36          GLT
        % 28 to 34          HXK
        % 57 to 60          PGI
        % 43 to 56 + 87     PFK
        % 11 to 15          FBA
    % line 2:   5 + 4 + 9 + 7 = 26
        % 79 to 82          TPI
        % 96 to 104         GPD1 (G3PDH)
        % 21 to 27          TDH (GAPDH)
        % 61 to 66          PGK
    % line 3:   4 + 4 + 6 + 4 + 10 = 28
        % 67 to 70          GPM1 (PGM)
        % 16 to 19          ENO
        % 71 to 77          PYK
        % 39 to 42          PDC
        % 1 to 10           ADH
        
% % % % %% 
% % % % load('parSet_99b.mat') % updated estimation Y3M1
% % % % load('parSet_104a.mat') % updated estimation Y3M1
% % % % load('parSet_105b.mat') % updated estimation Y3M1 (PFK and PYK regularized)
% % % % load('x16d.mat','x16d_final'); x16d_final(35) = x16d_final(38); % estimation Y3M2
% % % % find(x99b - xFinal ~= 0);
% % % % %%
% % % % x99b([35 38]);
% % % % xFinal([35 38]);
% % % % %%
% % % % x = xFinal;
% % % % % plot_InVitro_vs_InVivo
% % % % x(43:56) = x105b(43:56);
% % % % x(71:77) = x105b(71:77);

%%
load('x108c.mat'),
x = x_3rd;


%% neat version plot in vitro vs in vivo parameters

% pars numbers
parsHXT = 35:36; %GLT
parsHXK = 28:34;  %HXK
parsPGI = 57:60;  %PGI
parsPFK = 43:56;  %PFK
parsALD = 11:15;  %ALD
% 
parsGPD = 96:104;  %GPD
parsTPI = 79:82;  %TPI
parsTDH = 21:27;  %GAPDH
parsPGK = 61:66;  %PGK
% 
parsPGM = 67:70;  %PGM
parsENO = 16:19;  %ENO
parsPYK = 71:77;  %PYK
parsPDC = 39:42;  %PDC
parsADH = 1:10;  %ADH

% tick labels
x_ticks = cell(1,14);
x_ticks{1} = parsHXT; % = 35:36; %GLT
x_ticks{2} = parsHXK; % = 28:34;  %HXK
x_ticks{3} = parsPGI; % = 57:60;  %PGI
x_ticks{4} = parsPFK; % = 43:56;  %PFK
x_ticks{5} = parsALD; % = 11:15;  %ALD
% 
x_ticks{6} = parsGPD; % = 96:104;  %GPD
x_ticks{7} = parsTPI; % = 79:82;  %TPI
x_ticks{8} = parsTDH; % = 21:27;  %GAPDH
x_ticks{9} = parsPGK; % = 61:66;  %PGK
% 
x_ticks{10} = parsPGM; % = 67:70;  %PGM
x_ticks{11} = parsENO; % = 16:19;  %ENO
x_ticks{12} = parsPYK; % = 71:77;  %PYK
x_ticks{13} = parsPDC; % = 39:42;  %PDC
x_ticks{14} = parsADH; % = 1:10;  %ADH

% enzyme names
enz_names = cell(1,14);
enz_names{1} = 'HXT'; % = 35:36; %GLT
enz_names{2} = 'HXK'; % = 28:34;  %HXK
enz_names{3} = 'PGI'; % = 57:60;  %PGI
enz_names{4} = 'PFK'; % = 43:56;  %PFK
enz_names{5} = 'ALD'; % = 11:15;  %ALD
% 
enz_names{6} = 'GPD'; % = 96:104;  %GPD
enz_names{7} = 'TPI'; % = 79:82;  %TPI
enz_names{8} = 'TDH'; % = 21:27;  %GAPDH
enz_names{9} = 'PGK'; % = 61:66;  %PGK
% 
enz_names{10} = 'PGM'; % = 67:70;  %PGM
enz_names{11} = 'ENO'; % = 16:19;  %ENO
enz_names{12} = 'PYK'; % = 71:77;  %PYK
enz_names{13} = 'PDC'; % = 39:42;  %PDC
enz_names{14} = 'ADH'; % = 1:10;  %ADH

% pretreatment legend names
[legenda] = legendaFull;
legenda.parameters2 = legenda.parameters;
for i = 1:length(legenda.parameters2)
    % take out the beginning
    if i < 10
        legenda.parameters2{i} = legenda.parameters2{i}(6:end);
    elseif i < 100
        if((i >= 28)&&(i <= 36))
            legenda.parameters2{i} = legenda.parameters2{i}(10:end);
        elseif((i >= 43)&&(i <= 56))
            legenda.parameters2{i} = legenda.parameters2{i}(10:end);
        elseif(i == 87)
            legenda.parameters2{i} = legenda.parameters2{i}(10:end);
        else
            legenda.parameters2{i} = legenda.parameters2{i}(7:end);
        end
    else
        legenda.parameters2{i} = legenda.parameters2{i}(8:end);
    end
    % erase the parts in the    end
    if isempty(legenda.parameters2{i}) ~= 1
        str1 = convertCharsToStrings(legenda.parameters2{i}); legenda.parameters2{i} = erase(str1,", mM s^{-1}");
        str2 = convertCharsToStrings(legenda.parameters2{i}); legenda.parameters2{i} = erase(str2,", mM");
        str3 = convertCharsToStrings(legenda.parameters2{i}); legenda.parameters2{i} = erase(str3,", s^{-1}");
        str4 = convertCharsToStrings(legenda.parameters2{i}); legenda.parameters2{i} = erase(str4,", []");
        str5 = convertCharsToStrings(legenda.parameters2{i}); legenda.parameters2{i} = erase(str5,",");
    end
end
% % %%
% for i = 1:120
%     disp(legenda.parameters2{i})
% end
% %%


%% plots all pars. Dotplot.
% figures
load('customColormap.mat','customColormap')
% 
if exist('f_h','var')
    clf(f_h)
end
% f_h = figure(20001);
f_h = figure('Name','All parameters dotplot');
% %
% sp_h(1) = subplot(3,5,1);
% for i = parsHXT
% %     temp
%     hold on
%     color=customColormap(int8(1+127*(max(abs(x(i)))/3)),:);
%     sl_h = semilogy(i,10.^x_reorder(i),'o','MarkerFaceColor',color,'MarkerEdgeColor',color);
%     line([],[1 1],'Color','k')
% end
for j = 1:length(x_ticks)
    % plotting the data of interest
    if j < 4
        sp_h(j) = subplot(3,5,j);
    elseif j == 4
        sp_h(j) = subplot(3,5,[j j+1]);
    elseif j > 4
        sp_h(j) = subplot(3,5,j+1);
    end
    % 
    for i = x_ticks{j}
    %     temp
        color=customColormap(int8(1+127*(max(abs(x(i)))/3)),:);
        sl_h = semilogy(i,10.^x(i),'o','MarkerFaceColor',color,'MarkerEdgeColor',color);
        hold on
        if i == x_ticks{j}(1)
%             area([1 min(x_ticks{j})-1; 1 max(x_ticks{j})+1], 'LineStyle', 'none', 'FaceColor ', [.9 .9 .9])
%             area([min(x_ticks{j})-1 max(x_ticks{j})+1], [-1 -1], 'LineStyle', 'none', 'FaceColor ', [.9 .9 .9])
            line([min(x_ticks{j})-1 max(x_ticks{j})+1], [1 1],'Color','k','Linestyle','-')
            line([min(x_ticks{j})-1 max(x_ticks{j})+1], [10 10],'Color','k','Linestyle',':')
            line([min(x_ticks{j})-1 max(x_ticks{j})+1], [.1 .1],'Color','k','Linestyle',':')
        end
        sl_h = semilogy(i,10.^x(i),'o','MarkerFaceColor',color,'MarkerEdgeColor',color);

% % % %         % adding parameters that change upon GP (empty black circles)
% % % %         if((j == 5)&&(i == 11)) % FBA 
% % % % %             semilogy(i,10.^x(147),'o','MarkerFaceColor','black','MarkerEdgeColor','black');
% % % %         elseif((j == 4)&&(i == 52)) % PFK 
% % % %             semilogy(i,10.^x(148),'o','MarkerFaceColor','black','MarkerEdgeColor','black');
% % % %         elseif((j == 3)&&(i == 59)) % PGI 
% % % %             semilogy(i,10.^x(149),'o','MarkerFaceColor','black','MarkerEdgeColor','black');
% % % %         elseif((j == 7)&&(i == 79)) % TPI 
% % % %             semilogy(i,10.^x(150),'o','MarkerFaceColor','black','MarkerEdgeColor','black');
% % % %         elseif((j == 8)&&(i == 22)) % GAPDH 
% % % %             semilogy(i,10.^x(151),'o','MarkerFaceColor','black','MarkerEdgeColor','black');
% % % %         end

% % % %         % adding parameters that change upon FF (filled black circles)
% % % %         % x16d_final
% % % %         if(j == 1) % GLT
% % % % %             semilogy(35,10.^x16d_final(35),'o','MarkerFaceColor','blue','MarkerEdgeColor','blue');
% % % %             semilogy(36,10.^x16d_final(36),'o','MarkerFaceColor','blue','MarkerEdgeColor','blue');
% % % %         elseif(j == 2) % HXK 
% % % % %             semilogy(28,10.^x16d_final(28),'o','MarkerFaceColor','blue','MarkerEdgeColor','blue');
% % % % %             semilogy(29,10.^x16d_final(29),'o','MarkerFaceColor','blue','MarkerEdgeColor','blue');
% % % % %             semilogy(30,10.^x16d_final(30),'o','MarkerFaceColor','blue','MarkerEdgeColor','blue');
% % % % %             semilogy(31,10.^x16d_final(31),'o','MarkerFaceColor','blue','MarkerEdgeColor','blue');
% % % %             semilogy(32,10.^x16d_final(32),'o','MarkerFaceColor','blue','MarkerEdgeColor','blue');
% % % %             semilogy(33,10.^x16d_final(33),'o','MarkerFaceColor','blue','MarkerEdgeColor','blue');
% % % %             semilogy(34,10.^x16d_final(34),'o','MarkerFaceColor','blue','MarkerEdgeColor','blue');
% % % %         end
    end
    
    % commong editing
    sp_h(j).XLim = [min(x_ticks{j})-1 max(x_ticks{j})+1]; % ylims
    sp_h(j).YLim = [10^-3 10^3]; % xlims
    sp_h(j).XTick = x_ticks{j}; % xtickaxes
%     sp_h(j).YTick = [1E-3 1E-2 1E-1 1E0 1E1 1E2 1E3]; % ytickaxes
%     sp_h(j).YTickLabel = {'10^{-3}','10^{-2}','10^{-1}','10^{0}','10^{1}','10^{2}','10^{3}'};
    sp_h(j).YTick = [1E-3 1E-1 1E0 1E1 1E3]; % ytickaxes
    sp_h(j).YTickLabel = {'10^{-3}','10^{-1}','10^{0}','10^{1}','10^{3}'};
    sp_h(j).FontSize = 12;%12;
    text_h = text(sp_h(j).XLim(1) + 0.25, 2E2, enz_names{j},...
        'FontSize',12); % textbox / title
    sp_h(j).XTickLabel = legenda.parameters2(x_ticks{j});
    xtickangle(45)
        
    % colorbar spawning after the last plot
    if j == length(x_ticks)
        N = 64;   % or whatever
        lin_fwd = linspace(0,1,64)';
        lin_rev = flip(lin_fwd);
        colmin_green = [0    0.4980         0];
        colmax_red = [1.0000    0.0100         0];
        Color1 = [colmin_green;
            zeros(N,3);
            colmax_red];
        Color1(2:65,1) = lin_fwd;
        Color1(2:65,2) = colmin_green(2)*lin_rev + colmax_red(2)*lin_fwd;
%         Color1 = [ones(N,1),(N-1:-1:0)'/(N-1),zeros(N,1)];
        colormap(Color1)
        cb_h = colorbar;
        cb_h.Ticks = [0 1/3 2/3 1];
        cb_h.TickLabels = {'0', '1', '2', '3'};
%         cb_h.Limits = [0 3];
        cb_h.Position = [0.02 + sp_h(9).Position(1) + sp_h(9).Position(3),...
            sp_h(9).Position(2),...
            cb_h.Position(3),...
            cb_h.Position(4)];
    end
end

% % resizing individual plots
% 
sp_h(5).Position = sp_h(5).Position + [0 0 -0.03 0]; % bit left
sp_h(6).Position = sp_h(6).Position + [-0.03 0 0.07 0]; % bit expansion both 
sp_h(7).Position = sp_h(7).Position + [0.04 0 -0.03 0]; % bit right
sp_h(8).Position = sp_h(8).Position + [0.01 0 0 0]; % bit right
sp_h(9).Position = sp_h(9).Position + [0.01 0 -0.01 0]; % bit right
% 
sp_h(10).Position = sp_h(10).Position + [0 0 -0.01 0]; % bit left
sp_h(11).Position = sp_h(11).Position + [-0.01 0 -0.01 0]; % bit right
sp_h(12).Position = sp_h(12).Position + [-0.02 0 0 0]; % bit right
sp_h(13).Position = sp_h(13).Position + [-0.02 0 -0.01 0]; % bit right
sp_h(14).Position = sp_h(14).Position + [-0.03 0 0.07 0]; % bit right
    sp_h(14).Position(3) = 0.1537;
% % 
% sp_h(4).Position([2 4]) = sp_h(3).Position([2 4]); % bit left

% resizing full plot
% set(f_h,'color','white')
hold off

% %%
% for i = 1:14
%     disp(sp_h(i).Position)
% end

% %% resize
% no idea why it has to be here
for i = 1:3
    sp_h(i).Position([2 4]) = sp_h(4).Position([2 4]); 
end
%
set(f_h, 'Position', [100 100 1920*0.9 1083*0.9])

%% printing
% set(f_h,'Units','inches');
% screenposition = get(f_h,'Position');
% set(f_h,...
%     'PaperPosition',[0 0 screenposition(3:4)],...
%     'PaperSize',[screenposition(3:4)]);
% print -dpdf -painters figure_9_pSet
% print -dpng -painters figure_9_pSet


%% plotting in dot style, but by category
% 
x_ticks_keq = cell(1,14);
    x_ticks_keq{1} = []; % = 35:36; %GLT
    x_ticks_keq{2} = 30; % = 28:34;  %HXK
    x_ticks_keq{3} = 57; % = 57:60;  %PGI
    x_ticks_keq{4} = []; % = 43:56;  %PFK
    x_ticks_keq{5} = 12; % = 11:15;  %ALD
    x_ticks_keq{6} = 99; % = 96:104;  %GPD
    x_ticks_keq{7} = 80; % = 79:82;  %TPI
    x_ticks_keq{8} = 21; % = 21:27;  %GAPDH
    x_ticks_keq{9} = 61; % = 61:66;  %PGK
    x_ticks_keq{10} = 69; % = 67:70;  %PGM
    x_ticks_keq{11} = 17; % = 16:19;  %ENO
    x_ticks_keq{12} = []; % = 71:77;  %PYK
    x_ticks_keq{13} = []; % = 39:42;  %PDC
    x_ticks_keq{14} = 1; % = 1:10;  %ADH
x_ticks_keq_all = [];
x_ticks_keq_names = cell(1,14);
for i = 1:length(x_ticks_keq)
    x_ticks_keq_all = [x_ticks_keq_all, x_ticks_keq{i}];
    x_ticks_keq_names{i} = enz_names{i};
end
x_ticks_keq_names(13) = [];
x_ticks_keq_names(12) = [];
x_ticks_keq_names(4) = [];
x_ticks_keq_names(1) = [];
% 
x_ticks_km = cell(1,14);
    x_ticks_km{1} = 35; % = 35:36; %GLT
    x_ticks_km{2} = [28 29 31 32]; % = 28:34;  %HXK
    x_ticks_km{3} = [58 59]; % = 57:60;  %PGI
    x_ticks_km{4} = [48 49 50 52]; % = 43:56;  %PFK
    x_ticks_km{5} = [11 13 14]; % = 11:15;  %ALD
    x_ticks_km{6} = [96 97 98 100 101 102 103]; % = 96:104;  %GPD
    x_ticks_km{7} = [79 81]; % = 79:82;  %TPI
    x_ticks_km{8} = [22 23 24 25 26]; % = 21:27;  %GAPDH
    x_ticks_km{9} = [62 63 64 65]; % = 61:66;  %PGK
    x_ticks_km{10} = [67 68]; % = 67:70;  %PGM
    x_ticks_km{11} = [16 18]; % = 16:19;  %ENO
    x_ticks_km{12} = [71 72 73 74]; % = 71:77;  %PYK
    x_ticks_km{13} = [39 40]; % = 39:42;  %PDC
    x_ticks_km{14} = [6 7 8 9]; % = 1:10;  %ADH
x_ticks_km_all = [];
% x_ticks_km_names = cell(1,14);
for i = 1:length(x_ticks_km)
    x_ticks_km_all = [x_ticks_km_all, x_ticks_km{i}];
%     x_ticks_km_names{i} = num2str(x_ticks_km{i});
end
% 
x_ticks_kcat_vm = cell(1,14);
    x_ticks_kcat_vm{1} = 36; % = 35:36; %GLT
    x_ticks_kcat_vm{2} = 34; % = 28:34;  %HXK
    x_ticks_kcat_vm{3} = 60; % = 57:60;  %PGI
    x_ticks_kcat_vm{4} = 56; % = 43:56;  %PFK
    x_ticks_kcat_vm{5} = 15; % = 11:15;  %ALD
    x_ticks_kcat_vm{6} = 104; % = 96:104;  %GPD
    x_ticks_kcat_vm{7} = 82; % = 79:82;  %TPI
    x_ticks_kcat_vm{8} = 27; % = 21:27;  %GAPDH
    x_ticks_kcat_vm{9} = 66; % = 61:66;  %PGK
    x_ticks_kcat_vm{10} = 70; % = 67:70;  %PGM
    x_ticks_kcat_vm{11} = 19; % = 16:19;  %ENO
    x_ticks_kcat_vm{12} = 77; % = 71:77;  %PYK
    x_ticks_kcat_vm{13} = 42; % = 39:42;  %PDC
    x_ticks_kcat_vm{14} = 10; % = 1:10;  %ADH
x_ticks_kcat_vm_all = [];
x_ticks_kcat_vm_names = cell(1,14);
for i = 1:length(x_ticks_kcat_vm)
    x_ticks_kcat_vm_all = [x_ticks_kcat_vm_all, x_ticks_kcat_vm{i}];
    x_ticks_kcat_vm_names{i} = enz_names{i};
end

% extra labels km
temp_leg = legenda.parameters(x_ticks_km_all);
legenda.parameters2 = cell(1,length(temp_leg));
for i = 1:length(temp_leg)
    temp1 = erase(temp_leg{i},", mM");
    if((i < 12)&&(i ~= 7)&&(i ~= 6))
        legenda.parameters2{i} = temp1(10:end);
    elseif i < 43
        legenda.parameters2{i} = temp1(7:end);
    else
        legenda.parameters2{i} = temp1(6:end);
    end
end


%% plots Keq-specific. Dotplot.
if exist('f_h2','var')
    clf(f_h2)
end
% f_h2 = figure(20002);
f_h2 = figure('Name','Keq dotplot');

% plotting Keq
sp_h(1) = subplot(1,1,1);
for i = 1:length(x_ticks_keq_all)
        color=customColormap(int8(1+127*(max(abs(x(x_ticks_keq_all(i))))/3)),:);
        semilogy(i,10.^x(x_ticks_keq_all(i)),'o','MarkerFaceColor',color,'MarkerEdgeColor',color);
        hold on
        if i == 1
            line([0 11], [1 1],'Color','k','Linestyle','-')
            line([0 11], [10 10],'Color','k','Linestyle',':')
            line([0 11], [.1 .1],'Color','k','Linestyle',':')
        end
        sl_h = semilogy(i,10.^x(x_ticks_keq_all(i)),'o','MarkerFaceColor',color,'MarkerEdgeColor',color);
end
% commong editing
sp_h(1).XLim = [0 11]; % ylims
sp_h(1).YLim = [10^-3 10^3]; % xlims
sp_h(1).XTick = 1:length(x_ticks_keq_all); % xtickaxes 
sp_h(1).YTick = [1E-3 1E-1 1E0 1E1 1E3]; % ytickaxes
sp_h(1).YTickLabel = {'10^{-3}','10^{-1}','10^{0}','10^{1}','10^{3}'};
sp_h(1).FontSize = 12;%12;
% sp_h(1).XTickLabel = legenda.parameters(x_ticks_keq_all);%x_ticks_keq_names;
sp_h(1).XTickLabel = x_ticks_keq_names;
xtickangle(45)


%% plots Km-specific. Dotplot.
if exist('f_h3','var')
    clf(f_h3)
end
% f_h3 = figure(20003);
f_h3 = figure('Name','Km dotplot');

% plotting Km
sp_h(1) = subplot(1,1,1);
for i = 1:length(x_ticks_km_all)
        color=customColormap(int8(1+127*(max(abs(x(x_ticks_km_all(i))))/3)),:);
        semilogy(i,10.^x(x_ticks_km_all(i)),'o','MarkerFaceColor',color,'MarkerEdgeColor',color);
        hold on
        if i == 1
            % horizontal lines (:)
            line([0 47], [1 1],'Color','k','Linestyle','-')
            line([0 47], [10 10],'Color','k','Linestyle',':')
            line([0 47], [.1 .1],'Color','k','Linestyle',':')
            % vertical lines (--)
            tempLen_2 = 0;
            for j = 1:length(x_ticks_km)
                tempLen_2 = tempLen_2 + length(x_ticks_km{j});
                % add line
                line(zeros(1,2) + 0.5 + tempLen_2, [1E-3 1E3], 'Color','k','Linestyle',':') % GLT
                % add enzyme name
                text_h = text(0 + 0.5 + tempLen_2 - length(x_ticks_km{j}) / 2, 2000, enz_names{j},...
                    'FontSize',12,'HorizontalAlignment','center'); % textbox / title
            end
        end
        sl_h = semilogy(i,10.^x(x_ticks_km_all(i)),'o','MarkerFaceColor',color,'MarkerEdgeColor',color);

end
% commong editing
sp_h(1).XLim = [0 47]; % ylims
sp_h(1).YLim = [10^-3 10^3]; % xlims
sp_h(1).XTick = 1:length(x_ticks_km_all); % xtickaxes 
sp_h(1).YTick = [1E-3 1E-1 1E0 1E1 1E3]; % ytickaxes
sp_h(1).YTickLabel = {'10^{-3}','10^{-1}','10^{0}','10^{1}','10^{3}'};
sp_h(1).FontSize = 12;%12;
sp_h(1).XTickLabel = legenda.parameters2;
xtickangle(45)
% %%
% set(f_h3, 'Position', [100 100 1920*0.9 1083*0.8])
set(f_h3, 'Position', [100 100 1920*0.9 460])


%% plots Kcat,Vm-specific. Dotplot.
if exist('f_h4','var')
    clf(f_h4)
end
% f_h4 = figure(20004);
f_h4 = figure('Name','Kcat, Vm, dotplot');

% plotting Keq
sp_h4(1) = subplot(1,1,1);
for i = 1:length(x_ticks_kcat_vm_all)
        color=customColormap(int8(1+127*(max(abs(x(x_ticks_kcat_vm_all(i))))/3)),:);
        semilogy(i,10.^x(x_ticks_kcat_vm_all(i)),'o','MarkerFaceColor',color,'MarkerEdgeColor',color);
        hold on
        if i == 1
            line([0 length(x_ticks_kcat_vm_all)+1], [1 1],'Color','k','Linestyle','-')
            line([0 length(x_ticks_kcat_vm_all)+1], [10 10],'Color','k','Linestyle',':')
            line([0 length(x_ticks_kcat_vm_all)+1], [.1 .1],'Color','k','Linestyle',':')
        end
        sl_h = semilogy(i,10.^x(x_ticks_kcat_vm_all(i)),'o','MarkerFaceColor',color,'MarkerEdgeColor',color);
end
% commong editing
sp_h4(1).XLim = [0 length(x_ticks_kcat_vm_all)+1]; % ylims
sp_h4(1).YLim = [10^-3 10^3]; % xlims
sp_h4(1).XTick = 1:length(x_ticks_kcat_vm_all); % xtickaxes 
sp_h4(1).YTick = [1E-3 1E-1 1E0 1E1 1E3]; % ytickaxes
sp_h4(1).YTickLabel = {'10^{-3}','10^{-1}','10^{0}','10^{1}','10^{3}'};
sp_h4(1).FontSize = 12;%12;
% sp_h4(1).XTickLabel = legenda.parameters(x_ticks_keq_all);%x_ticks_keq_names;
sp_h4(1).XTickLabel = x_ticks_kcat_vm_names;
xtickangle(45)


%% plots Kcat,Vm-specific. Dotplot.
if exist('f_h4b','var')
    clf(f_h4b)
end
% f_h4 = figure(20004);
f_h4b = figure('Name','Kcat, Vm, barplot comparison');

% plotting Keq
sp_h4b(1) = subplot(1,1,1);
for i = 1:length(x_ticks_kcat_vm_all)
        color=customColormap(int8(1+127*(max(abs(x(x_ticks_kcat_vm_all(i))))/3)),:);
        semilogy(i,10.^x(x_ticks_kcat_vm_all(i)),'o','MarkerFaceColor',color,'MarkerEdgeColor',color);
        hold on
        sl_h = semilogy(i,10.^x(x_ticks_kcat_vm_all(i)),'o','MarkerFaceColor',color,'MarkerEdgeColor',color);
end
%%
% commong editing
sp_h4b(1).XLim = [0 length(x_ticks_kcat_vm_all)+1]; % ylims
sp_h4b(1).YLim = [10^-3 10^3]; % xlims
sp_h4b(1).XTick = 1:length(x_ticks_kcat_vm_all); % xtickaxes 
sp_h4b(1).YTick = [1E-3 1E-1 1E0 1E1 1E3]; % ytickaxes
sp_h4b(1).YTickLabel = {'10^{-3}','10^{-1}','10^{0}','10^{1}','10^{3}'};
sp_h4b(1).FontSize = 12;%12;
sp_h4b(1).XTickLabel = x_ticks_kcat_vm_names;
xtickangle(45)


%% Saving plots
% fh_5001
fh_1 = figure(1);
fh_2 = figure(2);
fh_3 = figure(3);
fh_4 = figure(4);
fh_5 = figure(5);
%%
% % % % % fh_1
% % % % set(fh_1,'Units','Inches');
% % % % pos = get(fh_1,'Position');
% % % % set(fh_1,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
% % % % print(fh_1,'E:\OneDrive - TU Eindhoven\Documents\ch4_Y3M1\results\figures\manuscript\SM_parameters_all','-dpdf','-r0')
% % % % % fh_2
% % % % set(fh_2,'Units','Inches');
% % % % pos = get(fh_2,'Position');
% % % % set(fh_2,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
% % % % print(fh_2,'E:\OneDrive - TU Eindhoven\Documents\ch4_Y3M1\results\figures\manuscript\SM_parameters_keq','-dpdf','-r0')
% % % % % fh_3
% % % % set(fh_3,'Units','Inches');
% % % % pos = get(fh_3,'Position');
% % % % set(fh_3,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
% % % % print(fh_3,'E:\OneDrive - TU Eindhoven\Documents\ch4_Y3M1\results\figures\manuscript\SM_parameters_km','-dpdf','-r0')
% % % % % fh_4
% % % % set(fh_4,'Units','Inches');
% % % % pos = get(fh_4,'Position');
% % % % set(fh_4,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
% % % % print(fh_4,'E:\OneDrive - TU Eindhoven\Documents\ch4_Y3M1\results\figures\manuscript\SM_parameters_vm','-dpdf','-r0')



%% memory safecopy
% % 
% fig_handle = figure('Name','Hello')
% fig_num = fig_handle.Number
% % 
% fig_num = double(fig_handle)


% % saving the plots
% savefig(1,'results\manuscript_supplementary_materials\supF_lit_params_all.fig')
% savefig(2,'results\manuscript_supplementary_materials\supF_lit_params_keq.fig')
% savefig(3,'results\manuscript_supplementary_materials\supF_lit_params_km.fig')
% savefig(4,'results\manuscript_supplementary_materials\supF_lit_params_kcat_vm.fig')











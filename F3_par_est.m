% % F3_PAR_EST.m
% This code reproduces the plots regarding the parameter estimation
% process.
% David Lao-Martil, 2022-05-23


%% (1) Plot part 1: regularization
if exist('fig_h1','var')
    close(fig_h1)
end
fig_h1 = figure('Name','regularization');
fig_h1.Position = [1965 670 355 295];
for plot_regularization = 1
    
    % loading SS D&C result enolase and metadata
    load('xENO_SSotput.mat');
    metadata;
    refferenceValues;
    setupParestENO; % ENO
    
    % setup options
    refpars     = setup.caseStudy.parameters;
    npars       = size(refpars);
    refvals     = setup.refVals.parameters(refpars);
    lambdalist  = setup.parEst.lambdalist;
    
    % setup options
    load('parameters_blank.mat');
    options = optimset('Display','iter');
    plength = length(setup.caseStudy.parameters);
    x_temp(plength)=0;
    lb = -3*ones(1,plength);
    ub =  3*ones(1,plength);
    ltot = length(setup.parEst.lambdalist);
    
    % data
    xres_array = outputENOss.regularization.xres_array;
    exitflag_array = outputENOss.regularization.exitflag_array;
    errorData_array = outputENOss.regularization.errorData_array;
    errorLambda_array = outputENOss.regularization.errorLambda_array;
    errorParams_array = outputENOss.regularization.errorParams_array;

    % plot
    % subplot(2,3,1)
    % regData
    [AX,H1,H2] = plotyy(lambdalist,errorParams_array,lambdalist,errorData_array,@loglog);
    left_color = [0, 0, 1];
    right_color = [1, 0, 0];
    set(H1,'Marker','o','Color','k','MarkerFaceColor',left_color)
    set(H2,'Marker','o','Color','k','MarkerFaceColor',right_color)
    set(AX(1),'yColor',left_color)
    set(AX(2),'yColor',right_color)
    set(AX,'defaultAxesColorOrder',[left_color; right_color]);
%     arrow('arrow',[0.5 0.5],[0 0.5])
    annot = annotation('arrow'); annot.X = [0.42 0.50]; annot.Y = [0.5 0.5];
    % lambda
    % lam = setup.parEst.lambda1;
    lam = 0.2;
    line([lam lam],[AX(1).YLim(1) AX(1).YLim(2)],'Color','black','LineStyle','--') 
    textlam = ['lambda = ', num2str(lam)];
    % style
    ylim = get(gca,'ylim');
    set(gca,'xtick',[])
    set(AX(1),'ytick',[])
    set(AX(2),'ytick',[])
    AX(1).XLim(2) = 5000;
    AX(2).XLim(2) = 5000;
    AX(1).YLim = [7 1000];
    AX(2).YLim = [0.1 100];
    ylabel(AX(1), 'Parameter Deviation');
    ylabel(AX(2), 'Model Error');
    xlabel('Regularization Factor')
    %title('regularization')
    
    box on
end
% 
% % % % print(gcf,'-dpdf','E:\OneDrive - TU Eindhoven\Documents\ch4_Y3M1\results\figures\manuscript\F3A');
%
figure

%% (2) Plot part 2: sampling (+regularization)
% plot ENO
for recall_data_plot_ENO = 1
    % recall data
    for i = 3
        % starting setup(1)
        load('xENO_SSotput.mat');
        metadata;
        refferenceValues;
        setupParestENO; % ENO
        solutions = outputENOss.MSparestNonReg.solutions;

        % starting setup(2)
        % early size selection
        ptot = length(setup.caseStudy.parameters);
        for i = 1:10
            j = i^2;
            %qcheck = exist q
            if j >= ptot
                %if(exist(q,'var'))
                exist q var;            % !! pretty ugly this part. Clean up needed. Use 'while' statement...
                k = ans; clear ans
                if k == 1   
                else
                    q = sqrt(j);
                end
            end
        end

        % getting values for the plot
        % parray = zeros(length(solutions),ptot); % commented to avoid leaving 0 values
        for i = 1:length(solutions)
            if solutions(i).Exitflag == 3
                    for j = 1:ptot
                        parray(i,j) = solutions(i).X(j);
                    end
            end
        end

        solutions2 = outputENOss.MSparestReg.solutions;
         for i = 1:length(solutions2)
            if solutions2(i).Exitflag == 3
                    for j = 1:ptot
                        parray2(i,j) = solutions2(i).X(j);
                    end
            end
         end
    end
end
% 
if exist('fig_h3','var')
    close(fig_h3)
end
% %%
fig_h3 = figure('Name','global_sampling_ENO_kcat');
fig_h3.Position = [2333 670 355 295];
for plot_ENO = 1
    i = 6;
    pv = 4;

    ylims = [0 50];
    area([-1 1], [ylims(2) ylims(2)],'LineStyle','none','FaceColor',[0.9 0.9 0.9])%,'FaceAlpha',0.5);
    hold on
    hh = histogram(parray(:,pv),[-3:0.15:3]);
    hh.FaceColor = [1 1 1];
    hold on
    hh2 = histogram(parray2(:,pv),[-3:0.15:3]);
    hh2.FaceColor = [0 0 0];
    axis([-3 3 0 50])

    set(gca,'ytick',[]), 

    xticks([-3 -1 0 1 3])
    leg_h1 = legend('accepted','reg.OFF','reg.ON'); % leg_h.Position = [0.7917    0.8144    0.1000    0.0883]
%     [leg_h1,ico] = legend('accepted','reg.OFF','reg.ON'); % leg_h.Position = [0.7917    0.8144    0.1000    0.0883]
%     leg_h1.Position(1) = 0.5175;
%     leg_h1.Position(2) = 0.8275;
    leg_h1.Position = [0.6871 0.7418 0.1786 0.1262];
%     leg_h1.Box = 'off';
    % sp3.YLim = [0 70];
%      icons_leg_h1(1).HandleVisibility = 'off';

    ylabel('Sample Number')
    xlabel('Deviation Literature Parameters (log-scale)')
%     title('histogram eno k_{cat}')
end
%
% % % % print(gcf,'-dpdf','E:\OneDrive - TU Eindhoven\Documents\ch4_Y3M1\results\figures\manuscript\F3B');

%%
% plot PFK
for recall_data_plot_PFK = 1
    % recall data
    for i = 11  
        % starting setup(1)
        load('xPFK_SSotput.mat');
        metadata;
        refferenceValues;
        setupParestPFK; % PFK
        solutions = outputPFKss.MSparestNonReg.solutions;

        % starting setup(2)
        % early size selection
        ptot = length(setup.caseStudy.parameters);
        for i = 1:10
            j = i^2;
            %qcheck = exist q
            if j >= ptot
                %if(exist(q,'var'))
                exist q var;            % !! pretty ugly this part. Clean up needed. Use 'while' statement...
                k = ans; clear ans
                if k == 1   
                else
                    q = sqrt(j);
                end
            end
        end

        % getting values for the plot
        % parray = zeros(length(solutions),ptot); % commented to avoid leaving 0 values
        for i = 1:length(solutions)
            if solutions(i).Exitflag == 3
                    for j = 1:ptot
                        parray(i,j) = solutions(i).X(j);
                    end
            end
        end

        solutions2 = outputPFKss.MSparestReg.solutions;
        for i = 1:length(solutions2)
            if solutions2(i).Exitflag == 3
                for j = 1:ptot
                    parray2(i,j) = solutions2(i).X(j);
                end
            end
         end   
        % 
    end
end
% 
if exist('fig_h4','var')
    close(fig_h4)
end
fig_h4 = figure('Name','global_sampling_PFK_kF6P');
fig_h4.Position = [2700 670 355 295];
% plotting
for plot_PFK = 1
    i = 20;
    % plotting
    pv = i-10;

%     sp3 = subplot(2,3,3);

    ylims = [0 50];
    area([-1 1], [ylims(2) ylims(2)],'LineStyle','none','FaceColor',[0.9 0.9 0.9])%,'FaceAlpha',0.5);
    hold on
    hh = histogram(parray(:,pv),[-3:0.15:3]);
    hh.FaceColor = [1 1 1];
    hold on
    hh2 = histogram(parray2(:,pv),[-3:0.15:3]);
    hh2.FaceColor = [0 0 0];
    axis([-3 3 0 50])
    set(gca,'ytick',[]), 

    xticks([-3 -1 0 1 3])
    leg_h2 = legend('accepted','reg.OFF','reg.ON'); % leg_h.Position = [0.7917    0.8144    0.1000    0.0883]
%     leg_h2.Position(1) = 0.7995;
%     leg_h2.Position(2) = 0.8275;
    leg_h2.Position = [0.6871 0.7418 0.1786 0.1262];
    % sp3.YLim = [0 70];

    ylabel('Sample Number')
    xlabel('Deviation Literature Parameters (log-scale)')
%     title('histogram pfk k_{f6p}')

end
%
% % % % print(gcf,'-dpdf','E:\OneDrive - TU Eindhoven\Documents\ch4_Y3M1\results\figures\manuscript\F3C');

%% (4/4) Plot part 3: weighting data types
clear xlim ylim
load('mF3sims.mat','error_PEPss','error_PEPgp','simRes_mF3','legenda','canelas_SS','data','setup','xDummy','namesHits')

%
if exist('fig_h5','var')
    close(fig_h5)
end
fig_h5 = figure('Name','pareto_front');
fig_h5.Position = [1965 275 355 295];
%
for plot_pareto = 1
    pt = plot(error_PEPss(7),error_PEPgp(7),'ko',...
        'MarkerSize', 10, 'LineWidth', 2); %,'markerfacecolor','red')
    % 
    hold on
    pt2 = plot(error_PEPss,error_PEPgp,'ko-');
    % 
    pt2.Color = [.5 .5 .5];
    pt2.MarkerFaceColor = [.7 .7 .7];
    % hold on
    xlabel('Model Error, SS_{PEP}')
    ylabel('Model Error, GP_{PEP}')
    xlim([-1 50])
    ylim([-1 50])
    % ticks
    xticks([0 10 20 30 40 50])
%     xticks([0 25 50])
    yticks([0 10 20 30 40 50])
    % 
    annot = annotation('arrow'); annot.X = [0.31 0.21]; annot.Y = [0.2475 0.2475];
    % 
    set(gca,'xtick',[])
    set(gca,'ytick',[])
    % 
end
%
% % % % print(gcf,'-dpdf','E:\OneDrive - TU Eindhoven\Documents\ch4_Y3M1\results\figures\manuscript\F3D');

%%
setup.plotResultsMode = 30;
setup.namesHits = namesHits(1,:);
[~] = plotAll_Y3M1(simRes_mF3(1,:),legenda,canelas_SS,data,setup,xDummy);
%%
% 
fig_SS_concentrations_h = figure(6);
d = fig_SS_concentrations_h.Children(22).Children(1).XData;
PEP_SS_exp = fig_SS_concentrations_h.Children(22).Children(1).YData;
PEP_SS_exp_negDelta = fig_SS_concentrations_h.Children(22).Children(1).YNegativeDelta;
PEP_SS_exp_posDelta = fig_SS_concentrations_h.Children(22).Children(1).YPositiveDelta;
PEP_SS_badFit = fig_SS_concentrations_h.Children(22).Children(2).YData;
PEP_SS_goodFit = fig_SS_concentrations_h.Children(22).Children(14).YData;
% PEP_SS_selectedFit = fig_SS_concentrations_h.Children(22).Children(8).YData;
PEP_SS_selectedFit = fig_SS_concentrations_h.Children(22).Children(15).YData;

% 
fig_GP_concentrations_h = figure(8);
time_exp = fig_GP_concentrations_h.Children(22).Children(1).XData;
PEP_GP_exp = fig_GP_concentrations_h.Children(22).Children(1).YData;
time_sim = fig_GP_concentrations_h.Children(22).Children(2).XData;
PEP_GP_goodFit = fig_GP_concentrations_h.Children(22).Children(2).YData;
PEP_GP_badFit = fig_GP_concentrations_h.Children(22).Children(14).YData;
% PEP_GP_selectedFit = fig_GP_concentrations_h.Children(22).Children(8).YData;
PEP_GP_selectedFit = fig_GP_concentrations_h.Children(22).Children(15).YData;
% 
close(6:11)


%% putting it back in the plot
% 
if exist('fig_h6','var')
    close(fig_h6)
end
fig_h6 = figure('Name','fit_SS_PEP');
fig_h6.Position = [2333 275 355 295];
for fit_SS_PEP = 1
    % GPbetter
    pt1 = plot(d,PEP_SS_badFit,'mo-');
    pt1.MarkerFaceColor = 'm';
    pt1.MarkerSize = 5;
    hold on
    % SSbetter
    pt2 = plot(d,PEP_SS_goodFit,'bo-');
    pt2.MarkerFaceColor = 'b';
    pt2.MarkerSize = 5;
    % selected
    pt3 = plot(d,PEP_SS_selectedFit,'ko-');
    pt3.MarkerFaceColor = 'k';
    pt3.MarkerSize = 5;
    % SSdata
    % plot(d,PEP_SS_exp,'r+')
    pt4 = errorbar(d,PEP_SS_exp,PEP_SS_exp_posDelta);
    pt4.Color = 'red';
    pt4.LineStyle = ':';
    % style
    legend('Model error,GP','Model error,SS','optimum','experimental')
    % title('fits SS PEP')
    xlim([0 0.4])
    ylim([0 3])
    xticks([0 0.1 0.2 0.3 0.4])
    yticks([0 1 2 3])
    % 
    set(gca,'xtick',[])
    set(gca,'ytick',[])
    xlabel('Dilution Rate (h^{-1})')
    ylabel('Concentration PEP_{SS}(mM)')
end
% 
% % % % print(gcf,'-dpdf','E:\OneDrive - TU Eindhoven\Documents\ch4_Y3M1\results\figures\manuscript\F3E');


%%
% 
if exist('fig_h7','var')
    close(fig_h7)
end
fig_h7 = figure('Name','fit_GP_PEP');
fig_h7.Position = [2700 275 355 295];
for fit_GP_PEP = 1
    % GPbetter
    pt1 = plot(time_sim,PEP_GP_goodFit,'m-');
    pt1.MarkerFaceColor = 'm';
    pt1.LineWidth = 1;
    hold on
    % SSbetter
    pt2 = plot(time_sim,PEP_GP_badFit,'b-');
    pt2.MarkerFaceColor = 'b';
    pt2.LineWidth = 1;
    % selected
    pt3 = plot(time_sim,PEP_GP_selectedFit,'k-');
    pt3.MarkerFaceColor = 'k';
    pt3.LineWidth = 1;
    % GP data
    pt4 = plot([-200,time_exp],[PEP_GP_exp(1), PEP_GP_exp],'r.:');
    pt4.MarkerSize = 10;
    % style
    legend('Model error,GP','Model error,SS','optimum','experimental')
    xlim([-100 340])
%     ylim([0 1.5])
    xticks([-100 0 100 200 300])
    yticks([0 0.5 1 1.5])
    % 
    set(gca,'xtick',[])
    set(gca,'ytick',[])
    xlabel('Time (s)')
    ylabel('Concentration PEP_{GP}(mM)')
end
% 
% % % % print(gcf,'-dpdf','E:\OneDrive - TU Eindhoven\Documents\ch4_Y3M1\results\figures\manuscript\F3F');


% % saving the plots
% savefig(1,'results\F3_regularization.fig')
% savefig(3,'results\F3_global_sampling_ENO_kcat.fig')
% savefig(4,'results\F3_global_sampling_PFK_kF6P.fig')
% savefig(5,'results\F3_pareto_front.fig')
% savefig(6,'results\F3_fit_SS_PEP.fig')
% savefig(7,'results\F3_fit_GP_PEP.fig')



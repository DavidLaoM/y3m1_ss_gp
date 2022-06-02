% % F6_LIT_PARAM.m
% This figure shows comparison with the parameters from the literature.
% David Lao-Martil, 2022-05-23


%% ENO: setup
% clear
set_paths
%     saveloc = [folder, '\results\variables\divide_and_conquer\xENO_SSotput.mat'];
%     savefigloc  = [folder, '\results\figures\ENO_SSstudy'];
dbstop if error
CanelasData2;
vHeerdenData;
metadata;
refferenceValues;
rng(1)
setupParestENO; % ENO
load('parameters_blank.mat');
Y3M1manus_labels;


% %% ENO: plots

% First plot
load('outputENOss.mat','outputENOss');
% mkdir(savefigloc);
setup.parEst.MSnonReg = 1;
setup.parEst.MSReg = 0;
setup.parEst.nonreg = 1;
setup.parEst.reg = 0;
%
tempdnCF = setup.parEst.drawnowCF;
setup.parEst.drawnowCF      = 1;
if setup.parEst.drawnowCF == 1
    setup.parEst.lambda = setup.parEst.lambda0;
    xres = outputENOss.costfunCheck.xres;
%     xres = zeros(size(xres));
    [error] = costfunENO(xres,canelas_SS,setup,x);

end
setup.parEst.drawnowCF      = tempdnCF;
clear tempdnCF

f_h = figure(1000);
% f2 = figure(101);
% ax2 = copyobj(f_h.Children(4),f2);
ENOsim_estPars = f_h.Children(4).Children(2).XData;
ENOexp_estPars = f_h.Children(4).Children(2).YData;


% Second plot
load('outputENOss.mat','outputENOss');
% mkdir(savefigloc);
setup.parEst.MSnonReg = 1;
setup.parEst.MSReg = 0;
setup.parEst.nonreg = 1;
setup.parEst.reg = 0;
%
tempdnCF = setup.parEst.drawnowCF;
setup.parEst.drawnowCF      = 1;
if setup.parEst.drawnowCF == 1
    setup.parEst.lambda = setup.parEst.lambda0;
    xres = outputENOss.costfunCheck.xres;
    xres = zeros(size(xres));
    [error] = costfunENO(xres,canelas_SS,setup,x);

end
setup.parEst.drawnowCF      = tempdnCF;
clear tempdnCF

f_h = figure(1000);
% f2 = figure(201);
% ax2 = copyobj(f_h.Children(4),f2);
ENOsim_litPars = f_h.Children(4).Children(2).XData;
ENOexp_litPars = f_h.Children(4).Children(2).YData;


%%
% figure
fig100 = figure(100);

% ENO
% 
maxVal = 1.1 * ceil(ENOexp_estPars(end)); 
sp8 = subplot(3,3,8);
line([0 maxVal],[0 maxVal],'Color','red','LineStyle','--')
hold on
% 
plot(ENOsim_estPars,ENOexp_estPars,...
    'o','MarkerSize',5,...
    'MarkerEdgeColor','k','MarkerFaceColor','k')
% 
plot(ENOsim_litPars,ENOexp_litPars,...
    'o','MarkerSize',5,...
    'MarkerEdgeColor','k','MarkerFaceColor',[.7 .7 .7])
% 
xlim([0 maxVal])
ylim([0 maxVal])
% title('ENO')
text(sp8.XLim(2) * 0.075, sp8.YLim(2) * 0.925,'ENO')
box on
% 
xlab = xlabel('simulated reaction rate (mmol L^{-1} s^{-1})');
% ylabel('experimental reaction rate (mmol L^{-1} s^{-1})')
% sp8.FontSize = 8;
xticks([0 1.5 3]), xlim([0 3.3])
yticks([0 1.5 3]), ylim([0 3.3])

% Last edits
% fig100.Position = [1921   257   1536   747];
fig100.Position = [1950 477 1463 497];

%% PGI
rng(1)
setupParestPGI;
load('parameters_blank.mat');

% early setup
load('outputPGIss.mat','outputPGIss');
% 
setup.parEst.MSnonReg = 1;
setup.parEst.MSReg = 0;
setup.parEst.nonreg = 1;
setup.parEst.reg = 0;
% 
tempdnCF = setup.parEst.drawnowCF;
setup.parEst.drawnowCF      = 1;
if setup.parEst.drawnowCF == 1
    setup.parEst.lambda = setup.parEst.lambda0;
    xres = outputPGIss.costfunCheck.xres;
    [error] = costfunPGI(xres,canelas_SS,setup,x);
end
setup.parEst.drawnowCF      = tempdnCF;
clear tempdnCF
% 
f_h = figure(1000);
PGIsim_estPars = f_h.Children(4).Children(2).XData;
PGIexp_estPars = f_h.Children(4).Children(2).YData;
% 
tempdnCF = setup.parEst.drawnowCF;
setup.parEst.drawnowCF      = 1;
if setup.parEst.drawnowCF == 1
    setup.parEst.lambda = setup.parEst.lambda0;
    xres = outputPGIss.costfunCheck.xres;
    xres = zeros(size(xres));
    [error] = costfunPGI(xres,canelas_SS,setup,x);
end
setup.parEst.drawnowCF      = tempdnCF;
clear tempdnCF
% 
f_h = figure(1000);
PGIsim_litPars = f_h.Children(4).Children(2).XData;
PGIexp_litPars = f_h.Children(4).Children(2).YData;

%% 
fig100 = figure(100);

% PGI
% 
maxVal = 1.1 * ceil(PGIexp_estPars(end)); 
sp2 = subplot(3,3,2);
line([0 maxVal],[0 maxVal],'Color','r','LineStyle','--')
hold on
% 
plot(PGIsim_estPars,PGIexp_estPars,...
    'o','MarkerSize',5,...
    'MarkerEdgeColor','k','MarkerFaceColor','k')
% 
plot(PGIsim_litPars,PGIexp_litPars,...
    'o','MarkerSize',5,...
    'MarkerEdgeColor','k','MarkerFaceColor',[.7 .7 .7])
% 
xlim([0 maxVal])
ylim([0 maxVal])
% title('PGI')
text(sp2.XLim(2) * 0.075, sp2.YLim(2) * 0.925, 'PGI')
box on
xticks([0 1 2]), xlim([0 2.2])
yticks([0 1 2]), ylim([0 2.2])

%% PFK
rng(1)
setupParestPFK; % PFK
load('parameters_blank.mat');


% early setup
load('outputPFKss.mat','outputPFKss');
setup.parEst.MSnonReg = 1;
setup.parEst.MSReg = 0;
setup.parEst.nonreg = 1;
setup.parEst.reg = 0;

% 
tempdnCF = setup.parEst.drawnowCF;
setup.parEst.drawnowCF      = 1;
if setup.parEst.drawnowCF == 1
    setup.parEst.lambda = setup.parEst.lambda0;
    xres = outputPFKss.costfunCheck.xres;
    [error] = costfunPFK(xres,canelas_SS,setup,x);
end
setup.parEst.drawnowCF      = tempdnCF;
clear tempdnCF
% 
f_h = figure(1000);
PFKsim_estPars = f_h.Children(4).Children(2).XData;
PFKexp_estPars = f_h.Children(4).Children(2).YData;
% 
tempdnCF = setup.parEst.drawnowCF;
setup.parEst.drawnowCF      = 1;
if setup.parEst.drawnowCF == 1
    setup.parEst.lambda = setup.parEst.lambda0;
    xres = outputPFKss.costfunCheck.xres;
    xres = zeros(size(xres));
    [error] = costfunPFK(xres,canelas_SS,setup,x);
end
setup.parEst.drawnowCF      = tempdnCF;
clear tempdnCF
% 
f_h = figure(1000);
PFKsim_litPars = f_h.Children(4).Children(2).XData;
PFKexp_litPars = f_h.Children(4).Children(2).YData;

%% 
fig100 = figure(100);

% 
maxVal = 1.1 * ceil(PFKexp_estPars(end)); 
sp3 = subplot(3,3,3);
line([0 maxVal],[0 maxVal],'Color','r','LineStyle','--')
hold on
% 
plot(PFKsim_estPars,PFKexp_estPars,...
    'o','MarkerSize',5,...
    'MarkerEdgeColor','k','MarkerFaceColor','k')
% 
plot(PFKsim_litPars,PFKexp_litPars,...
    'o','MarkerSize',5,...
    'MarkerEdgeColor','k','MarkerFaceColor',[.7 .7 .7])
% 
xlim([0 maxVal])
ylim([0 maxVal])
% title('PFK')
text(sp3.XLim(2) * 0.075, sp3.YLim(2) * 0.925,'PFK')
box on
xticks([0 1 2]), xlim([0 2.2])
yticks([0 1 2]), ylim([0 2.2])

%% TPI
rng(1)
setupParestTPI; % TPI
load('parameters_blank.mat');

% early setup
load('outputTPIss.mat','outputTPIss');
setup.parEst.MSnonReg = 1;
setup.parEst.MSReg = 0;
setup.parEst.nonreg = 1;
setup.parEst.reg = 0;

% estpars
tempdnCF = setup.parEst.drawnowCF;
setup.parEst.drawnowCF      = 1;
if setup.parEst.drawnowCF == 1
    setup.parEst.lambda = setup.parEst.lambda0;
    xres = outputTPIss.costfunCheck.xres;
    [error] = costfunTPI(xres,canelas_SS,setup,x);
end
setup.parEst.drawnowCF      = tempdnCF;
clear tempdnCF
% 
f_h = figure(1000);
TPIsim_estPars = f_h.Children(4).Children(2).XData;
TPIexp_estPars = f_h.Children(4).Children(2).YData;

% litpars
tempdnCF = setup.parEst.drawnowCF;
setup.parEst.drawnowCF      = 1;
if setup.parEst.drawnowCF == 1
    setup.parEst.lambda = setup.parEst.lambda0;
    xres = outputTPIss.costfunCheck.xres;
    xres = zeros(size(xres));
    [error] = costfunTPI(xres,canelas_SS,setup,x);
end
setup.parEst.drawnowCF      = tempdnCF;
clear tempdnCF
% 
f_h = figure(1000);
TPIsim_litPars = f_h.Children(4).Children(2).XData;
TPIexp_litPars = f_h.Children(4).Children(2).YData;

%% 
fig100 = figure(100);
% 
maxVal = 1.1 * ceil(TPIexp_estPars(end)); 
sp5 = subplot(3,3,5);
line([0 maxVal],[0 maxVal],'Color','r','LineStyle','--')
hold on
% 
plot(TPIsim_estPars,TPIexp_estPars,...
    'o','MarkerSize',5,...
    'MarkerEdgeColor','k','MarkerFaceColor','k')
% 
plot(TPIsim_litPars,TPIexp_litPars,...
    'o','MarkerSize',5,...
    'MarkerEdgeColor','k','MarkerFaceColor',[.7 .7 .7])
% 
xlim([0 maxVal])
ylim([0 maxVal])
% title('TPI')
text(sp5.XLim(2) * 0.075, sp5.YLim(2) * 0.925,'TPI')
box on
xticks([0 1 2]), xlim([0 2.2])
yticks([0 1 2]), ylim([0 2.2])

%% PYK
rng(1)
setupParestPYK; % PYK
load('parameters_blank.mat');

% 
load('outputPYKss.mat','outputPYKss');
setup.parEst.MSnonReg = 1;
setup.parEst.MSReg = 0;
setup.parEst.nonreg = 1;
setup.parEst.reg = 0;
% 
tempdnCF = setup.parEst.drawnowCF;
setup.parEst.drawnowCF      = 1;
if setup.parEst.drawnowCF == 1
    setup.parEst.lambda = setup.parEst.lambda0;
    xres = outputPYKss.costfunCheck.xres;
    [error] = costfunPYK(xres,canelas_SS,setup,x);
end
setup.parEst.drawnowCF      = tempdnCF;
clear tempdnCF
% 
f_h = figure(1000);
PYKsim_estPars = f_h.Children(4).Children(2).XData;
PYKexp_estPars = f_h.Children(4).Children(2).YData;

% 
load('outputPYKss.mat','outputPYKss');
setup.parEst.MSnonReg = 1;
setup.parEst.MSReg = 0;
setup.parEst.nonreg = 1;
setup.parEst.reg = 0;
% 
tempdnCF = setup.parEst.drawnowCF;
setup.parEst.drawnowCF      = 1;
if setup.parEst.drawnowCF == 1
    setup.parEst.lambda = setup.parEst.lambda0;
    xres = outputPYKss.costfunCheck.xres;
    xres = zeros(size(xres));
    [error] = costfunPYK(xres,canelas_SS,setup,x);
end
setup.parEst.drawnowCF      = tempdnCF;
clear tempdnCF
% 
f_h = figure(1000);
PYKsim_litPars = f_h.Children(4).Children(2).XData;
PYKexp_litPars = f_h.Children(4).Children(2).YData;

%% 
fig100 = figure(100);
% 
maxVal = 1.1 * ceil(PYKexp_estPars(end)); 
sp9 = subplot(3,3,9);
line([0 maxVal],[0 maxVal],'Color','r','LineStyle','--')
hold on
% 
plot(PYKsim_estPars,PYKexp_estPars,...
    'o','MarkerSize',5,...
    'MarkerEdgeColor','k','MarkerFaceColor','k')
% 
plot(PYKsim_litPars,PYKexp_litPars,...
    'o','MarkerSize',5,...
    'MarkerEdgeColor','k','MarkerFaceColor',[.7 .7 .7])
% 
xlim([0 maxVal])
ylim([0 maxVal])
% title('PYK')
text(sp9.XLim(2) * 0.075, sp9.YLim(2) * 0.925,'PYK')
box on
xticks([0 1.5 3]), xlim([0 3.3])
yticks([0 1.5 3]), ylim([0 3.3])

%% PGM
rng(1)
setupParestPGM; % pgm
load('parameters_blank.mat');

% early setup
load('outputPGMss.mat','outputPGMss');
setup.parEst.MSnonReg = 1;
setup.parEst.MSReg = 0;
setup.parEst.nonreg = 1;
setup.parEst.reg = 0;

% 
tempdnCF = setup.parEst.drawnowCF;
setup.parEst.drawnowCF      = 1;
if setup.parEst.drawnowCF == 1
    setup.parEst.lambda = setup.parEst.lambda0;
    xres = outputPGMss.costfunCheck.xres;
    [error] = costfunPGM(xres,canelas_SS,setup,x);
end
setup.parEst.drawnowCF      = tempdnCF;
clear tempdnCF
% 
f_h = figure(1000);
PGMsim_estPars = f_h.Children(4).Children(2).XData;
PGMexp_estPars = f_h.Children(4).Children(2).YData;

% early setup
load('outputPGMss.mat','outputPGMss');
setup.parEst.MSnonReg = 1;
setup.parEst.MSReg = 0;
setup.parEst.nonreg = 1;
setup.parEst.reg = 0;

% 
tempdnCF = setup.parEst.drawnowCF;
setup.parEst.drawnowCF      = 1;
if setup.parEst.drawnowCF == 1
    setup.parEst.lambda = setup.parEst.lambda0;
    xres = outputPGMss.costfunCheck.xres;
    xres = zeros(size(xres));
    [error] = costfunPGM(xres,canelas_SS,setup,x);
end
setup.parEst.drawnowCF      = tempdnCF;
clear tempdnCF
% 
f_h = figure(1000);
PGMsim_litPars = f_h.Children(4).Children(2).XData;
PGMexp_litPars = f_h.Children(4).Children(2).YData;

%% 
fig100 = figure(100);
% 
maxVal = 1.1 * ceil(PGMexp_estPars(end)); 
sp7 = subplot(3,3,7);
line([0 maxVal],[0 maxVal],'Color','r','LineStyle','--')
hold on
% 
plot(PGMsim_estPars,PGMexp_estPars,...
    'o','MarkerSize',5,...
    'MarkerEdgeColor','k','MarkerFaceColor','k')
% 
plot(PGMsim_litPars,PGMexp_litPars,...
    'o','MarkerSize',5,...
    'MarkerEdgeColor','k','MarkerFaceColor',[.7 .7 .7])
% 
xlim([0 maxVal])
ylim([0 maxVal])
% title('PGM')
text(sp7.XLim(2) * 0.075, sp7.YLim(2) * 0.925,'PGM')
box on
xticks([0 1.5 3]), xlim([0 3.3])
yticks([0 1.5 3]), ylim([0 3.3])

%% ALD/FBA
rng(1)
setupParestFBA; % FBA
load('parameters_blank.mat');

% early setup
load('outputFBAss.mat','outputFBAss');
setup.parEst.MSnonReg = 1;
setup.parEst.MSReg = 0;
setup.parEst.nonreg = 1;
setup.parEst.reg = 0;

% 
tempdnCF = setup.parEst.drawnowCF;
setup.parEst.drawnowCF      = 1;
if setup.parEst.drawnowCF == 1
    setup.parEst.lambda = setup.parEst.lambda0;
    xres = outputFBAss.costfunCheck.xres;
    [error] = costfunFBA(xres,canelas_SS,setup,x);
end
setup.parEst.drawnowCF      = tempdnCF;
clear tempdnCF
% 
f_h = figure(1000);
FBAsim_estPars = f_h.Children(4).Children(2).XData;
FBAexp_estPars = f_h.Children(4).Children(2).YData;

% 
load('outputFBAss.mat','outputFBAss');
setup.parEst.MSnonReg = 1;
setup.parEst.MSReg = 0;
setup.parEst.nonreg = 1;
setup.parEst.reg = 0;

% 
tempdnCF = setup.parEst.drawnowCF;
setup.parEst.drawnowCF      = 1;
if setup.parEst.drawnowCF == 1
    setup.parEst.lambda = setup.parEst.lambda0;
    xres = outputFBAss.costfunCheck.xres;
    xres = zeros(size(xres));
    [error] = costfunFBA(xres,canelas_SS,setup,x);
end
setup.parEst.drawnowCF      = tempdnCF;
clear tempdnCF
% 
f_h = figure(1000);
FBAsim_litPars = f_h.Children(4).Children(2).XData;
FBAexp_litPars = f_h.Children(4).Children(2).YData;

%% 
fig100 = figure(100);
% 
maxVal = 1.1 * ceil(FBAexp_estPars(end)); 
sp4 = subplot(3,3,4);
line([0 maxVal],[0 maxVal],'Color','r','LineStyle','--')
hold on
% 
plot(FBAsim_estPars,FBAexp_estPars,...
    'o','MarkerSize',5,...
    'MarkerEdgeColor','k','MarkerFaceColor','k')
% 
plot(FBAsim_litPars,FBAexp_litPars,...
    'o','MarkerSize',5,...
    'MarkerEdgeColor','k','MarkerFaceColor',[.7 .7 .7])
% 
xlim([0 maxVal])
ylim([0 maxVal])
% title('ALD')
text(sp4.XLim(2) * 0.075, sp4.YLim(2) * 0.925, 'ALD')
box on
xticks([0 1 2]), xlim([0 2.2])
yticks([0 1 2]), ylim([0 2.2])

%% GLT + HXK
rng(1)
setupParestGLT_HXK_old;
load('parameters_blank.mat');
% adjusted literature setup
setup.adj_lit = 0;

% early setup
load('outputGLT_HXKss_adjLit_oldcheck.mat','outputGLT_HXKss_adjLit_oldcheck');    
setup.parEst.MSnonReg = 1;
setup.parEst.MSReg = 0;
setup.parEst.nonreg = 1;
setup.parEst.reg = 0;

%     costfunCheck: [1×1 struct]
tempdnCF = setup.parEst.drawnowCF;
setup.parEst.drawnowCF      = 1;
if setup.parEst.drawnowCF == 1
    ltot = outputGLT_HXKss_adjLit_oldcheck.regularization.ltot;
    lambda1 = setup.parEst.lambda1;
    k = 1;
    xres = outputGLT_HXKss_adjLit_oldcheck.regularization.xres_array(k,:);
    [error] = costfunGLT_HXK(xres,canelas_SS,setup,x);
end
setup.parEst.drawnowCF      = tempdnCF;
clear tempdnCF
% 
f_h = figure(1000);
GLTsim_estPars = f_h.Children(3).Children(2).XData;
GLTexp_estPars = f_h.Children(3).Children(2).YData;
HXKsim_estPars = f_h.Children(2).Children(2).XData;
HXKexp_estPars = f_h.Children(2).Children(2).YData;

%     costfunCheck: [1×1 struct]
tempdnCF = setup.parEst.drawnowCF;
setup.parEst.drawnowCF      = 1;
if setup.parEst.drawnowCF == 1
    ltot = outputGLT_HXKss_adjLit_oldcheck.regularization.ltot;
    lambda1 = setup.parEst.lambda1;
    k = 1;
    xres = outputGLT_HXKss_adjLit_oldcheck.regularization.xres_array(k,:);
    xres = zeros(size(xres));
    [error] = costfunGLT_HXK(xres,canelas_SS,setup,x);
end
setup.parEst.drawnowCF      = tempdnCF;
clear tempdnCF
% 
f_h = figure(1000);
GLTsim_litPars = f_h.Children(3).Children(2).XData;
GLTexp_litPars = f_h.Children(3).Children(2).YData;
HXKsim_litPars = f_h.Children(2).Children(2).XData;
HXKexp_litPars = f_h.Children(2).Children(2).YData;


%% 
fig100 = figure(100);
% % % % % GLT
% % % % maxVal = 1.1 * ceil(GLTexp_estPars(end)); 
% % % % subplot(3,4,8)
% % % % line([0 maxVal],[0 maxVal],'Color','r','LineStyle','--')
% % % % hold on
% % % % % 
% % % % plot(GLTsim_estPars,GLTexp_estPars,...
% % % %     'o','MarkerSize',5,...
% % % %     'MarkerEdgeColor','k','MarkerFaceColor','k')
% % % % % 
% % % % plot(GLTsim_litPars,GLTexp_litPars,...
% % % %     'o','MarkerSize',5,...
% % % %     'MarkerEdgeColor','k','MarkerFaceColor',[.7 .7 .7])
% % % % % 
% % % % xlim([0 maxVal])
% % % % ylim([0 maxVal])
% % % % title('GLT')
% % % % box on
% HXK
maxVal = 1.1 * ceil(HXKexp_estPars(end)); 
sp1 = subplot(3,3,1);
line([0 maxVal],[0 maxVal],'Color','r','LineStyle','--')
hold on
% 
plot(HXKsim_estPars,HXKexp_estPars,...
    'o','MarkerSize',5,...
    'MarkerEdgeColor','k','MarkerFaceColor','k')
% 
plot(HXKsim_litPars,HXKexp_litPars,...
    'o','MarkerSize',5,...
    'MarkerEdgeColor','k','MarkerFaceColor',[.7 .7 .7])
% 
xlim([0 maxVal])
ylim([0 maxVal])
% title('HXK')
text(sp1.XLim(2) * 0.075, sp1.YLim(2) * 0.925,'HXK')
box on
xticks([0 1 2]), xlim([0 2.2])
yticks([0 1 2]), ylim([0 2.2])

%% GAPDH + PGK
rng(1)
setupParestGAPDH_PGK;
load('parameters_blank.mat');

% early setup
load('outputGAPDH_PGKss.mat','outputGAPDH_PGKss');
setup.parEst.MSnonReg = 1;
setup.parEst.MSReg = 0;
setup.parEst.nonreg = 1;
setup.parEst.reg = 0;

setup.parEst.MSnonReg = 1;
setup.parEst.MSReg = 0;
setup.parEst.nonreg = 1;
setup.parEst.reg = 0;

% 
tempdnCF = setup.parEst.drawnowCF;
setup.parEst.drawnowCF      = 1;
if setup.parEst.drawnowCF == 1
    setup.parEst.lambda = setup.parEst.lambda0;
    xres = outputGAPDH_PGKss.regularization.xres_array(1,:);
    [error] = costfunGAPDH_PGK(xres,canelas_SS,setup,x);
end
setup.parEst.drawnowCF      = tempdnCF;
clear tempdnCF
% 
f_h = figure(1000);
GAPDHsim_estPars = f_h.Children(3).Children(2).XData;
GAPDHexp_estPars = f_h.Children(3).Children(2).YData;
PGKsim_estPars = f_h.Children(2).Children(2).XData;
PGKexp_estPars = f_h.Children(2).Children(2).YData;
% 
tempdnCF = setup.parEst.drawnowCF;
setup.parEst.drawnowCF      = 1;
if setup.parEst.drawnowCF == 1
    setup.parEst.lambda = setup.parEst.lambda0;
    xres = outputGAPDH_PGKss.regularization.xres_array(1,:);
    xres = zeros(size(xres));
    [error] = costfunGAPDH_PGK(xres,canelas_SS,setup,x);
end
setup.parEst.drawnowCF      = tempdnCF;
clear tempdnCF
% 
f_h = figure(1000);
GAPDHsim_litPars = f_h.Children(3).Children(2).XData;
GAPDHexp_litPars = f_h.Children(3).Children(2).YData;
PGKsim_litPars = f_h.Children(2).Children(2).XData;
PGKexp_litPars = f_h.Children(2).Children(2).YData;


%% 
fig100 = figure(100);
% GAPDH
maxVal = 1.1 * ceil(GAPDHexp_estPars(end)); 
sp6 = subplot(3,3,6);
line([0 maxVal],[0 maxVal],'Color','r','LineStyle','--')
hold on
% 
plot(GAPDHsim_estPars,GAPDHexp_estPars,...
    'o','MarkerSize',5,...
    'MarkerEdgeColor','k','MarkerFaceColor','k')
% 
plot(GAPDHsim_litPars,GAPDHexp_litPars,...
    'o','MarkerSize',5,...
    'MarkerEdgeColor','k','MarkerFaceColor',[.7 .7 .7])
% 
xlim([0 maxVal])
ylim([0 maxVal])
% title('GAPDH')
text(sp6.XLim(2) * 0.075, sp6.YLim(2) * 0.925,'GAPDH')
box on
ylab = ylabel('experimental reaction rate (mmol L^{-1} s^{-1})');
ylab.Position = [-0.4072 3.8956 -1.0000];
xticks([0 1.5 3]), xlim([0 3.3])
yticks([0 1.5 3]), ylim([0 3.3])

% % % % % PGK
% % % % maxVal = 1.1 * ceil(PGKexp_estPars(end)); 
% % % % subplot(3,4,11)
% % % % line([0 maxVal],[0 maxVal],'Color','r','LineStyle','--')
% % % % hold on
% % % % % 
% % % % plot(PGKsim_estPars,PGKexp_estPars,...
% % % %     'o','MarkerSize',5,...
% % % %     'MarkerEdgeColor','k','MarkerFaceColor','k')
% % % % % 
% % % % plot(PGKsim_litPars,PGKexp_litPars,...
% % % %     'o','MarkerSize',5,...
% % % %     'MarkerEdgeColor','k','MarkerFaceColor',[.7 .7 .7])
% % % % % 
% % % % xlim([0 maxVal])
% % % % ylim([0 maxVal])
% % % % title('PGK')
% % % % box on


% %% PDC (+ ADH)
% rng(1)
% setupParestPDC_ADH_branches;
% load('parameters_blank.mat');
% 
% % early setup
% load('outputPDC_ADH_branchesss.mat','outputPDC_ADH_branchesss');
% % 
% setup.parEst.MSnonReg = 1;
% setup.parEst.MSReg = 0;
% setup.parEst.nonreg = 1;
% setup.parEst.reg = 0;
% 
% % 
% tempdnCF = setup.parEst.drawnowCF;
% setup.parEst.drawnowCF      = 1;
% if setup.parEst.drawnowCF == 1
%     setup.parEst.lambda = setup.parEst.lambda0;
% %     xres = outputPDC_ADH_branchesss.costfunCheck.xres;
%     xres = outputPDC_ADH_branchesss.regularization.xres_array(1,:);
% %     [error] = costfunPDC_ADH_branches(xres,canelas_SS,setup,x);
%     [error] = costfunPDC_ADH_BRANCHES(xres,canelas_SS,setup,x);
% end
% setup.parEst.drawnowCF      = tempdnCF;
% clear tempdnCF
% % 
% f_h = figure(1000);
% PDCsim_estPars = f_h.Children(3).Children(2).XData;
% PDCexp_estPars = f_h.Children(3).Children(2).YData;
% % ADHsim_estPars = f_h.Children(3).Children(2).XData;
% % ADHexp_estPars = f_h.Children(3).Children(2).YData;
% 
% % 
% tempdnCF = setup.parEst.drawnowCF;
% setup.parEst.drawnowCF      = 1;
% if setup.parEst.drawnowCF == 1
%     setup.parEst.lambda = setup.parEst.lambda0;
% %     xres = outputPDC_ADH_branchesss.costfunCheck.xres;
%     xres = outputPDC_ADH_branchesss.regularization.xres_array(1,:);
%     xres = zeros(size(xres));
% %     [error] = costfunPDC_ADH_branches(xres,canelas_SS,setup,x);
%     [error] = costfunPDC_ADH_BRANCHES(xres,canelas_SS,setup,x);
% end
% setup.parEst.drawnowCF      = tempdnCF;
% clear tempdnCF
% % 
% f_h = figure(1000);
% PDCsim_litPars = f_h.Children(3).Children(2).XData;
% PDCexp_litPars = f_h.Children(3).Children(2).YData;
% % ADHsim_estPars = f_h.Children(3).Children(2).XData;
% % ADHexp_estPars = f_h.Children(3).Children(2).YData;
% 
% %% 
% fig100 = figure(100);
% % PDC
% maxVal = 1.1 * ceil(PDCexp_estPars(end)); 
% sp10 = subplot(3,3,10);
% line([0 maxVal],[0 maxVal],'Color','r','LineStyle','--')
% hold on
% % 
% plot(PDCsim_estPars,PDCexp_estPars,...
%     'o','MarkerSize',5,...
%     'MarkerEdgeColor','k','MarkerFaceColor','k')
% % 
% plot(PDCsim_litPars,PDCexp_litPars,...
%     'o','MarkerSize',5,...
%     'MarkerEdgeColor','k','MarkerFaceColor',[.7 .7 .7])
% % 
% xlim([0 maxVal])
% ylim([0 maxVal])
% % title('PDC')
% text(sp10.XLim(2) * 0.075, sp10.YLim(2) * 0.925,'PDC')
% box on
% xticks([0 1.5 3]), xlim([0 3.3])
% yticks([0 1.5 3]), ylim([0 3.3])



% % % % %% saving
% % % % savefig_loc = 'E:\OneDrive - TU Eindhoven\Documents\ch4_Y3M1\results\figures\manuscript\';
% % % % saveName = [savefig_loc, 'fig4_DC_litParameters']; savefig(fig100, saveName);
% % % % 
% % % % % savefig(1, 'tempfig21.fig')
% % % % % specs printing (method 3)
% % % % set(gcf,'Units','inches');
% % % % screenposition = get(gcf,'Position');
% % % % set(gcf,...
% % % %     'PaperPosition',[0 0 screenposition(3:4)],...
% % % %     'PaperSize',[screenposition(3:4)]);
% % % % print -dpdf -painters fig4_DC_litParameters
% % % % print -dpng -painters fig4_DC_litParameters

% % % % %%
% % % % fig100.Position = [1945 369 731 610];
% % % % print(gcf,'-dpdf','E:\OneDrive - TU Eindhoven\Documents\ch4_Y3M1\results\figures\manuscript\F6AI');


% % saving the plots
% savefig(100,'results\F6_lit_param_sims.fig')

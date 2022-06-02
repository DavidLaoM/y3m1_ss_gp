% function x=parEstPDC_ADH_branches
%% parameter estimation PDC_ADH_branches

%% (0) Initial setup
clc, clear, close all
set_paths
% % % %     saveloc = [folder, '\results\variables\divide_and_conquer\xPDC_ADH_branches_SSotput.mat'];
% % % %     savefigloc  = [folder, '\results\figures\PDC_ADH_branches_SSstudy'];
dbstop if error
CanelasData2;
vHeerdenData;
metadata;
refferenceValues;
rng(1)
setupParestPDC_ADH_branches;
load('parameters_blank.mat');
% setup.startParamsOFF = 0;
setup.treCycleON = 0;
setup.importPiON = 0;
setup.conditionsSS = 0;
setup.conditionsGP = 0;


%% Parameter Estimation
% % %% (4) Parameter estimation. Study objective funtion.
% % % parameter estimation using lsqnonlin
% % [xres,resnorm,residual,exitflag] = runLSQsimple(canelas_SS,setup);
% % outputPDC_ADH_branchesss.costfunCheck.xres        = xres;
% % outputPDC_ADH_branchesss.costfunCheck.resnorm     = resnorm;
% % outputPDC_ADH_branchesss.costfunCheck.residual    = residual;
% % outputPDC_ADH_branchesss.costfunCheck.exitflag    = exitflag;
% % 
% % 
% % %% (5) Parameter estimation. Non-regularized. Multi Start.
% % [b,fval,exitflag,output,solutions1] = MSparest(canelas_SS,setup);
% % outputPDC_ADH_branchesss.MSparestNonReg.b              = b;
% % outputPDC_ADH_branchesss.MSparestNonReg.fval           = fval;
% % outputPDC_ADH_branchesss.MSparestNonReg.exitflag       = exitflag;
% % outputPDC_ADH_branchesss.MSparestNonReg.output         = output;
% % outputPDC_ADH_branchesss.MSparestNonReg.solutions      = solutions1;
% % 
% %     
% % %% (6) Parameter estimation. Non-regularized. Identifiability. PLA.
% % [h,PLAresult] = runPLA(b,x,canelas_SS,setup);
% % outputPDC_ADH_branchesss.PLAnonReg                       = PLAresult;


%% (7) Parameter estimation. Regularization.
[h,output] = regularization(canelas_SS,setup);
outputPDC_ADH_branchesss.regularization           = output;


% % %% (8) Parameter estimation. Regularized. Multi Start.
% % setup.parEst.MSnonReg         = 0;
% % setup.parEst.MSReg            = 1;
% % [b,fval,exitflag,output,solutions2] = MSparest(canelas_SS,setup);
% % outputPDC_ADH_branchesss.MSparestReg.b              = b;
% % outputPDC_ADH_branchesss.MSparestReg.fval           = fval;
% % outputPDC_ADH_branchesss.MSparestReg.exitflag       = exitflag;
% % outputPDC_ADH_branchesss.MSparestReg.output         = output;
% % outputPDC_ADH_branchesss.MSparestReg.solutions      = solutions2;
% %   
% % 
% % %% (9) Parameter estimation. Regularized. Identifiability. PLA + mPLA
% % setup.parEst.nonreg     = 0;
% % setup.parEst.reg        = 1;
% % [h,PLAresult] = runPLA(b,x,canelas_SS,setup);
% % outputPDC_ADH_branchesss.PLAreg                       = PLAresult;


%% (10) Final save
% save(saveloc,'outputPDC_ADH_branchesss');
% % % % save('temp_outputPDC_ADH_branchesss.mat','outputPDC_ADH_branchesss')


%% visualization of results
setup.parEst.drawnow = 1; %1;

if setup.parEst.drawnow == 1
    % early setup
    load('outputPDC_ADH_branchesss.mat','outputPDC_ADH_branchesss');
% % % %     save(saveloc,'outputPDC_ADH_branchesss');
% % % %     mkdir(savefigloc);
    setup.parEst.MSnonReg = 1;
    setup.parEst.MSReg = 0;
    setup.parEst.nonreg = 1;
    setup.parEst.reg = 0;
    
%     %       costfunCheck: [1×1 struct] % the temporary CF is activated for
%     %       this run
%     tempdnCF = setup.parEst.drawnowCF;
%     setup.parEst.drawnowCF      = 1;
%     if setup.parEst.drawnowCF == 1
%         setup.parEst.lambda = setup.parEst.lambda0;
%         xres = outputPDC_ADH_branchesss.costfunCheck.xres;
%         [error] = costfunPDC_ADH_branches(xres,canelas_SS,setup,x);
%         savefig(1000,[savefigloc, '\PDC_ADH_branches_SScosfunNonreg']);
%     end
%     setup.parEst.drawnowCF      = tempdnCF;
%     clear tempdnCF
%     
%     %     MSparestNonReg: [1×1 struct]
%     setup.parEst.MSnonReg = 1;
%     setup.parEst.MSReg = 0;
%     [h] = histMS(solutions1,setup);
%     savefig(h,[savefigloc, '\PDC_ADH_branches_SSMSnonreg']);
%     
%     %          PLAnonReg: [1×1 struct]
%     [h] = plotPLA(outputPDC_ADH_branchesss,setup);
%     savefig(h,[savefigloc, '\PDC_ADH_branches_PLAnonreg.fig']);

    setup.parEst.MSnonReg = 0;
    setup.parEst.MSReg = 1;
    setup.parEst.nonreg = 0;
    setup.parEst.reg = 1;
    
    %     regularization: [1×1 struct]
    [h] = plotRegularization(outputPDC_ADH_branchesss,setup);
% % % %     savefig(h,[savefigloc, '\PDC_ADH_branches_regularization.fig']);
    
%     %        MSparestReg: [1×1 struct]
%     [h] = histMS(solutions2,setup);
%     savefig(h,[savefigloc, '\PDC_ADH_branches_SSMSreg']);
% 
%     %             PLAreg: [1×1 struct]
%     [h] = plotPLA(outputPDC_ADH_branchesss,setup);
%     savefig(h,[savefigloc, '\PDC_ADH_branches_PLAreg.fig']);
end


%% memoryDump
% %% (5) Parameter estimation. Non-regularized. Multi Start.
% setupParestPDC_ADH_branches;
% 
% % figure
% [b,fval,exitflag,output,solutions] = MSparest(canelas_SS,setup);
% save('D:\PROJECT\2 Working files\Chapter 3 Consensus model\yeast_consensus_v1\results\variables\PDC_ADH_branches_SSstudy\xPDC_ADH_branches_MSnonreg.mat','b','fval','exitflag','output','solutions')
% 
% % visualization in a histogram plot
% [h] = histMS(solutions,setup);
% savefig(h,'D:\PROJECT\2 Working files\Chapter 3 Consensus model\yeast_consensus_v1\results\figures\PDC_ADH_branches_SSstudy\PDC_ADH_branches_MSnonreg.fig')
% close all
% % b = [3.0000    2.9206    3.0000    3.0000    0.9139   -1.6944    2.9942   -0.0957    0.7140];    
%     
% %% (6) Parameter estimation. Non-regularized. Identifiability. PLA.
% setupParestPDC_ADH_branches;
% setup.parEst.nonreg     = 1;
% setup.parEst.reg        = 0;
% load('parameters_blank.mat');
% load('xPDC_ADH_branches_MSnonreg.mat');
% % b = zeros(1,13);
% [h] = runPLA(b,x,canelas_SS,setup);
% savefig(99,'D:\PROJECT\2 Working files\Chapter 3 Consensus model\yeast_consensus_v1\results\figures\PDC_ADH_branches_SSstudy\PDC_ADH_branches_PLAnonreg.fig');
% close all

% %% (7) Parameter estimation. Regularization.
% setupParestPDC_ADH_branches;
% %     setup.parEst.lambdalist     = [1E-5, 3E-5, 1E-4, 3E-4, ...
% %                                 1E-3, 3E-3, 1E-2, 3E-2, ...
% %                                 1E-1, 3E-1, 1E0, 3E0, ...
% %                                 1E1, 3E1, 1E2, 3E2, ...
% %                                 1E3, 3E3];
% saveloc     = [folder, '\results\variables\PDC_ADH_branches_SSstudy\xPDC_ADH_branches_outputReg.mat'];
% savefigloc  = [folder, '\results\figures\PDC_ADH_branches_SSstudy\PDC_ADH_branches_regularization.fig'];
% [h,output] = regularization(canelas_SS,setup);
% save(saveloc,'output')
% savefig(100,savefigloc)
% % x = output;
% close all

% %% (8) Parameter estimation. Regularized. Multi Start.
% setupParestPDC_ADH_branches;
% setup.parEst.MSnonReg         = 0;
% setup.parEst.MSReg            = 1;
% [b,fval,exitflag,output,solutions] = MSparest(canelas_SS,setup);
% save('D:\PROJECT\2 Working files\Chapter 3 Consensus model\yeast_consensus_v1\results\variables\PDC_ADH_branches_SSstudy\xPDC_ADH_branches_MSreg.mat','b','fval','exitflag','output','solutions')
% % b =
% 
% % [h] = histMSreg(solutions,setup);
% [h] = histMS(solutions,setup);
% savefig(h,'D:\PROJECT\2 Working files\Chapter 3 Consensus model\yeast_consensus_v1\results\figures\PDC_ADH_branches_SSstudy\PDC_ADH_branches_MSreg.fig')
% close all
%     
% %% (9) Parameter estimation. Regularized. Identifiability. PLA + mPLA
% setupParestPDC_ADH_branches;
% setup.parEst.nonreg     = 0;
% setup.parEst.reg        = 1;
% load('parameters_blank.mat');
% load('xPDC_ADH_branches_MSreg.mat'); 
% % b = zeros(1,13);
% [h] = runPLA(b,x,canelas_SS,setup);
% savefig(99,'D:\PROJECT\2 Working files\Chapter 3 Consensus model\yeast_consensus_v1\results\figures\PDC_ADH_branches_SSstudy\PDC_ADH_branches_PLAreg.fig')
% close all

%% memoryDump

%     A2 = 0.20 < A1 < 0.25; 
%     A2 = (0.20 < A1) && (A1 < 0.25);    % it seems that plots here had the same behavior
% %     A2 = 0.205 < A1;
% %     A3 = A2 < 0.25;
% %     A4 = [A1, A2, A3]; disp(A4);
%     A2 = A(A & )
%     A2 = A1(A1<0.25 & A1>0.205);
%         A2 = 
%         d
%         V.rand{1}(end,3)
%     A3 = [A1, A2]; disp(A3);    


    % Show several strong maxima
%     s
    % Show several moderate maxima
%     s
    % Display a similar behavior as the ideal fit
    
%         A1 = zeros(setup.subsystem.numEvals,1);
%         for i = 1:setup.subsystem.numEvals
%             A1(i) = V.rand{i}(end,3);
%         end
%         A2 = find(A1<0.25 & A1>0.20);%0.215);
%         setup.subsystem.iPGI = A2;


% % plotAll_GP_vH
% figure
% for i = 1:30
%     subplot(3,10,i)
%     plot(T,Y(:,i))
%     xlim([-100 340])
%     title(legenda.metabolites{i})
% end
% 
% figure
% for i = 1:42
%     subplot(4,11,i)
%     plot(T,V(:,i))
%     xlim([-100 340])
%     title(legenda.fluxes{i})
% end

% %% clustering
%     A1 = zeros(setup.subsystem.numEvals,1);
%     for i = 1:setup.subsystem.numEvals
%         A1(i) = V.rand{i}(end,3);
%     end
%     A2 = find(A1<0.25 & A1>0.20);%0.215);
%     setup.subsystem.iPGI = A2;

% %% figure color
% clf(2)
% figure(2)
% plot([0 1], [0 0.1], 'Color', [0 0 1])
% hold on
% plot([0 1], [0 0.2], 'Color', [0.25 0.75 1])


% % p_sol1 = output_x(:,8);
% p_sol1 = output_x(:,17);


% old cost function part

% % % select setup options
% % options = optimset('Display','iter');
% % plength = length(setup.indexParameters.PGI);
% % x_temp(plength)=0;
% % % x_temp=[0.0319    0.2053    3.0000    1.9139];
% % lb = -3*ones(1,plength);
% % ub =  3*ones(1,plength);
% % setup.parEst.lambda = setup.parEst.lambda0;
% % 
% % % run lsqnonlinwith objective function
% % [xres,resnorm,residual,exitflag] = lsqnonlin(@costfunPGI,x_temp,lb,ub,options,canelas_SS,setup,x);

%% parameter estimation GAPDH enzyme


%% (0) Initial setup
clc, clear, close all
set_paths
% % % %     saveloc = [folder, '\results\variables\divide_and_conquer\xGAPDH_PGK_SSotput.mat'];
% % % %     savefigloc  = [folder, '\results\figures\GAPDH_PGK_SSstudy'];
dbstop if error
CanelasData2;
vHeerdenData;
metadata;
refferenceValues;
rng(1)
setupParestGAPDH_PGK;
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
% % outputGAPDH_PGKss.costfunCheck.xres        = xres;
% % outputGAPDH_PGKss.costfunCheck.resnorm     = resnorm;
% % outputGAPDH_PGKss.costfunCheck.residual    = residual;
% % outputGAPDH_PGKss.costfunCheck.exitflag    = exitflag;
% % 
% % 
% % %% (5) Parameter estimation. Non-regularized. Multi Start.
% % [b,fval,exitflag,output,solutions1] = MSparest(canelas_SS,setup);
% % outputGAPDH_PGKss.MSparestNonReg.b              = b;
% % outputGAPDH_PGKss.MSparestNonReg.fval           = fval;
% % outputGAPDH_PGKss.MSparestNonReg.exitflag       = exitflag;
% % outputGAPDH_PGKss.MSparestNonReg.output         = output;
% % outputGAPDH_PGKss.MSparestNonReg.solutions      = solutions1;
% % 
% %     
% % %% (6) Parameter estimation. Non-regularized. Identifiability. PLA.
% % [h,PLAresult] = runPLA(b,x,canelas_SS,setup);
% % outputGAPDH_PGKss.PLAnonReg                       = PLAresult;


%% (7) Parameter estimation. Regularization.
[h,output] = regularization(canelas_SS,setup);
outputGAPDH_PGKss.regularization           = output;


% % %% (8) Parameter estimation. Regularized. Multi Start.
% % setup.parEst.MSnonReg         = 0;
% % setup.parEst.MSReg            = 1;
% % [b,fval,exitflag,output,solutions2] = MSparest(canelas_SS,setup);
% % outputGAPDH_PGKss.MSparestReg.b              = b;
% % outputGAPDH_PGKss.MSparestReg.fval           = fval;
% % outputGAPDH_PGKss.MSparestReg.exitflag       = exitflag;
% % outputGAPDH_PGKss.MSparestReg.output         = output;
% % outputGAPDH_PGKss.MSparestReg.solutions      = solutions2;
% %   
% % 
% % %% (9) Parameter estimation. Regularized. Identifiability. PLA + mPLA
% % setup.parEst.nonreg     = 0;
% % setup.parEst.reg        = 1;
% % [h,PLAresult] = runPLA(b,x,canelas_SS,setup);
% % outputGAPDH_PGKss.PLAreg                       = PLAresult;


%% (10) Final save
% save(saveloc,'outputGAPDH_PGKss');
% % % % save('outputGAPDH_PGKss.mat','outputGAPDH_PGKss')


%% visualization of results
setup.parEst.drawnow = 1; %1;

if setup.parEst.drawnow == 1
    % early setup
    load('outputGAPDH_PGKss.mat','outputGAPDH_PGKss');
% % % %     save(saveloc,'outputGAPDH_PGKss');
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
%         xres = outputGAPDH_PGKss.costfunCheck.xres;
%         [error] = costfunGAPDH_PGK(xres,canelas_SS,setup,x);
%         savefig(1000,[savefigloc, '\GAPDH_PGK_SScosfunNonreg']);
%     end
%     setup.parEst.drawnowCF      = tempdnCF;
%     clear tempdnCF
    
%     %     MSparestNonReg: [1×1 struct]
%     setup.parEst.MSnonReg = 1;
%     setup.parEst.MSReg = 0;
%     [h] = histMS(solutions1,setup);
%     savefig(h,[savefigloc, '\GAPDH_PGK_SSMSnonreg']);
%     
%     %          PLAnonReg: [1×1 struct]
%     [h] = plotPLA(outputGAPDH_PGKss,setup);
%     savefig(h,[savefigloc, '\GAPDH_PGK_PLAnonreg.fig']);

    setup.parEst.MSnonReg = 0;
    setup.parEst.MSReg = 1;
    setup.parEst.nonreg = 0;
    setup.parEst.reg = 1;
    
    %     regularization: [1×1 struct]
    [h] = plotRegularization(outputGAPDH_PGKss,setup);
% % % %     savefig(h,[savefigloc, '\GAPDH_PGK_regularization.fig']);
    
%     %        MSparestReg: [1×1 struct]
%     [h] = histMS(solutions2,setup);
%     savefig(h,[savefigloc, '\GAPDH_PGK_SSMSreg']);
% 
%     %             PLAreg: [1×1 struct]
%     [h] = plotPLA(outputGAPDH_PGKss,setup);
%     savefig(h,[savefigloc, '\GAPDH_PGK_PLAreg.fig']);
end


%% Network analysis
% % %% (1a) Literature parameters. Dynamics GLT. sPSA
% % setupParestGLT; % HXK
% % [xdata, ydata] = sPSA(setup, canelas_SS);
% % [h] = plotSPSA(xdata,ydata,setup);
% % % savefig(1,'D:\PROJECT\2 Working files\Chapter 3 Consensus model\yeast_consensus_v1\results\figures\GAPDH_PGK_SSstudy\GLT_sPSA.fig')
% % 
% % 
% % %% (1b) Literature parameters. Dynamics HXK. sPSA
% % setupParestHXK; % HXK
% % [xdata, ydata] = sPSA(setup, canelas_SS);
% % [h] = plotSPSA(xdata,ydata,setup);
% % % savefig(1,'D:\PROJECT\2 Working files\Chapter 3 Consensus model\yeast_consensus_v1\results\figures\GAPDH_PGK_SSstudy\HXK_sPSA.fig')
% % 
% % 
% % %% (3a) Sensitivity. Steady state conditions. GLT
% % setupParestGLT; % HXK
% % [T,Y,V,parray] = subsystemDynamics(setup,canelas_SS,vHeerden_GP);
% % [ind] = clusteringSensitivity(T,Y,V,setup);
% % plotSubsystem(T,Y,V,legenda,setup,vHeerden_GP,parray,ind)
% % % savefig(4,'D:\PROJECT\2 Working files\Chapter 3 Consensus model\yeast_consensus_v1\results\figures\GAPDH_PGK_SSstudy\GLTsensitivity_dependencies.fig')
% % 
% % 
% % %% (3b) Sensitivity. Steady state conditions. HXK
% % setupParestHXK; % HXK
% % [T,Y,V,parray] = subsystemDynamics(setup,canelas_SS,vHeerden_GP);
% % [ind] = clusteringSensitivity(T,Y,V,setup);
% % plotSubsystem(T,Y,V,legenda,setup,vHeerden_GP,parray,ind)
% % % savefig(5,'D:\PROJECT\2 Working files\Chapter 3 Consensus model\yeast_consensus_v1\results\figures\GAPDH_PGK_SSstudy\HXKsensitivity_dependencies025.fig')
% % % savefig(4,'D:\PROJECT\2 Working files\Chapter 3 Consensus model\yeast_consensus_v1\results\figures\GAPDH_PGK_SSstudy\HXKsensitivity_dependencies_kcat025.fig')
% % 
% % 
% % % %% (4) Parameter estimation. Study objective funtion.
% % % 
% % % % parameter estimation using lsqnonlin
% % % [xres,resnorm,residual,exitflag] = runLSQsimple(canelas_SS,setup);
% % % % savefig(1000,'D:\PROJECT\2 Working files\Chapter 3 Consensus model\yeast_consensus_v1\results\figures\GAPDH_PGK_SSstudy\GAPDH_PGK_costfunCheck.fig')
% % % % save('D:\PROJECT\2 Working files\Chapter 3 Consensus model\yeast_consensus_v1\results\variables\GAPDH_PGK_SSstudy\xGAPDH_PGK_SSf.mat','xres')
% % % 
% % % 

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

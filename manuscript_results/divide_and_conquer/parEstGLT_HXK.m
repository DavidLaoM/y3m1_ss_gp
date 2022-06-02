%% parameter estimation GLT_HXK enzyme


%% (0) Initial setup
clc, clear, close all
set_paths
%     savefigloc  = [folder, '\results\figures\step1_divide_and_conquer\GLT_HXK_SSstudy'];
%     saveloc = [folder, '\results\variables\step1_divide_and_conquer\outputGLT_HXKss.mat'];
% % % %     savefigloc  = [folder, '\results\figures\step1_divide_and_conquer\GLT_HXK_SSstudy_adjLit'];
% % % %     saveloc = [folder, '\results\variables\step1_divide_and_conquer\outputGLT_HXKss_adjLit.mat'];
dbstop if error
CanelasData2;
vHeerdenData;
metadata;
refferenceValues;
rng(1)
setupParestGLT_HXK_old;
load('parameters_blank.mat');
% adjusted literature setup
setup.adj_lit = 1;
% setup.startParamsOFF = 0;
setup.treCycleON = 0;
setup.importPiON = 0;
setup.conditionsSS = 0;
setup.conditionsGP = 0;


% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
%% Parameter Estimation (low growth rates)
setup.parEst.costfun       = 1;


% % %% (4) Parameter estimation. Study objective funtion.
% % % parameter estimation using lsqnonlin
% % [xres,resnorm,residual,exitflag] = runLSQsimple(canelas_SS,setup);
% % outputGLT_HXK_lowMUss.costfunCheck.xres        = xres;
% % outputGLT_HXK_lowMUss.costfunCheck.resnorm     = resnorm;
% % outputGLT_HXK_lowMUss.costfunCheck.residual    = residual;
% % outputGLT_HXK_lowMUss.costfunCheck.exitflag    = exitflag;
% % 
% % 
% % %% (5) Parameter estimation. Non-regularized. Multi Start.
% % [b,fval,exitflag,output,solutions1] = MSparest(canelas_SS,setup);
% % outputGLT_HXK_lowMUss.MSparestNonReg.b              = b;
% % outputGLT_HXK_lowMUss.MSparestNonReg.fval           = fval;
% % outputGLT_HXK_lowMUss.MSparestNonReg.exitflag       = exitflag;
% % outputGLT_HXK_lowMUss.MSparestNonReg.output         = output;
% % outputGLT_HXK_lowMUss.MSparestNonReg.solutions      = solutions1;
% % 
% %     
% % %% (6) Parameter estimation. Non-regularized. Identifiability. PLA.
% % [h,PLAresult] = runPLA(b,x,canelas_SS,setup);
% % outputGLT_HXK_lowMUss.PLAnonReg                       = PLAresult;


%% (7) Parameter estimation. Regularization.
[h,output] = regularization(canelas_SS,setup);
outputGLT_HXKss.regularization           = output;


% % %% (8) Parameter estimation. Regularized. Multi Start.
% % setup.parEst.MSnonReg         = 0;
% % setup.parEst.MSReg            = 1;
% % [b,fval,exitflag,output,solutions2] = MSparest(canelas_SS,setup);
% % outputGLT_HXK_lowMUss.MSparestReg.b              = b;
% % outputGLT_HXK_lowMUss.MSparestReg.fval           = fval;
% % outputGLT_HXK_lowMUss.MSparestReg.exitflag       = exitflag;
% % outputGLT_HXK_lowMUss.MSparestReg.output         = output;
% % outputGLT_HXK_lowMUss.MSparestReg.solutions      = solutions2;
% %   
% % 
% % %% (9) Parameter estimation. Regularized. Identifiability. PLA + mPLA
% % setup.parEst.nonreg     = 0;
% % setup.parEst.reg        = 1;
% % [h,PLAresult] = runPLA(b,x,canelas_SS,setup);
% % outputGLT_HXK_lowMUss.PLAreg                       = PLAresult;


%% (10) Final save
outputGLT_HXKss_adjLit = outputGLT_HXKss; clear outputGLT_HXKss
% % % % mkdir(savefigloc);
% % % % save(saveloc,'outputGLT_HXKss_adjLit');
% save('outputGLT_HXK_lowMUss.mat','outputGLT_HXK_lowMUss');
% save('outputGLT_HXK_lowMU_14_58ss.mat','outputGLT_HXK_lowMUss');
% save('outputGLT_HXK_lowMU_13_48ss.mat','outputGLT_HXK_lowMUss');


%% visualization of results
setup.parEst.drawnow = 1; %1;
setup.parEst.costfun = 1;
setup.parEst.lambda1 = 0.1; %same profile for all growth rates
% setup.parEst.lambda1 = 0.01; %same profile for all growth rates
% setup.parEst.lambda1 = 0.00001; %same profile for all growth rates
% setup.parEst.lambda1 = 0.002;
% setup.parEst.lambda1 = 0.005; %14_58
% setup.parEst.lambda1 = 0.002; %13_48
setup.parEst.lambda = setup.parEst.lambda1;

if setup.parEst.drawnow == 1
    % early setup
    load('outputGLT_HXKss_adjLit.mat','outputGLT_HXKss_adjLit');    
% %     load('outputGLT_HXK_lowMUss.mat','outputGLT_HXK_lowMUss');
%     load('outputGLT_HXK_lowMU_14_58ss.mat','outputGLT_HXK_lowMUss');
% %     load('outputGLT_HXK_lowMU_13_48ss.mat','outputGLT_HXK_lowMUss');
%     save(saveloc,'outputGLT_HXK_lowMUss');
%     mkdir(savefigloc);
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
%         xres = outputGLT_HXK_lowMUss.costfunCheck.xres;
%         [error] = costfunGLT_HXK(xres,canelas_SS,setup,x);
%         savefig(1000,[savefigloc, '\GLT_HXK_lowMU_SScosfunNonreg']);
%     end
%     setup.parEst.drawnowCF      = tempdnCF;
%     clear tempdnCF
%     
%     %     MSparestNonReg: [1×1 struct]
%     setup.parEst.MSnonReg = 1;
%     setup.parEst.MSReg = 0;
%     [h] = histMS(solutions1,setup);
%     savefig(h,[savefigloc, '\GLT_HXK_lowMU_SSMSnonreg']);
%     
%     %          PLAnonReg: [1×1 struct]
%     [h] = plotPLA(outputGLT_HXK_lowMUss,setup);
%     savefig(h,[savefigloc, '\GLT_HXK_lowMU_PLAnonreg.fig']);

    setup.parEst.MSnonReg = 0;
    setup.parEst.MSReg = 1;
    setup.parEst.nonreg = 0;
    setup.parEst.reg = 1;
    
    %     regularization: [1×1 struct]
    [h] = plotRegularization(outputGLT_HXKss_adjLit,setup);
% % % %     savefig(h,[savefigloc, '\GLT_HXK_regularization_adjLit.fig']);
    
    %     costfunCheck: [1×1 struct]
    tempdnCF = setup.parEst.drawnowCF;
    setup.parEst.drawnowCF      = 1;
    if setup.parEst.drawnowCF == 1
        ltot = outputGLT_HXKss_adjLit.regularization.ltot;
        lambda1 = setup.parEst.lambda1;
        k = find(ltot==lambda1);
%         k = k - 3;
%         k = 1;
        xres = outputGLT_HXKss_adjLit.regularization.xres_array(k,:);
        [error] = costfunGLT_HXK(xres,canelas_SS,setup,x);
% % % %         savefig(1000,[savefigloc, '\GLT_HXK_SScosfunNonreg_adjLit.fig']);
    end
    setup.parEst.drawnowCF      = tempdnCF;
    clear tempdnCF
%     
%     %        MSparestReg: [1×1 struct]
%     [h] = histMS(solutions2,setup);
%     savefig(h,[savefigloc, '\GLT_HXK_lowMU_SSMSreg']);
%     
%     %             PLAreg: [1×1 struct]
%     [h] = plotPLA(outputGLT_HXK_lowMUss,setup);
%     savefig(h,[savefigloc, '\GLT_HXK_lowMU_PLAreg.fig']);
%     
end

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% % % % %% Parameter Estimation (high growth rates)
% % % % setup.parEst.costfun       = 2;
% % % % 
% % % % 
% % % % % % %% (4) Parameter estimation. Study objective funtion.
% % % % % % % parameter estimation using lsqnonlin
% % % % % % [xres,resnorm,residual,exitflag] = runLSQsimple(canelas_SS,setup);
% % % % % % outputGLT_HXK_highMUss.costfunCheck.xres        = xres;
% % % % % % outputGLT_HXK_highMUss.costfunCheck.resnorm     = resnorm;
% % % % % % outputGLT_HXK_highMUss.costfunCheck.residual    = residual;
% % % % % % outputGLT_HXK_highMUss.costfunCheck.exitflag    = exitflag;
% % % % % % 
% % % % % % 
% % % % % % %% (5) Parameter estimation. Non-regularized. Multi Start.
% % % % % % [b,fval,exitflag,output,solutions1] = MSparest(canelas_SS,setup);
% % % % % % outputGLT_HXK_highMUss.MSparestNonReg.b              = b;
% % % % % % outputGLT_HXK_highMUss.MSparestNonReg.fval           = fval;
% % % % % % outputGLT_HXK_highMUss.MSparestNonReg.exitflag       = exitflag;
% % % % % % outputGLT_HXK_highMUss.MSparestNonReg.output         = output;
% % % % % % outputGLT_HXK_highMUss.MSparestNonReg.solutions      = solutions1;
% % % % % % 
% % % % % %     
% % % % % % %% (6) Parameter estimation. Non-regularized. Identifiability. PLA.
% % % % % % [h,PLAresult] = runPLA(b,x,canelas_SS,setup);
% % % % % % outputGLT_HXK_highMUss.PLAnonReg                       = PLAresult;
% % % % 
% % % % 
% % % % %% (7) Parameter estimation. Regularization.
% % % % [h,output] = regularization(canelas_SS,setup);
% % % % outputGLT_HXK_highMUss.regularization           = output;
% % % % 
% % % % 
% % % % % % %% (8) Parameter estimation. Regularized. Multi Start.
% % % % % % setup.parEst.MSnonReg         = 0;
% % % % % % setup.parEst.MSReg            = 1;
% % % % % % [b,fval,exitflag,output,solutions2] = MSparest(canelas_SS,setup);
% % % % % % outputGLT_HXK_highMUss.MSparestReg.b              = b;
% % % % % % outputGLT_HXK_highMUss.MSparestReg.fval           = fval;
% % % % % % outputGLT_HXK_highMUss.MSparestReg.exitflag       = exitflag;
% % % % % % outputGLT_HXK_highMUss.MSparestReg.output         = output;
% % % % % % outputGLT_HXK_highMUss.MSparestReg.solutions      = solutions2;
% % % % % %   
% % % % % % 
% % % % % % %% (9) Parameter estimation. Regularized. Identifiability. PLA + mPLA
% % % % % % setup.parEst.nonreg     = 0;
% % % % % % setup.parEst.reg        = 1;
% % % % % % [h,PLAresult] = runPLA(b,x,canelas_SS,setup);
% % % % % % outputGLT_HXK_highMUss.PLAreg                       = PLAresult;
% % % % 
% % % % 
% % % % %% (10) Final save
% % % % % save(saveloc,'outputGLT_HXKss');
% % % % save('outputGLT_HXK_highMUss.mat','outputGLT_HXK_highMUss')
% % % % % save('outputGLT_HXK_highMU_14_58ss.mat','outputGLT_HXK_highMUss')
% % % % % save('outputGLT_HXK_highMU_13_48ss.mat','outputGLT_HXK_highMUss')
% % % % 
% % % % 
% % % % %% visualization of results
% % % % setup.parEst.drawnow = 1;
% % % % setup.parEst.costfun = 2;
% % % % setup.parEst.lambda1 = 0.0005; %14_58
% % % % % setup.parEst.lambda1 = 0.002; %13_48
% % % % setup.parEst.lambda = setup.parEst.lambda1;
% % % % 
% % % % if setup.parEst.drawnow == 1
% % % %     % early setup
% % % % %     load('outputGLT_HXK_highMUss.mat','outputGLT_HXK_highMUss');
% % % %     load('outputGLT_HXK_highMU_14_58ss.mat','outputGLT_HXK_highMUss');
% % % % %     load('outputGLT_HXK_highMU_13_48ss.mat','outputGLT_HXK_highMUss');
% % % %     save(saveloc,'outputGLT_HXK_highMUss');
% % % %     mkdir(savefigloc);
% % % %     setup.parEst.MSnonReg = 1;
% % % %     setup.parEst.MSReg = 0;
% % % %     setup.parEst.nonreg = 1;
% % % %     setup.parEst.reg = 0;
% % % %     
% % % % %     %       costfunCheck: [1×1 struct] % the temporary CF is activated for
% % % % %     %       this run
% % % % %     tempdnCF = setup.parEst.drawnowCF;
% % % % %     setup.parEst.drawnowCF      = 1;
% % % % %     if setup.parEst.drawnowCF == 1
% % % % %         setup.parEst.lambda = setup.parEst.lambda0;
% % % % %         xres = outputGLT_HXK_highMUss.costfunCheck.xres;
% % % % %         [error] = costfunGLT_HXK(xres,canelas_SS,setup,x);
% % % % %         savefig(1000,[savefigloc, '\GLT_HXK_highMU_SScosfunNonreg']);
% % % % %     end
% % % % %     setup.parEst.drawnowCF      = tempdnCF;
% % % % %     clear tempdnCF
% % % % %     
% % % % %     %     MSparestNonReg: [1×1 struct]
% % % % %     setup.parEst.MSnonReg = 1;
% % % % %     setup.parEst.MSReg = 0;
% % % % %     [h] = histMS(solutions1,setup);
% % % % %     savefig(h,[savefigloc, '\GLT_HXK_highMU_SSMSnonreg']);
% % % % %     
% % % % %     %          PLAnonReg: [1×1 struct]
% % % % %     [h] = plotPLA(outputGLT_HXK_highMUss,setup);
% % % % %     savefig(h,[savefigloc, '\GLT_HXK_highMU_PLAnonreg.fig']);
% % % % 
% % % %     setup.parEst.MSnonReg = 0;
% % % %     setup.parEst.MSReg = 1;
% % % %     setup.parEst.nonreg = 0;
% % % %     setup.parEst.reg = 1;
% % % %     
% % % %     %     regularization: [1×1 struct]
% % % %     [h] = plotRegularization(outputGLT_HXK_highMUss,setup);
% % % %     savefig(h,[savefigloc, '\GLT_HXK_highMU_regularization.fig']);
% % % %     %     costfunCheck: [1×1 struct]
% % % %     tempdnCF = setup.parEst.drawnowCF;
% % % %     setup.parEst.drawnowCF      = 1;
% % % %     if setup.parEst.drawnowCF == 1
% % % %         ltot = outputGLT_HXK_highMUss.regularization.ltot;
% % % %         lambda1 = setup.parEst.lambda1;
% % % %         k = find(ltot==lambda1); 
% % % %         k = k + 6;
% % % %         xres = outputGLT_HXK_highMUss.regularization.xres_array(k,:);
% % % %         [error] = costfunGLT_HXK(xres,canelas_SS,setup,x);
% % % %         savefig(1000,[savefigloc, '\GLT_HXK_lowMU_SScosfunNonreg']);
% % % %     end
% % % %     setup.parEst.drawnowCF      = tempdnCF;
% % % %     clear tempdnCF
% % % %     
% % % %     
% % % % %     %        MSparestReg: [1×1 struct]
% % % % %     [h] = histMS(solutions2,setup);
% % % % %     savefig(h,[savefigloc, '\GLT_HXK_highMU_SSMSreg']);
% % % % % 
% % % % %     %             PLAreg: [1×1 struct]
% % % % %     [h] = plotPLA(outputGLT_HXK_highMUss,setup);
% % % % %     savefig(h,[savefigloc, '\GLT_HXK_highMU_PLAreg.fig']);
% % % % end
% % % % 
% % % % 
% % % % %% Network analysis
% % % % % %% (1a) Literature parameters. Dynamics GLT. sPSA
% % % % % setupParestGLT; % HXK
% % % % % [xdata, ydata] = sPSA(setup, canelas_SS);
% % % % % [h] = plotSPSA(xdata,ydata,setup);
% % % % % % savefig(1,'D:\PROJECT\2 Working files\Chapter 3 Consensus model\yeast_consensus_v1\results\figures\GLT_HXK_SSstudy\GLT_sPSA.fig')
% % % % % 
% % % % % 
% % % % % %% (1b) Literature parameters. Dynamics HXK. sPSA
% % % % % setupParestHXK; % HXK
% % % % % [xdata, ydata] = sPSA(setup, canelas_SS);
% % % % % [h] = plotSPSA(xdata,ydata,setup);
% % % % % % savefig(1,'D:\PROJECT\2 Working files\Chapter 3 Consensus model\yeast_consensus_v1\results\figures\GLT_HXK_SSstudy\HXK_sPSA.fig')
% % % % % 
% % % % % 
% % % % % %% (3a) Sensitivity. Steady state conditions. GLT
% % % % % setupParestGLT; % HXK
% % % % % [T,Y,V,parray] = subsystemDynamics(setup,canelas_SS,vHeerden_GP);
% % % % % [ind] = clusteringSensitivity(T,Y,V,setup);
% % % % % plotSubsystem(T,Y,V,legenda,setup,vHeerden_GP,parray,ind)
% % % % % % savefig(4,'D:\PROJECT\2 Working files\Chapter 3 Consensus model\yeast_consensus_v1\results\figures\GLT_HXK_SSstudy\GLTsensitivity_dependencies.fig')
% % % % % 
% % % % % 
% % % % % %% (3b) Sensitivity. Steady state conditions. HXK
% % % % % setupParestHXK; % HXK
% % % % % [T,Y,V,parray] = subsystemDynamics(setup,canelas_SS,vHeerden_GP);
% % % % % [ind] = clusteringSensitivity(T,Y,V,setup);
% % % % % plotSubsystem(T,Y,V,legenda,setup,vHeerden_GP,parray,ind)
% % % % % % savefig(5,'D:\PROJECT\2 Working files\Chapter 3 Consensus model\yeast_consensus_v1\results\figures\GLT_HXK_SSstudy\HXKsensitivity_dependencies025.fig')
% % % % % % savefig(4,'D:\PROJECT\2 Working files\Chapter 3 Consensus model\yeast_consensus_v1\results\figures\GLT_HXK_SSstudy\HXKsensitivity_dependencies_kcat025.fig')
% % % % % 
% % % % % 
% % % % % % %% (4) Parameter estimation. Study objective funtion.
% % % % % % 
% % % % % % % parameter estimation using lsqnonlin
% % % % % % [xres,resnorm,residual,exitflag] = runLSQsimple(canelas_SS,setup);
% % % % % % % savefig(1000,'D:\PROJECT\2 Working files\Chapter 3 Consensus model\yeast_consensus_v1\results\figures\GLT_HXK_SSstudy\GLT_HXK_costfunCheck.fig')
% % % % % % % save('D:\PROJECT\2 Working files\Chapter 3 Consensus model\yeast_consensus_v1\results\variables\GLT_HXK_SSstudy\xGLT_HXK_SSf.mat','xres')
% % % % 
% % % % 
% % % % %% memoryDump
% % % % 
% % % % %     A2 = 0.20 < A1 < 0.25; 
% % % % %     A2 = (0.20 < A1) && (A1 < 0.25);    % it seems that plots here had the same behavior
% % % % % %     A2 = 0.205 < A1;
% % % % % %     A3 = A2 < 0.25;
% % % % % %     A4 = [A1, A2, A3]; disp(A4);
% % % % %     A2 = A(A & )
% % % % %     A2 = A1(A1<0.25 & A1>0.205);
% % % % %         A2 = 
% % % % %         d
% % % % %         V.rand{1}(end,3)
% % % % %     A3 = [A1, A2]; disp(A3);    
% % % % 
% % % % 
% % % %     % Show several strong maxima
% % % % %     s
% % % %     % Show several moderate maxima
% % % % %     s
% % % %     % Display a similar behavior as the ideal fit
% % % %     
% % % % %         A1 = zeros(setup.subsystem.numEvals,1);
% % % % %         for i = 1:setup.subsystem.numEvals
% % % % %             A1(i) = V.rand{i}(end,3);
% % % % %         end
% % % % %         A2 = find(A1<0.25 & A1>0.20);%0.215);
% % % % %         setup.subsystem.iPGI = A2;
% % % % 
% % % % 
% % % % % % plotAll_GP_vH
% % % % % figure
% % % % % for i = 1:30
% % % % %     subplot(3,10,i)
% % % % %     plot(T,Y(:,i))
% % % % %     xlim([-100 340])
% % % % %     title(legenda.metabolites{i})
% % % % % end
% % % % % 
% % % % % figure
% % % % % for i = 1:42
% % % % %     subplot(4,11,i)
% % % % %     plot(T,V(:,i))
% % % % %     xlim([-100 340])
% % % % %     title(legenda.fluxes{i})
% % % % % end
% % % % 
% % % % % %% clustering
% % % % %     A1 = zeros(setup.subsystem.numEvals,1);
% % % % %     for i = 1:setup.subsystem.numEvals
% % % % %         A1(i) = V.rand{i}(end,3);
% % % % %     end
% % % % %     A2 = find(A1<0.25 & A1>0.20);%0.215);
% % % % %     setup.subsystem.iPGI = A2;
% % % % 
% % % % % %% figure color
% % % % % clf(2)
% % % % % figure(2)
% % % % % plot([0 1], [0 0.1], 'Color', [0 0 1])
% % % % % hold on
% % % % % plot([0 1], [0 0.2], 'Color', [0.25 0.75 1])
% % % % 
% % % % 
% % % % % % p_sol1 = output_x(:,8);
% % % % % p_sol1 = output_x(:,17);
% % % % 
% % % % 
% % % % % old cost function part
% % % % 
% % % % % % % select setup options
% % % % % % options = optimset('Display','iter');
% % % % % % plength = length(setup.indexParameters.PGI);
% % % % % % x_temp(plength)=0;
% % % % % % % x_temp=[0.0319    0.2053    3.0000    1.9139];
% % % % % % lb = -3*ones(1,plength);
% % % % % % ub =  3*ones(1,plength);
% % % % % % setup.parEst.lambda = setup.parEst.lambda0;
% % % % % % 
% % % % % % % run lsqnonlinwith objective function
% % % % % % [xres,resnorm,residual,exitflag] = lsqnonlin(@costfunPGI,x_temp,lb,ub,options,canelas_SS,setup,x);


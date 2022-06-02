%% parameter estimation PFK enzyme


%% (0) Initial setup
clc, clear, close all
set_paths
% % % %     saveloc = [folder, '\results\variables\divide_and_conquer\xPFK_SSotput.mat'];
% % % %     savefigloc  = [folder, '\results\figures\PFK_SSstudy'];
dbstop if error
CanelasData2;
vHeerdenData;
metadata;
refferenceValues;
rng(1)
setupParestPFK; % PFK
load('parameters_blank.mat');
% setup.startParamsOFF = 0;
setup.treCycleON = 0;
setup.importPiON = 0;
setup.conditionsSS = 0;
setup.conditionsGP = 0;



%% Parameter Estimation
%% (4) Parameter estimation. Study objective funtion.
% parameter estimation using lsqnonlin
[xres,resnorm,residual,exitflag] = runLSQsimple(canelas_SS,setup);
outputPFKss.costfunCheck.xres        = xres;
outputPFKss.costfunCheck.resnorm     = resnorm;
outputPFKss.costfunCheck.residual    = residual;
outputPFKss.costfunCheck.exitflag    = exitflag;


%% (5) Parameter estimation. Non-regularized. Multi Start.
[b,fval,exitflag,output,solutions1] = MSparest(canelas_SS,setup);
outputPFKss.MSparestNonReg.b              = b;
outputPFKss.MSparestNonReg.fval           = fval;
outputPFKss.MSparestNonReg.exitflag       = exitflag;
outputPFKss.MSparestNonReg.output         = output;
outputPFKss.MSparestNonReg.solutions      = solutions1;

    
%% (6) Parameter estimation. Non-regularized. Identifiability. PLA.
[h,PLAresult] = runPLA(b,x,canelas_SS,setup);
outputPFKss.PLAnonReg                       = PLAresult;


%% (7) Parameter estimation. Regularization.
[h,output] = regularization(canelas_SS,setup);
outputPFKss.regularization           = output;


%% (8) Parameter estimation. Regularized. Multi Start.
setup.parEst.MSnonReg         = 0;
setup.parEst.MSReg            = 1;
[b,fval,exitflag,output,solutions2] = MSparest(canelas_SS,setup);
outputPFKss.MSparestReg.b              = b;
outputPFKss.MSparestReg.fval           = fval;
outputPFKss.MSparestReg.exitflag       = exitflag;
outputPFKss.MSparestReg.output         = output;
outputPFKss.MSparestReg.solutions      = solutions2;
  

%% (9) Parameter estimation. Regularized. Identifiability. PLA + mPLA
setup.parEst.nonreg     = 0;
setup.parEst.reg        = 1;
[h,PLAresult] = runPLA(b,x,canelas_SS,setup);
outputPFKss.PLAreg                       = PLAresult;


%% (10) Final save
% save(saveloc,'outputPFKss');
% % % % save('outputPFKss.mat','outputPFKss')


%% visualization of results
setup.parEst.drawnow = 1; %1;

if setup.parEst.drawnow == 1
    % early setup
    load('outputPFKss.mat','outputPFKss');
% % % %     save(saveloc,'outputPFKss');
% % % %     mkdir(savefigloc);
    setup.parEst.MSnonReg = 1;
    setup.parEst.MSReg = 0;
    setup.parEst.nonreg = 1;
    setup.parEst.reg = 0;
    
    %       costfunCheck: [1×1 struct] % the temporary CF is activated for
    %       this run
    tempdnCF = setup.parEst.drawnowCF;
    setup.parEst.drawnowCF      = 1;
    if setup.parEst.drawnowCF == 1
        setup.parEst.lambda = setup.parEst.lambda0;
        xres = outputPFKss.costfunCheck.xres;
        [error] = costfunPFK(xres,canelas_SS,setup,x);
% % % %         savefig(1000,[savefigloc, '\PFK_SScosfunNonreg']);
    end
    setup.parEst.drawnowCF      = tempdnCF;
    clear tempdnCF
    
    %     MSparestNonReg: [1×1 struct]
    setup.parEst.MSnonReg = 1;
    setup.parEst.MSReg = 0;
    solutions1 = outputPFKss.MSparestNonReg.solutions;
    [h] = histMS(solutions1,setup);
% % % %     savefig(h,[savefigloc, '\PFK_SSMSnonreg']);
    
    %          PLAnonReg: [1×1 struct]
    [h] = plotPLA(outputPFKss,setup);
% % % %     savefig(h,[savefigloc, '\PFK_PLAnonreg.fig']);

    setup.parEst.MSnonReg = 0;
    setup.parEst.MSReg = 1;
    setup.parEst.nonreg = 0;
    setup.parEst.reg = 1;
    
    %     regularization: [1×1 struct]
    [h] = plotRegularization(outputPFKss,setup);
% % % %     savefig(h,[savefigloc, '\PFK_regularization.fig']);
    
    %        MSparestReg: [1×1 struct]
    solutions2 = outputPFKss.MSparestReg.solutions;
    [h] = histMS(solutions2,setup);
% % % %     savefig(h,[savefigloc, '\PFK_SSMSreg']);

    %             PLAreg: [1×1 struct]
    [h] = plotPLA(outputPFKss,setup);
% % % %     savefig(h,[savefigloc, '\PFK_PLAreg.fig']);
end


%% Network analysis
% %% (1) Literature Parameters. Dynamics
% [xdata, ydata] = sPSA(setup, canelas_SS);
% [h] = plotSPSA(xdata,ydata,setup);
% % savefig(1,'D:\PROJECT\2 Working files\Chapter 3 Consensus model\yeast_consensus_v1\results\figures\PFK_SSstudy\PFK_sPSA.fig')
% 
% 
% %% (2) Sensitivity in steady conditions. Random parameter sampling. Analysis of parameter effect and dependencies
% 
% 
% %% (3) Sensitivity in dynamic conditions. Random parameter sampling. Analysis of parameter effect and dependencies
% [T,Y,V,parray] = subsystemDynamics(setup,canelas_SS,vHeerden_GP);
% [ind] = clusteringSensitivity(T,Y,V,setup);
% plotSubsystem(T,Y,V,legenda,setup,vHeerden_GP,parray,ind)
% % savefig(4,'D:\PROJECT\2 Working files\Chapter 3 Consensus model\yeast_consensus_v1\results\figures\PFK_SSstudy\PFK_GPsensitivities.fig')
% 
% 
% % %% (4) Parameter estimation. Study objective funtion.
% % 
% % % parameter estimation using lsqnonlin
% % [xres,resnorm,residual,exitflag] = runLSQsimple(canelas_SS,setup);
% % % savefig(1000,'D:\PROJECT\2 Working files\Chapter 3 Consensus model\yeast_consensus_v1\results\figures\PFK_SSstudy\PFK_costfunCheck.fig')
% % % save('D:\PROJECT\2 Working files\Chapter 3 Consensus model\yeast_consensus_v1\results\variables\PFK_SSstudy\xPFK_SSf.mat','xres')
% % 
% % 


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

%% parameter estimation ALD/FBA enzyme
    % (0) Initial setup
    % (1) ParEst pipeline
    % (2) Visualization
    % (3) Sensitivity

%% (0) Initial setup
clc, clear, close all
set_paths
% % % %     saveloc = [folder, '\results\variables\divide_and_conquer\xFBA_SSotput.mat'];
% % % %     savefigloc  = [folder, '\results\figures\FBA_SSstudy'];
dbstop if error
CanelasData2;
vHeerdenData;
metadata;
refferenceValues;
rng(1)
setupParestFBA; % FBA
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
outputFBAss.costfunCheck.xres        = xres;
outputFBAss.costfunCheck.resnorm     = resnorm;
outputFBAss.costfunCheck.residual    = residual;
outputFBAss.costfunCheck.exitflag    = exitflag;


%% (5) Parameter estimation. Non-regularized. Multi Start.
[b,fval,exitflag,output,solutions1] = MSparest(canelas_SS,setup);
outputFBAss.MSparestNonReg.b              = b;
outputFBAss.MSparestNonReg.fval           = fval;
outputFBAss.MSparestNonReg.exitflag       = exitflag;
outputFBAss.MSparestNonReg.output         = output;
outputFBAss.MSparestNonReg.solutions      = solutions1;

    
%% (6) Parameter estimation. Non-regularized. Identifiability. PLA.
[h,PLAresult] = runPLA(b,x,canelas_SS,setup);
outputFBAss.PLAnonReg                       = PLAresult;


%% (7) Parameter estimation. Regularization.
[h,output] = regularization(canelas_SS,setup);
outputFBAss.regularization           = output;


%% (8) Parameter estimation. Regularized. Multi Start.
setup.parEst.MSnonReg         = 0;
setup.parEst.MSReg            = 1;
[b,fval,exitflag,output,solutions2] = MSparest(canelas_SS,setup);
outputFBAss.MSparestReg.b              = b;
outputFBAss.MSparestReg.fval           = fval;
outputFBAss.MSparestReg.exitflag       = exitflag;
outputFBAss.MSparestReg.output         = output;
outputFBAss.MSparestReg.solutions      = solutions2;
  

%% (9) Parameter estimation. Regularized. Identifiability. PLA + mPLA
setup.parEst.nonreg     = 0;
setup.parEst.reg        = 1;
[h,PLAresult] = runPLA(b,x,canelas_SS,setup);
outputFBAss.PLAreg                       = PLAresult;


%% (10) Final save
% save(saveloc,'outputFBAss');
% % % % save('outputFBAss.mat','outputFBAss')


%% visualization of results
setup.parEst.drawnow = 1; %1;

if setup.parEst.drawnow == 1
    % early setup
    load('outputFBAss.mat','outputFBAss');
% % % %     save(saveloc,'outputFBAss');
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
        xres = outputFBAss.costfunCheck.xres;
        [error] = costfunFBA(xres,canelas_SS,setup,x);
% % % %         savefig(1000,[savefigloc, '\FBA_SScosfunNonreg']);
    end
    setup.parEst.drawnowCF      = tempdnCF;
    clear tempdnCF
    
    %     MSparestNonReg: [1×1 struct]
    setup.parEst.MSnonReg = 1;
    setup.parEst.MSReg = 0;
    solutions1 = outputFBAss.MSparestNonReg.solutions;
    [h] = histMS(solutions1,setup);
% % % %     savefig(h,[savefigloc, '\FBA_SSMSnonreg']);
    
    %          PLAnonReg: [1×1 struct]
    [h] = plotPLA(outputFBAss,setup);
% % % %     savefig(h,[savefigloc, '\FBA_PLAnonreg.fig']);

    setup.parEst.MSnonReg = 0;
    setup.parEst.MSReg = 1;
    setup.parEst.nonreg = 0;
    setup.parEst.reg = 1;
    
    %     regularization: [1×1 struct]
    [h] = plotRegularization(outputFBAss,setup);
% % % %     savefig(h,[savefigloc, '\FBA_regularization.fig']);
    
    %        MSparestReg: [1×1 struct]
    solutions2 = outputFBAss.MSparestReg.solutions;
    [h] = histMS(solutions2,setup);
% % % %     savefig(h,[savefigloc, '\FBA_SSMSreg']);

    %             PLAreg: [1×1 struct]
    [h] = plotPLA(outputFBAss,setup);
% % % %     savefig(h,[savefigloc, '\FBA_PLAreg.fig']);
end


%% Network analysis
% %% (1) Literature Parameters. Dynamics
% % Sensitivity analysis of the different equation contents over the
% % simulation.
% %         setup.refVals.metabolites(12) = setup.refVals.metabolites(12) * 0.25; % product needs to be decreased for reaction in the right direction
% [xdata, ydata] = sPSA(setup, canelas_SS);
% [h] = plotSPSA(xdata,ydata,setup);
% % savefig(1,'D:\PROJECT\2 Working files\Chapter 3 Consensus model\yeast_consensus_v1\results\figures\FBA_SSstudy\FBA_sPSA.fig')
% 
% %!!! plot in a semilogarithmic manner, with all parameters within the same
% %order of difference of change (so they fall within the real examination
% %region, wqith -3 to +3 orders of magnitude).
% 
% 
% %% (2) Sensitivity in steady conditions. Random parameter sampling. Analysis of parameter effect and dependencies
% 
% 
% %% (3) Sensitivity in dynamic conditions. Random parameter sampling. Analysis of parameter effect and dependencies
% % A simulation of the GP is use to study the effect of the parameters.
% % In here, the model intially received has been used.
%     % !!! Here we are using the initla model provided by joep. Once the entire
%     % system ODE has been updated for this version, change it. 
% 
% % Simulation
% % (in subsystemDynamics)parray      = - 0.25 + (0.25+0.25) * rand(numEvals,npars);
% [T,Y,V,parray] = subsystemDynamics(setup,canelas_SS,vHeerden_GP);
% 
% % Manual clustering criteria (upon later visualization. Simulations that:
% [ind] = clusteringSensitivity(T,Y,V,setup);
% 
% % Plotting
% plotSubsystem(T,Y,V,legenda,setup,vHeerden_GP,parray,ind)
% % savefig(4,'D:\PROJECT\2 Working files\Chapter 3 Consensus model\yeast_consensus_v1\results\figures\FBA_SSstudy\FBA_GPsensitivity_dependencies_limits1035.fig')
% % savefig(1,'D:\PROJECT\2 Working files\Chapter 3 Consensus model\yeast_consensus_v1\results\figures\FBA_SSstudy\FBA_GPsensitivity_profile.fig')
% 
% 
% % % %% (4) Parameter estimation. Study objective funtion.
% % % 
% % % % parameter estimation using lsqnonlin
% % % [xres,resnorm,residual,exitflag] = runLSQsimple(canelas_SS,setup);
% % % % savefig(1000,'D:\PROJECT\2 Working files\Chapter 3 Consensus model\yeast_consensus_v1\results\figures\FBA_SSstudy\FBA_costfunCheck.fig')
% % % % save('D:\PROJECT\2 Working files\Chapter 3 Consensus model\yeast_consensus_v1\results\variables\FBA_SSstudy\xFBA_SSf.mat','xres')
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

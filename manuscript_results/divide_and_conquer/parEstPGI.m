%% parameter estimation PGI enzyme
% topology has been kept the same since Teusink model
% parameter values have been updated through time and now come from
% Smallbone, vanEunen and updated from Teusink

% in this script:
% 	Subsystem dynamics
% 	1. Literature parameters. Dynamics.
% 	2. Random parameter sampling: Reproduce the glucose pulse step with different parameter combinations.
% 	3. Random parameter sampling: Study dependencies and try to identify the most relevant parameter for the response.
% 	Parameter estimation
% 	4. Parameter estimation. Study objective funtion.
% 	5. Parameter estimation. Non-regularized. Multi Start.
% 	6. Parameter estimation. Non-regularized. Identifiability. PLA.
% 	7. Parameter estimation. Regularization.
% 	8. Parameter estimation. Regularized. Multi Start.
% 	9. Parameter estimation. Regularized. Identifiability. PLA + mPLA


%% (0) Initial setup
clc, clear, close all
% clc, clear, clf([1,2,3,4])%, close all
set_paths
% % % %     saveloc = [folder, '\results\variables\divide_and_conquer\xPGI_SSotput.mat'];
% % % %     savefigloc  = [folder, '\results\figures\PGI_SSstudy'];
dbstop if error
% load('parameters_blank.mat');
CanelasData2;
vHeerdenData;
metadata; % !!! add here the setup metadata thing
refferenceValues
rng(1)
setupParestPGI;
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
outputPGIss.costfunCheck.xres        = xres;
outputPGIss.costfunCheck.resnorm     = resnorm;
outputPGIss.costfunCheck.residual    = residual;
outputPGIss.costfunCheck.exitflag    = exitflag;


%% (5) Parameter estimation. Non-regularized. Multi Start.
[b,fval,exitflag,output,solutions1] = MSparest(canelas_SS,setup);
outputPGIss.MSparestNonReg.b              = b;
outputPGIss.MSparestNonReg.fval           = fval;
outputPGIss.MSparestNonReg.exitflag       = exitflag;
outputPGIss.MSparestNonReg.output         = output;
outputPGIss.MSparestNonReg.solutions      = solutions1;

    
%% (6) Parameter estimation. Non-regularized. Identifiability. PLA.
[h,PLAresult] = runPLA(b,x,canelas_SS,setup);
outputPGIss.PLAnonReg                       = PLAresult;


%% (7) Parameter estimation. Regularization.
[h,output] = regularization(canelas_SS,setup);
outputPGIss.regularization           = output;


%% (8) Parameter estimation. Regularized. Multi Start.
setup.parEst.MSnonReg         = 0;
setup.parEst.MSReg            = 1;
[b,fval,exitflag,output,solutions2] = MSparest(canelas_SS,setup);
outputPGIss.MSparestReg.b              = b;
outputPGIss.MSparestReg.fval           = fval;
outputPGIss.MSparestReg.exitflag       = exitflag;
outputPGIss.MSparestReg.output         = output;
outputPGIss.MSparestReg.solutions      = solutions2;
  

%% (9) Parameter estimation. Regularized. Identifiability. PLA + mPLA
setup.parEst.nonreg     = 0;
setup.parEst.reg        = 1;
[h,PLAresult] = runPLA(b,x,canelas_SS,setup);
outputPGIss.PLAreg                       = PLAresult;


%% (10) Final save
% save(saveloc,'outputPGIss');
% % % % save('outputPGIss.mat','outputPGIss')


%% visualization of results
setup.parEst.drawnow = 1; %1;

if setup.parEst.drawnow == 1
    % early setup
    load('outputPGIss.mat','outputPGIss');
% % % %     save(saveloc,'outputPGIss');
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
        xres = outputPGIss.costfunCheck.xres;
        [error] = costfunPGI(xres,canelas_SS,setup,x);
% % % %         savefig(1000,[savefigloc, '\PGI_SScosfunNonreg']);
    end
    setup.parEst.drawnowCF      = tempdnCF;
    clear tempdnCF
    
    %     MSparestNonReg: [1×1 struct]
    setup.parEst.MSnonReg = 1;
    setup.parEst.MSReg = 0;
    solutions1 = outputPGIss.MSparestNonReg.solutions;
    [h] = histMS(solutions1,setup);
% % % %     savefig(h,[savefigloc, '\PGI_SSMSnonreg']);
    
    %          PLAnonReg: [1×1 struct]
    [h] = plotPLA(outputPGIss,setup);
% % % %     savefig(h,[savefigloc, '\PGI_PLAnonreg.fig']);

    setup.parEst.MSnonReg = 0;
    setup.parEst.MSReg = 1;
    setup.parEst.nonreg = 0;
    setup.parEst.reg = 1;
    
    %     regularization: [1×1 struct]
    [h] = plotRegularization(outputPGIss,setup);
% % % %     savefig(h,[savefigloc, '\PGI_regularization.fig']);
    
    %        MSparestReg: [1×1 struct]
    solutions2 = outputPGIss.MSparestReg.solutions;
    [h] = histMS(solutions2,setup);
% % % %     savefig(h,[savefigloc, '\PGI_SSMSreg']);

    %             PLAreg: [1×1 struct]
    [h] = plotPLA(outputPGIss,setup);
% % % %     savefig(h,[savefigloc, '\PGI_PLAreg.fig']);
end


%% Network analysis
%     %% (1) Literature Parameters. Dynamics
%     % Sensitivity analysis of the different equation contents over the
%     % simulation.
%         setup.refVals.metabolites(4) = setup.refVals.metabolites(4) * 0.25; % product needs to be decreased for reaction in the right direction
%     [xdata, ydata] = sPSA(setup, canelas_SS);
%     [h] = plotSPSA(xdata,ydata,setup);
%     % savefig(h,'D:\PROJECT\2 Working files\Chapter 3 Consensus model\yeast_consensus_v1\results\figures\PGI_SSstudy\PGI_sPSA.fig')
%         % !!! sPSA --> runPSA: To simplify how the reaction rates are
%         % calculated, create a 3rd dimension that will have 2 options, one for
%         % changing variables and the other for the setpoint. When the initial i
%         % and j are obtained, a k value can be calculated to set when to fit
%         % one, an when the other. If the value is not the one of assay, give it
%         % a 1 (sp matrix) if not, get from dimension 2 (changing matrix).
% 
%     %% (2) Sensitivity in steady conditions. Random parameter sampling. Analysis of parameter effect and dependencies
% 
%     %% (3) Sensitivity in dynamic conditions. Random parameter sampling. Analysis of parameter effect and dependencies
%     % A simulation of the GP is use to study the effect of the parameters.
%     % In here, the model intially received has been used.
%         % !!! Here we are using the inital model provided by joep. Once the entire
%         % system ODE has been updated for this version, change it. 
% 
%     % Simulation
%     % Initializing setup
%     % setup.subsystem.pPGI = - 1 + (1+1) * rand(setup.subsystem.numEvals,4);
%     [T,Y,V,parray] = subsystemDynamics(setup,canelas_SS,vHeerden_GP);
% 
%     % Manual clustering criteria (upon later visualization. Simulations that:
%     [ind] = clusteringSensitivity(T,Y,V,setup);
% 
%     % Plotting
%     plotSubsystem(T,Y,V,legenda,setup,vHeerden_GP,parray,ind)
%     % h1 = figure(1);
%     % h2 = figure(2);
%     % h3 = figure(3);
%     % h4 = figure(4);
%     % savefig(h1,'D:\PROJECT\2 Working files\Chapter 3 Consensus model\yeast_consensus_v1\results\figures\PGI_GPsensitivity_vPGI.fig')
%     % savefig(h2,'D:\PROJECT\2 Working files\Chapter 3 Consensus model\yeast_consensus_v1\results\figures\PGI_GPsensitivity_G6P.fig')
%     % savefig(h3,'D:\PROJECT\2 Working files\Chapter 3 Consensus model\yeast_consensus_v1\results\figures\PGI_GPsensitivity_F6P.fig')
%     % savefig(h4,'D:\PROJECT\2 Working files\Chapter 3 Consensus model\yeast_consensus_v1\results\figures\PGI_GPsensitivity_dependencies.fig')
%     % savefig(1,'D:\PROJECT\2 Working files\Chapter 3 Consensus model\yeast_consensus_v1\results\figures\PGI_SSstudy\PGI_GPsensitivity_dependencies2.fig')
% 
%     % !!!savefig(h,'subsystemDynamicsPGI.fig') % (1) Needs to change the path
%     % where is is saved from C to D drive and (2) group all figures (hq to h4
%     % here) into a single h. All inside.


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


% % % high and double maxima
% % A1 = zeros(setup.subsystem.numEvals,1);
% % for i = 1:setup.subsystem.numEvals
% %     A1(i) = V.rand{i}(2017,3);
% % end
% % iTOP    = find(A1>0.5);
% % setup.subsystem.iTOP = iTOP; clear iTOP;
% % 
% % % slow start
% % A1 = zeros(setup.subsystem.numEvals,1);
% % for i = 1:setup.subsystem.numEvals
% % %     A1(i) = V.rand{i}(2002,3);
% %     A1(i) = V.rand{i}(end,3);
% % end
% % % iSLOW   = find(A1<0.05 & A1>0.025);
% % iSLOW   = find(A1<0.21 & A1>0.1499);
% % setup.subsystem.iSLOW = iSLOW; clear iSLOW;
% % 
% % % very slow start
% % A1 = zeros(setup.subsystem.numEvals,1);
% % for i = 1:setup.subsystem.numEvals
% % %     A1(i) = V.rand{i}(2002,3);
% %     A1(i) = V.rand{i}(end,3);
% % end
% % % iMSLOW   = find(A1<0.025);
% % iMSLOW   = find(A1<0.1499);
% % 
% % setup.subsystem.iMSLOW = iMSLOW; clear iSLOW;

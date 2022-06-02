%% parameter estimation ENO enzyme
    % (0) Initial setup
    % (1) ParEst pipeline
    % (2) Visualization
    % (3) Sensitivity

%% (0) Initial setup
clc, clear, close all
set_paths
% % % %     saveloc = [folder, '\results\variables\divide_and_conquer\xENO_SSotput.mat'];
% % % %     savefigloc  = [folder, '\results\figures\ENO_SSstudy'];
%     saveloc     = [folder, '\results\variables\PDC_ADH_SSstudy\xPDC_ADH_outputReg.mat'];
%     savefigloc  = [folder, '\results\figures\PDC_ADH_SSstudy\PDC_ADH_regularization.fig'];
%     [h,output] = regularization(canelas_SS,setup);
%     save(saveloc,'output')
%     savefig(100,savefigloc)
dbstop if error
CanelasData2;
vHeerdenData;
metadata;
refferenceValues;
rng(1)
setupParestENO; % ENO
load('parameters_blank.mat');


%% Parameter estimation
%% (4) Parameter estimation. Study objective funtion.
% parameter estimation using lsqnonlin
[xres,resnorm,residual,exitflag] = runLSQsimple(canelas_SS,setup);
outputENOss.costfunCheck.xres        = xres;
outputENOss.costfunCheck.resnorm     = resnorm;
outputENOss.costfunCheck.residual    = residual;
outputENOss.costfunCheck.exitflag    = exitflag;

%% testing the output function's results
its = exitflag.history.iteration;

figure

subplot(3,5,1), %[18×4 double]
plot(its, exitflag.history.x,'.-'), title('x')

subplot(3,5,3), %[18×1 double]
plot(its, exitflag.history.iteration,'.-'), title('iteration')
subplot(3,5,4), %[18×1 double]
plot(its, exitflag.history.funccount,'.-'), title('funccount')
subplot(3,5,5), %[18×1 double]
plot(its, exitflag.history.stepsize,'.-'), title('stepsize')
subplot(3,5,6), %[18×4 double]
plot(its, exitflag.history.gradient,'.-'), title('gradient')

subplot(3,5,9), %0
plot(its, exitflag.history.positivedefinite,'.-'), title('positivedefinite')
subplot(3,5,10), %[18×1 double]
plot(its, exitflag.history.ratio,'.-'), title('ratio')
subplot(3,5,11), %[18×1 double]
plot(its, exitflag.history.degenerate,'.-'), title('degenerate')
subplot(3,5,12), %[18×1 double]
plot(its, exitflag.history.trustregionradius,'.-'), title('trustregionradius')
subplot(3,5,13), %[18×12 double]
plot(its, exitflag.history.residual,'.-'), title('residual')
subplot(3,5,14), %[18×1 double]
plot(its, exitflag.history.resnorm,'.-'), title('resnorm')
subplot(3,5,15), %[18×1 double]
plot(its, exitflag.history.cgiterations,'.-'), title('cgiterations')


%% (5) Parameter estimation. Non-regularized. Multi Start.
[b,fval,exitflag,output,solutions1] = MSparest(canelas_SS,setup);
outputENOss.MSparestNonReg.b              = b;
outputENOss.MSparestNonReg.fval           = fval;
outputENOss.MSparestNonReg.exitflag       = exitflag;
outputENOss.MSparestNonReg.output         = output;
outputENOss.MSparestNonReg.solutions      = solutions1;
  

%% (6) Parameter estimation. Non-regularized. Identifiability. PLA.
[h,PLAresult] = runPLA(b,x,canelas_SS,setup);
outputENOss.PLAnonReg                       = PLAresult;


%% (7) Parameter estimation. Regularization.
[h,output] = regularization(canelas_SS,setup);
outputENOss.regularization           = output;


%% (8) Parameter estimation. Regularized. Multi Start.
setup.parEst.MSnonReg         = 0;
setup.parEst.MSReg            = 1;
[b,fval,exitflag,output,solutions2] = MSparest(canelas_SS,setup);
outputENOss.MSparestReg.b              = b;
outputENOss.MSparestReg.fval           = fval;
outputENOss.MSparestReg.exitflag       = exitflag;
outputENOss.MSparestReg.output         = output;
outputENOss.MSparestReg.solutions      = solutions2;
  

%% (9) Parameter estimation. Regularized. Identifiability. PLA + mPLA
setup.parEst.nonreg     = 0;
setup.parEst.reg        = 1;
[h,PLAresult] = runPLA(b,x,canelas_SS,setup);
outputENOss.PLAreg                       = PLAresult;


%% (10) Final save
% save(saveloc,'outputENOss');
% % % % save('outputENOss.mat','outputENOss')


%% visualization of results
setup.parEst.drawnow = 1; %1;

if setup.parEst.drawnow == 1
    % early setup
    load('outputENOss.mat','outputENOss');
%     mkdir(saveloc);
% % % %     save(saveloc,'outputENOss');
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
        xres = outputENOss.costfunCheck.xres;
        [error] = costfunENO(xres,canelas_SS,setup,x);
% % % %         savefig(1000,[savefigloc, '\ENO_SScosfunNonreg']);
    end
    setup.parEst.drawnowCF      = tempdnCF;
    clear tempdnCF
    
    %     MSparestNonReg: [1×1 struct]
    setup.parEst.MSnonReg = 1;
    setup.parEst.MSReg = 0;
    solutions1 = outputENOss.MSparestNonReg.solutions;
    [h] = histMS(solutions1,setup);
% % % %     savefig(h,[savefigloc, '\ENO_SSMSnonreg']);
    
    %          PLAnonReg: [1×1 struct]
    [h] = plotPLA(outputENOss,setup);
% % % %     savefig(h,[savefigloc, '\ENO_PLAnonreg.fig']);

    setup.parEst.MSnonReg = 0;
    setup.parEst.MSReg = 1;
    setup.parEst.nonreg = 0;
    setup.parEst.reg = 1;
    
    %     regularization: [1×1 struct]
    [h] = plotRegularization(outputENOss,setup);
% % % %     savefig(h,[savefigloc, '\ENO_regularization.fig']);
        % % extra plot
        % tempdnCF = setup.parEst.drawnowCF;
        % setup.parEst.drawnowCF      = 1;
        % if setup.parEst.drawnowCF == 1
        %     xres = outputENOss.regularization.xres_array(17,:);
        %     setup.parEst.lambda = setup.parEst.lambda1;
        %     [error] = costfunENO(xres,canelas_SS,setup,x);
        %     % savefig(1000,[savefigloc, '\ENO_SScosfunReg']);
        % end 
        % setup.parEst.drawnowCF      = tempdnCF;
        % clear tempdnCF
    
    %        MSparestReg: [1×1 struct]
    solutions2 = outputENOss.MSparestReg.solutions;
    [h] = histMS(solutions2,setup);
% % % %     savefig(h,[savefigloc, '\ENO_SSMSreg']);

    %             PLAreg: [1×1 struct]
    [h] = plotPLA(outputENOss,setup);
% % % %     savefig(h,[savefigloc, '\ENO_PLAreg.fig']);
    
end



%% Network analysis
%     %% (1) Literature Parameters. Dynamics
%     % Sensitivity analysis of the different equation contents over the
%     % simulation.
%     setup.refVals.metabolites(12) = setup.refVals.metabolites(12) * 0.25; % product needs to be decreased for reaction in the right direction
%     [xdata, ydata] = sPSA(setup, canelas_SS);
%     [h] = plotSPSA(xdata,ydata,setup);
%     % savefig(1,'D:\PROJECT\2 Working files\Chapter 3 Consensus model\yeast_consensus_v1\results\figures\ENO_SSstudy\ENO_sPSA.fig')
% 
%     %% (2) Sensitivity in steady conditions. Random parameter sampling. Analysis of parameter effect and dependencies
% 
% 
%     %% (3) Sensitivity in dynamic conditions. Random parameter sampling. Analysis of parameter effect and dependencies
%     % A simulation of the GP is use to study the effect of the parameters.
%     % In here, the model intially received has been used.
%         % !!! Here we are using the initla model provided by joep. Once the entire
%         % system ODE has been updated for this version, change it. 
% 
%     % Simulation
%     [T,Y,V,parray] = subsystemDynamics(setup,canelas_SS,vHeerden_GP);
% 
%     % Manual clustering criteria (upon later visualization. Simulations that:
%     [ind] = clusteringSensitivity(T,Y,V,setup);
% 
%     % Plotting
%     plotSubsystem(T,Y,V,legenda,setup,vHeerden_GP,parray,ind)
%     % savefig(4,'D:\PROJECT\2 Working files\Chapter 3 Consensus model\yeast_consensus_v1\results\figures\ENO_SSstudy\ENO_dependencies.fig')






%% memoryDump

% setup.parEst.drawnow = setup.parEst.drawnow; %
% 
% if setup.parEst.drawnow == 1
%     %       costfunCheck: [1×1 struct]
%     %     MSparestNonReg: [1×1 struct]
%     %          PLAnonReg: [1×1 struct]
%     %     regularization: [1×1 struct]
%     %        MSparestReg: [1×1 struct]
%     %             PLAreg: [1×1 struct]
% end



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

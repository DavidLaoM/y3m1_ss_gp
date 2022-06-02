function [simResults] = simulateY3M1_separateGPSS_safecopy(xAll, canelas_SS, data, dataset, setup)
% simulation of teh experimental datasets with a given parameter array (if
% npSets = 1) or a series of parameter arrays (if npSets > 1).

[npSets,~] = size(xAll);
simResults = cell(npSets,1);
% disp message
if setup.plotResultsMode == 30
    f = waitbar(0,'starting simulations...');
    nbytes = fprintf('Simulation 0 of %d', npSets);
    s = clock;
end
if setup.plotResultsMode == 81
    nbytes = fprintf('Simulation 0 of %d', npSets);
    s = clock;
end

for i = 1:npSets
    x = xAll(i,:);
    % (2021/11/07) Added option to test if all parameters can be unified
    if(isfield(setup,'sameGLTHXK4all')&&(setup.sameGLTHXK4all == 1))
        x(144:146) = x([35 36 38]); % GLT fixing
        x(137:143) = x(28:34); % HXK fixing
    end
    % simulations
    if setup.runSScanelas  == 1
        [Tss,Yss,Vss] = simulateSScanelas_Y3M1(x,canelas_SS,data,setup);
    end
    if setup.runGSvanHeerden == 1
        x2 = x;
        x2(11) = x2(147);     %     x(147) = x(11);     % fba  (11:15)
        x2(52) = x2(148);     %     x(148) = x(52);     % pfk  (43:56)
        x2(59) = x2(149);     %     x(149) = x(59);     % pgi  (57:60)
        x2(79) = x2(150);     %     x(150) = x(79);     % tpi  (79:82)
        x2(22) = x2(151);     %     x(151) = x(22);     % gapdh(21:27)
        [Tgs,Ygs,Vgs] = simulateGSvanHeerden_Y3M1(x2,canelas_SS,data,setup);
    end
    % save results
    if setup.runSScanelas  == 1
        simResult.ss.T = Tss;
        simResult.ss.Y = Yss;
        simResult.ss.V = Vss;
    end
    if setup.runGSvanHeerden == 1
        simResult.gs.T = Tgs;
        simResult.gs.Y = Ygs;
        simResult.gs.V = Vgs;
    end
    simResults{i,1} = simResult;
%     fprintf('Iteration through p.set %4.2f done \n',i);
    % disp message
    if setup.plotResultsMode == 30
%         fprintf(repmat('\b',1,nbytes))
%         nbytes = fprintf('Simulation %d of %d\n', i, npSets);
        if i == 1
            is = etime(clock,s);
            esttime = is * npSets;
        end
        fprintf(repmat('\b',1,nbytes))
        remTime = ['remaining time = ',num2str((esttime-etime(clock,s)),'%4.1f'),'seconds'];
        nbytes = fprintf(['Simulation %d of %d\n',remTime], i, npSets);
        %waitbar
        waitbarText1 = sprintf('Simulation %d of %d',i,npSets);
        waitbarText2 = ['remaining time = ',num2str((esttime-etime(clock,s)),'%4.1f'),'seconds'];
        waitbarText = [waitbarText1,' ', waitbarText2];
        waitbar(i/npSets,f,waitbarText);      
    end
    if setup.plotResultsMode == 81
        if i == 1
            is = etime(clock,s);
            esttime = is * npSets;
        end
        fprintf(repmat('\b',1,nbytes))
        remTime = ['remaining time = ',num2str((esttime-etime(clock,s)),'%4.1f'),'seconds'];
        nbytes = fprintf(['Simulation %d of %d\n',remTime], i, npSets);
    end
    
end
% plotting results
if setup.plotResultsMode ~= 0
    [legenda] = legendaFull;
    [~] = plotAll_Y3M1(simResults,legenda,canelas_SS,data,setup,xAll);
end

end
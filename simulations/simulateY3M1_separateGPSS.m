function [simResults] = simulateY3M1_separateGPSS(xAll, canelas_SS, data, dataset, setup)
% simulation of teh experimental datasets with a given parameter array (if
% npSets = 1) or a series of parameter arrays (if npSets > 1).

[npSets,~] = size(xAll);
simResults = cell(npSets,1);
for i = 1:npSets
    % recall paramter set
    x = xAll(i,:);
    % simulations
    if setup.runSScanelas  == 1
        [Tss,Yss,Vss] = simulateSScanelas_Y3M1(x,canelas_SS,data,setup);
    end
    if setup.runGSvanHeerden == 1
        x2 = x; % slowed upper glycolysis upon glucose perturbation
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
end

% plotting results
if setup.plotResultsMode ~= 0
    [legenda] = legendaFull;
    [~] = plotAll_Y3M1(simResults,legenda,canelas_SS,data,setup,xAll);
end

end
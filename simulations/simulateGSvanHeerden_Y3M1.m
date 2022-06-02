function [T,Y,V]=simulateGSvanHeerden_Y3M1(x,canelas_SS,data,setup)
% simulations of the van Heerden dataset. Glucose pulse of 110 mM at
% dilution rate of 0.1 h^{-1}. BioScope bioreactor.

% for GP simulations
setup.conditionsSS = 0;
setup.conditionsGP = 1;

% simulation 1: fixed growth rate value (0.1 h-1)
% setup 
experiment = 3;
d = canelas_SS.mtlD.D(experiment);
PvHoek_EnzymeExpressionData;
setParameterStructure_Y3M1; 
InitCond_SS_0_1_Y3M1;
IC(9) = 4.43./2; %ATP % AXP initial concentrations in (van Heerden et al., 2014).
IC(15) = 1.06./2; %ADP
IC(16) = 0.0693./2; %AMP
tspan=[0:1:3000]; 
f.GLCo=canelas_SS.mtlD.Cs(experiment); 
options=odeset('RelTol',1e-4,'RelTol',1e-4 ,'NonNegative',ones(1,length(IC)));
if setup.clamp10.IXP == 1 % if IXP is clamped, model structure is adjusted from here on
    setup.conditionsSS = 1;
    setup.conditionsGP = 0;
end
% actual simulation
p.ATPase_ratio2 = 10.^x(87); % ATPase activity is different before and during the 110 mM glucose perturbation
setup.sim1_ixp_off = 0; % IXP fluxes are kept off before the 110 mM glucose perturbation
[T_ic2ss,Y_ic2ss] = ode15s(@ODE_model_Y3M1,tspan,IC,options, p,f,d,canelas_SS,data,experiment,setup);
setup.sim1_ixp_off = 1;
% structure output
Y = Y_ic2ss;
calcFluxes_consensus_Y3M1;
V_ic2ss = v;
tspan = [0:1:340];

% simulation 2: extracellular glucose 110 mM concentration perturbation
f.GLCo = 110; % increase extracellular glucose concnetration (step, fixed) upon 110 mM glucose perturbation
Y_ic2ss(end,27) = 25; % Increased cytosolic inorganic phosphate upon 110 mMglucose perturbation
setup.increase_Katpase = 1; % ATPase activity is increased
setParameterStructure_Y3M1; 
[T_gs,Y_gs] = ode15s(@ODE_model_Y3M1,tspan,Y_ic2ss(end,:),options, p,f,d,canelas_SS,data,experiment,setup);
    
% pre-structure output
if setup.clamp10.TRE == 1 % in case of clamped parts of the model, this is needed for calculating reaction rates
    Y_gs(:,21) = interp1(data.time_metabolites,data.metabolites(:,13),T_gs,'pchip','extrap');
    Y_gs(:,24) = interp1(data.time_metabolites,data.metabolites(:,14),T_gs,'pchip','extrap');
    Y_gs(:,26) = interp1(data.time_metabolites,data.metabolites(:,15),T_gs,'pchip','extrap');
    Y_gs(:,25) = interp1(data.time_metabolites,data.metabolites(:,16),T_gs,'pchip','extrap');
end
if setup.clamp10.AXP == 1
    Y_gs(:,9) = interp1(data.time_nucleotides,data.nucleotides(:,4),T_gs,'pchip','extrap') + 1.0810;
    Y_gs(:,15) = interp1(data.time_nucleotides,data.nucleotides(:,3),T_gs,'pchip','extrap') + 0.3613;
    Y_gs(:,16) = interp1(data.time_nucleotides,data.nucleotides(:,2),T_gs,'pchip','extrap') + 0.2238;
end
if setup.clamp10.Pi == 1 % if clamped, no need for the initial increase
    Y_gs(:,27) = 10*ones(size(Y_gs(:,27)));
end
T = T_gs;
Y = Y_gs;
calcFluxes_consensus_Y3M1;
V_gs = v;

% structure output
T = [T_ic2ss-3000; T_gs];
Y = [Y_ic2ss; Y_gs];
V = [V_ic2ss; V_gs];
end

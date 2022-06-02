function [T,Y,V] = simulateSScanelas_Y3M1(x,canelas_SS,data,setup)
% simulations of the canelas dataset. Steady state concentrations at growth
% rates from 0.020 to 0.375 h^{-1}.
setup.conditionsSS = 0;
setup.conditionsGP = 0;

% simulation 1: fixed growth rate value (0.1 h-1)
% setup
experiment = 3;
d = canelas_SS.mtlD.D(experiment);
PvHoek_EnzymeExpressionData;
setParameterStructure_Y3M1;
InitCond_SS_0_1_Y3M1;
tspan=[0:1:3000];    
f.GLCo=canelas_SS.mtlD.Cs(experiment); 
options=odeset('RelTol',1e-4,'RelTol',1e-4 ,'NonNegative',ones(1,length(IC)));
% actual simulation
[T_ic2ss,Y_ic2ss] = ode15s(@ODE_model_Y3M1,tspan,IC,options,p,f,d,canelas_SS,data,experiment,setup);
% structure output
Y = Y_ic2ss;
calcFluxes_consensus_Y3M1;
V_ic2ss = v;
[~,ny] = size(Y_ic2ss);
[~,nv] = size(V_ic2ss);
Yt = zeros(length(tspan), ny);
Vt = zeros(length(tspan), nv);
Ytt = cell(length(canelas_SS.mtlD.D),1); %Yc{:} = Yce;
Vtt = cell(length(canelas_SS.mtlD.D),1); %Vc{:} = Vce; 
for i = 1:length(canelas_SS.mtlD.D)
    Ytt{i} = Yt;
    Vtt{i} = Vt;
end

% simulation 2: growth rate dependent (0.025 - 0.375 h-1)
for experiment = 1:length(canelas_SS.mtlD.D)
    
    % adjusting growth rate dependency
    % HXT growth rate dependency
    if experiment == 1, safecopy_x36 = x(36); end
    if experiment == 1, x(36) = -0.05; %x(36) = safecopy_x36-1-0.1-0.1; 
    elseif experiment == 2, x(36) = 0.8; %x(36) = safecopy_x36-0.1-0.2-0.2;%-0.2; 
    elseif experiment == 3, x(36) = safecopy_x36; %x(36) = safecopy_x36-0.1-0.2;%-0.1;
    elseif experiment == 4, x(36) = safecopy_x36; %x(36) = safecopy_x36;
    elseif experiment == 5, x(36) = 1.15; %x(36) = safecopy_x36+1;
    elseif experiment == 6, x(36) = 1.5; %x(36) = safecopy_x36+1;
    elseif experiment == 7, x(36) = 1.26; %x(36) = safecopy_x36+1; %+2; %1
    elseif experiment == 8, x(36) = 1.2; %x(36) = safecopy_x36+1; %+2; %1
    else, disp('something went wrong, no -experiment- being selected')
    end
    setParameterStructure_Y3M1; 
    % possibility to consider PFK growth rate dependency
    updateParameters_PFK_Y3M1;
    % 
    Yss2 = Y_ic2ss;
    d=canelas_SS.mtlD.D(experiment);
    PvHoek_EnzymeExpressionData;
    if(isfield(setup,'vanHoek_off')&&(setup.vanHoek_off == 1)) % option to shut down growth rate dependent glycolytic protein abundance
        PvHoek_EnzymeExpressionData_dummy;
    end
    f.GLCo=canelas_SS.mtlD.Cs(experiment);
    
    % option to clamp some modules in the model
    % (when simulating steady states, IXP and Trecycle are kept clamped)
    %PI
    if setup.clamp10.Pi == 1
        Yss2(end,27)    = 10;
    end
    %NADX
    if setup.clamp10.NADX == 1
        Yss2(end,7)     = canelas_SS.mtlD.NAD(experiment);     %NAD
        Yss2(end,8)     = canelas_SS.mtlD.NADH(experiment);    %NADH
    end
    % AXP
    if setup.clamp10.AXP == 1
        Yss2(end,9)     = canelas_SS.mtlD.ATP(experiment);     %ATP
        Yss2(end,15)    = canelas_SS.mtlD.ADP(experiment);     %ADP
        Yss2(end,16)    = canelas_SS.mtlD.AMP(experiment);     %AMP
    elseif((setup.clamp10.AXP == 0))%&&(setup.clamp10.IXP == 1))
        Yss2(end,9)     = canelas_SS.mtlD.ATP(experiment);     %ATP
        Yss2(end,15)    = canelas_SS.mtlD.ADP(experiment);     %ADP
        Yss2(end,16)    = canelas_SS.mtlD.AMP(experiment);     %AMP
    end
    %TRE
    Yss2(end,26)    = canelas_SS.mtlD.T6P(experiment);    %T6P
    Yss2(end,21)    = canelas_SS.mtlD.G1P(experiment);    %G1P
    Yss2(end,25)    = canelas_SS.mtlD.TRE(experiment);    %TRE
    %IXP
    Yss2(end,28)    = 0.1;%IMP
    Yss2(end,29)    = 0.1;%INO
    Yss2(end,30)    = 1.5;%HYP
    
    % actual simulation
    [T_ss2can,Y_ss2can] = ode15s(@ODE_model_Y3M1,tspan,Yss2(end,:),options,p,f,d,canelas_SS,data,experiment,setup);

    % Pre-structure Output
    Y = Y_ss2can;
    calcFluxes_consensus_Y3M1;
    V_ss2can = v;
    Ytt{experiment} = Y_ss2can;
    Vtt{experiment} = V_ss2can;  
end

% Structure Output
T.ic2ss     = T_ic2ss;
T.ss2can    = T_ss2can;
    clear Y
Y.ic2ss     = Y_ic2ss;
Y.ss2can    = Ytt;
V.ic2ss     = V_ic2ss;
V.ss2can    = Vtt;
end

function v = ODE_model_Y3M1(t,IC,p,f,d,canelas_SS,data,experiment,setup)
% system of ODE equations in the Y3M1 model. The clamps affect the
% calculation of reaction rates in rateEquations_Y3M1.
% Since some reactions in glycolysis do also affect the cofactors mass
% balance, there is an additional effect of the clamps at the end of this
% file.

% recall metabolite concentrations
ACE=IC(1);
BPG=IC(2);
F16BP=IC(3);
F6P=IC(4);
G6P=IC(5);
GLCi=IC(6);
NAD=IC(7);
NADH=IC(8);
ATP=IC(9); 
P2G=IC(10);
P3G=IC(11);
PEP=IC(12);
PYR=IC(13);
GLYCERAL3P=IC(14);
ADP=IC(15);
AMP=IC(16); 
DHAP=IC(17);
GLYC3P=IC(18);
GLYCEROL=IC(19);
ETOH=IC(20);
G1P=IC(21);
UTP=IC(22);
UDP=IC(23);
UDP_GLC=IC(24);
TRE=IC(25);
T6P=IC(26);
PI=IC(27);
IMP=IC(28);
INO=IC(29);
HYP=IC(30);
ETOHec=IC(31);
GLYCec=IC(32);

% clamps, concentrations
if((setup.conditionsSS == 0)&&(setup.conditionsGP == 1))% second simulation GP
    if setup.clamp10.TRE == 1 % trehalose cycle
        UDP_GLC=interp1(data.time_metabolites,data.metabolites(:,14),t,'pchip','extrap');
        TRE=interp1(data.time_metabolites,data.metabolites(:,16),t,'pchip','extrap');
        T6P=interp1(data.time_metabolites,data.metabolites(:,15),t,'pchip','extrap');
        G1P=interp1(data.time_metabolites,data.metabolites(:,13),t,'pchip','extrap');
    end
    if setup.clamp10.AXP == 1 % axp
        ATP = interp1(data.time_nucleotides,data.nucleotides(:,4),t,'linear','extrap');
        ADP = interp1(data.time_nucleotides,data.nucleotides(:,3),t,'linear','extrap');
        AMP = interp1(data.time_nucleotides,data.nucleotides(:,2),t,'linear','extrap');
    end
    if setup.clamp10.Pi == 1 % pi
        PI     = 10; %PI
    end
end

% reaction rates
rateEquations_Y3M1;

% clamps, reaction rates
if setup.clamp10.Pi == 1
    v_vacuolePi = zeros(size(v_vacuolePi));
end
if setup.clamp10.NADX == 1
    v_mitoNADH = zeros(size(v_mitoNADH));
end    
if setup.clamp10.AXP == 1
    v_ADK1 = zeros(size(v_ADK1));
    v_ATPase = zeros(size(v_ATPase));
    v_mito = zeros(size(v_mito));
end  
% clamps, trehalose and inosine salvage pathways are now modelled during steady state simulations.
if(((setup.conditionsSS == 0)&&(setup.conditionsGP == 0))||((setup.conditionsSS == 1)&&(setup.conditionsGP == 0)))
    v_PGM1 = zeros(size(v_PGM1));
    v_TPS1 = zeros(size(v_TPS1));
    v_TPS2 = zeros(size(v_TPS2));
    v_NTH1 = zeros(size(v_NTH1));
    v_Amd1 = zeros(size(v_Amd1));
    v_Ade13_v_Ade12 = zeros(size(v_Ade13_v_Ade12));
    v_Isn1 = zeros(size(v_Isn1));
    v_Pnp1 = zeros(size(v_Pnp1));
    v_Hpt1 = zeros(size(v_Hpt1));
    dAXPdt = 0; % option used meanwhile developing AXP and IXP, not anymore.
    dAXPdD = 0; % option used meanwhile developing AXP and IXP, not anymore.
elseif((setup.conditionsSS == 0)&&(setup.conditionsGP == 1))% second simulation GP
    if setup.clamp10.TRE == 1
        v_PGM1 = interp1(data.time_fluxes,(data.fluxes(:,6)-data.fluxes(:,7)),t,'pchip','extrap');
        v_TPS1 = interp1(data.time_fluxes,data.fluxes(:,9),t,'pchip','extrap');
        v_TPS2 = interp1(data.time_fluxes,data.fluxes(:,10),t,'pchip','extrap');
        v_NTH1 = interp1(data.time_fluxes,data.fluxes(:,11),t,'pchip','extrap');
    elseif setup.clamp10.TRE == 2
        v_PGM1 = zeros(size(v_PGM1));
        v_TPS1 = zeros(size(v_TPS1));
        v_TPS2 = zeros(size(v_TPS2));
        v_NTH1 = zeros(size(v_NTH1));
    end
    dAXPdt = 0; % option used meanwhile developing AXP and IXP, not anymore.
    dAXPdD = 0; % option used meanwhile developing AXP and IXP, not anymore.
    if setup.clamp10.IXP == 1
        v_Amd1 = zeros(size(v_Amd1));
        v_Ade13_v_Ade12 = zeros(size(v_Ade13_v_Ade12));
        v_Isn1 = zeros(size(v_Isn1));
        v_Pnp1 = zeros(size(v_Pnp1));
        v_Hpt1 = zeros(size(v_Hpt1));
        dAXPdt = interp1(data.time_dAXPdt,data.dAXPdt,t,'pchip','extrap');
    end
end
if(isfield(setup,'sim1_ixp_off')&&(setup.sim1_ixp_off == 0)) % IXP fluxes are kept off before the 110 mM glucose perturbation
	v_Amd1 = zeros(size(v_Amd1));
    v_Ade13_v_Ade12 = zeros(size(v_Ade13_v_Ade12));
    v_Isn1 = zeros(size(v_Isn1));
	v_Pnp1 = zeros(size(v_Pnp1));
    v_Hpt1 = zeros(size(v_Hpt1));
end

% mass balances, system of ODEs
v(1)= + v_PDC - v_ADH + v_sinkACE;% Using sinkACE (- v_ACE)
v(2)= + v_GAPDH - v_PGK;%BPG
v(3)= +v_PFK - v_ALD;%F16BP
v(4)= -v_PFK + v_PGI + v_sinkF6P;%F6P
v(5)= +v_GLK - v_PGI + v_sinkG6P;%G6P
v(6)= -v_GLK + v_GLT;%GLCi
v(7)= + v_G3PDH - v_GAPDH + v_ADH + v_mitoNADH;%NAD
v(8)= - v_G3PDH + v_GAPDH - v_ADH - v_mitoNADH;%NADH
v(9)=  +v_ADK1 -v_GLK - v_ATPase - v_PFK +v_PGK + v_PYK - v_TPS1 + v_mito + dAXPdt.*(ATP./(ATP+ADP+AMP)) + dAXPdD.*(ATP./(ATP+ADP+AMP)); %ATP
v(10)= +v_PGM - v_ENO; %P2G
v(11)= +v_PGK - v_PGM + v_sinkP3G; %P3G
v(12)= +v_ENO - v_PYK + v_sinkPEP; % PEP
v(13)= +v_PYK - v_PDC + v_sinkPYR; % PYR
v(14)= + v_ALD - v_GAPDH + v_TPI1 + v_sinkGAP; %GLYCERAL3P
v(15)=  -2.*v_ADK1 +v_GLK + v_ATPase + v_PFK - v_PGK - v_PYK + v_TPS1 - v_mito + dAXPdt.*(ADP./(ATP+ADP+AMP)) + dAXPdD.*(ADP./(ATP+ADP+AMP));% ADP
v(16)=  +v_ADK1 - v_Amd1 + v_Ade13_v_Ade12 + dAXPdt.*(AMP./(ATP+ADP+AMP)) + dAXPdD.*(AMP./(ATP+ADP+AMP));% AMP
v(17)= + v_ALD - v_TPI1 - v_G3PDH; % DHAP
v(18)=  + v_G3PDH - v_HOR2 - v_RHR2; % GLYC3P
v(19)= + v_HOR2 + v_RHR2 - v_GLYCEROLtransport; %GLYCEROL
v(20)= + v_ADH - v_ETOHtransport; %ETOH
v(21)= - v_PGM1 - v_TPS1; %G1P
v(22) = 0; %UTP
v(23) = 0; %UDP
v(24) = 0; %UDP_GLC
v(25) = + v_TPS2 - v_NTH1; %TRE 
v(26) = + v_TPS1 - v_TPS2; %T6P 
v(27) = - v_GAPDH + v_ATPase + v_HOR2 + v_RHR2 +2.* v_TPS1 +1.* v_TPS2 - v_mito + v_Hpt1- v_Isn1 + v_vacuolePi... %PI
        - v_sinkG6P - v_sinkF6P - v_sinkGAP - v_sinkP3G - v_sinkPEP; %PI (added effect of sinks)
v(28) = +v_Amd1-v_Ade13_v_Ade12+v_Hpt1-v_Isn1; %IMP
v(29) = v_Isn1-v_Pnp1; %INO
v(30) = +v_Pnp1-v_Hpt1; %HYP
v(31) = v_ETOHtransport * canelas_SS.mtlD.ECbiom(experiment)*0.002 - ETOHec/3600*canelas_SS.mtlD.D(experiment); % Product Ethanol
v(32) = v_GLYCEROLtransport * canelas_SS.mtlD.ECbiom(experiment)*0.002 - GLYCec/3600*canelas_SS.mtlD.D(experiment); % Byproduct Glycerol

% trehalose cycle active during 110 mM glucose perturbation
if((setup.conditionsSS == 0)&&(setup.conditionsGP == 1)) 
    v(5)= +v_GLK - v_PGI + v_sinkG6P + v_PGM1 -v_TPS1; %G6P
    v(6)= -v_GLK + v_GLT + 2.*v_NTH1; %GLCi    
end 

% clamps
if setup.clamp10.Pi == 1
    v(27) = 0; %Pi
end
if setup.clamp10.AXP == 1
    v(9)    = 0; %ATP
    v(15)   = 0; %ADP
    v(16)   = 0; %AMP
end    
if(((setup.conditionsSS == 0)&&(setup.conditionsGP == 0))||((setup.conditionsSS == 1)&&(setup.conditionsGP == 0)))% setup for canelas_SS simulations
    % Tre fixed for SS simulations
    v(21)   = 0; %G1P
    v(25)   = 0; %TRE 
    v(26)   = 0; %T6P
    % IXP fixes for GP simulations
    v(28)   = 0; %IMP
    v(29)   = 0; %INO
    v(30)   = 0; %HYP 
    if setup.clamp10.NADX == 1 % NADX
        v(7) = 0; %NAD
        v(8) = 0; %NADH
    end   
elseif((setup.conditionsSS == 0)&&(setup.conditionsGP == 1))% second simulation GP
    if setup.clamp10.TRE == 1 % T6P + G1P
        v(26)   = 0; %T6P
        v(25)   = 0; %TRE
        v(21)   = 0; %G1P
    end
    if setup.clamp10.IXP == 1 % IXP
        v(28)   = 0; %IMP
        v(29)   = 0; %INO
        v(30)   = 0; %HYP  
    end
end

% structure output
v=v';
end
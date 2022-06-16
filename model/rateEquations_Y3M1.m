% All reactions are calculated. The effect of the clamps is added outside 
% this code.

% accounting for growth rate dependency here
if setup.biomass == 1
    % the difference in the sink of pyruvate is understood as a proxy of
    % the effect of the mitochondrial carbon uptake (towards the TCA
    % cycle). This 'relative activity' is used to correct the other
    % cofactors.
    km_sinkPYR      = 0.001; 
    dexp = canelas_SS.mtlD.D(experiment);
    dref = 0.1;
    polyexp = - 8.4853e+03 * dexp.^6 + 9.4027e+03 * dexp.^5 - 3.8027e+03 * dexp.^4 + 700.5 * dexp.^3 - 60.26 * dexp.^2 + 0.711 * dexp - 0.0356; % % INITIAL FIT
    polyref = - 8.4853e+03 * dref.^6 + 9.4027e+03 * dref.^5 - 3.8027e+03 * dref.^4 + 700.5 * dref.^3 - 60.26 * dref.^2 + 0.711 * dref - 0.0356; % % INITIAL FIT
    valexp = -polyexp * (PYR./(PYR + km_sinkPYR));
    valref = -polyref * (PYR./(PYR + km_sinkPYR));

    % ratioATPmito
    ratioATP = 10.^(setup.x133.*(valexp-valref));
    p.mitoVmax = setup.ratioATP.*ratioATP;%100;

    % ratioNADX
    ratioNADH = 10.^(setup.x132.*(valexp-valref));
    p.mitoNADHVmax = setup.ratioNADH.*ratioNADH;%100;

    % rationATPase
    ratioATPase = 10.^(setup.x134.*(valexp-valref));
    p.ATPaseK = setup.ratioATPase.*ratioATPase;%100;

    % ratioADK1_k
    ratioADK1_k = 10.^(setup.x135.*(valexp-valref));
    p.ADK1_k = setup.ratioADK1_k.*ratioADK1_k;

    % ratioADK1_Keq
    ratioADK1_Keq = 10.^(setup.x136.*(valexp-valref));
    p.ADK1_Keq = setup.ratioADK1_Keq.*ratioADK1_Keq;

end


%% UPPER GLYCOLYSIS

%% (vGLT origin?)
% kinetics
v_GLT=p.GLT.VmGLT.*(f.GLCo-GLCi./p.GLT.KeqGLT)./...
    (p.GLT.KmGLTGLCo.*(1+f.GLCo/p.GLT.KmGLTGLCo+GLCi / p.GLT.KmGLTGLCi + ...
    0.91.*f.GLCo.*GLCi./(p.GLT.KmGLTGLCi.*p.GLT.KmGLTGLCo)));

%% HXK1 consensus
% ATP + D-glucose = ADP + "D-glucose 6-phosphate" + H+ 
% % HXK_pH = (-41.76886.*pH.^6 + 1733.52327.*pH.^5 - 29934.91149.*pH.^4 + 275321.55826.*pH.^3 - 1422599.46274.*pH.^2 + 3915935.67919.*pH - 4486910.34625)./100;
v_GLK=p.HXK_ExprsCor.*(((p.HXK1_kcat.*(f.HXK1+f.HXK2))./(p.HXK1_Katp.*p.HXK1_Kglc).*(ATP.*GLCi-((ADP.*G6P)./p.HXK1_Keq)))./...
    ((1+ATP./p.HXK1_Katp+ADP./p.HXK1_Kadp).*(1+GLCi./p.HXK1_Kglc+G6P./p.HXK1_Kg6p+T6P./p.HXK1_Kt6p)));

%% PGI1 consensus
% D-glucose 6-phosphate = D-fructose 6-phosphate 
v_PGI=p.PGI_ExprsCor.*((((p.PGI1_kcat*f.PGI1)./p.PGI1_Kg6p).*(G6P-(F6P./p.PGI1_Keq)))./...
    (1+G6P./p.PGI1_Kg6p+1+F6P./p.PGI1_Kf6p-1));

%% PFK consensus

%% added missing pfk regulation
% used in supplementary materials figure
if ((isfield(setup,'missing_regulation_PFK_option'))&&(setup.missing_regulation_PFK_option == 1)&&(isfield(p,'PFK_factor1')))
    % selecting the factor (deppending on experiment value)
    if setup.experiment2 == 1, p.PFK_factor = p.PFK_factor1; 
    elseif setup.experiment2 == 2, p.PFK_factor = p.PFK_factor2; 
    elseif setup.experiment2 == 3, p.PFK_factor = p.PFK_factor3; 
    elseif setup.experiment2 == 4, p.PFK_factor = p.PFK_factor4; 
    elseif setup.experiment2 == 5, p.PFK_factor = p.PFK_factor5; 
    elseif setup.experiment2 == 6, p.PFK_factor = p.PFK_factor6; 
    elseif setup.experiment2 == 7, p.PFK_factor = p.PFK_factor7; 
    elseif setup.experiment2 == 8, p.PFK_factor = p.PFK_factor8; 
    end
    % specific parameter to change
    if setup.missing_regulation_PFK_parameterNum == 43
        p.PFK_Camp = p.PFK_Camp * p.PFK_factor;
    elseif setup.missing_regulation_PFK_parameterNum == 44
        p.PFK_Catp = p.PFK_Catp * p.PFK_factor;
    elseif setup.missing_regulation_PFK_parameterNum == 45
        p.PFK_Cf16bp = p.PFK_Cf16bp * p.PFK_factor;
    elseif setup.missing_regulation_PFK_parameterNum == 46
        p.PFK_Cf26bp = p.PFK_Cf26bp * p.PFK_factor;
    elseif setup.missing_regulation_PFK_parameterNum == 47
        p.PFK_Ciatp = p.PFK_Ciatp * p.PFK_factor;
    elseif setup.missing_regulation_PFK_parameterNum == 48
        p.PFK_Kamp = p.PFK_Kamp * p.PFK_factor;
    elseif setup.missing_regulation_PFK_parameterNum == 49
        p.PFK_Katp = p.PFK_Katp * p.PFK_factor;
    elseif setup.missing_regulation_PFK_parameterNum == 50
        p.PFK_Kf16bp = p.PFK_Kf16bp * p.PFK_factor;
    elseif setup.missing_regulation_PFK_parameterNum == 51
        p.PFK_Kf26bp = p.PFK_Kf26bp * p.PFK_factor;
    elseif setup.missing_regulation_PFK_parameterNum == 52
        p.PFK_Kf6p = p.PFK_Kf6p * p.PFK_factor;
    elseif setup.missing_regulation_PFK_parameterNum == 53
        p.PFK_Kiatp = p.PFK_Kiatp * p.PFK_factor;
    elseif setup.missing_regulation_PFK_parameterNum == 54
        p.PFK_L = p.PFK_L * p.PFK_factor;
    elseif setup.missing_regulation_PFK_parameterNum == 55
        p.PFK_gR = p.PFK_gR * p.PFK_factor;
    elseif setup.missing_regulation_PFK_parameterNum == 56
        p.PFK_kcat = p.PFK_kcat * p.PFK_factor;
    elseif setup.missing_regulation_PFK_parameterNum == 87
        p.PFK.F26BP = p.PFK.F26BP * p.PFK_factor;
    end    
end
% 
F26BP=p.PFK.F26BP;
PFK_nom=(p.PFK_kcat.*f.PFK.*p.PFK_gR.*(F6P./p.PFK_Kf6p).*(ATP./p.PFK_Katp).*(1+(F6P./p.PFK_Kf6p)+(ATP./p.PFK_Katp)+p.PFK_gR.*((F6P./p.PFK_Kf6p).*(ATP./p.PFK_Katp))));
PFK_denom=(1+F6P./p.PFK_Kf6p+ATP./p.PFK_Katp+(p.PFK_gR.*(F6P./p.PFK_Kf6p).*(ATP./p.PFK_Katp))).^2+...
    p.PFK_L.*...
    ((1+p.PFK_Ciatp.*(ATP./p.PFK_Kiatp))./(1+ATP./p.PFK_Kiatp)).^2.*...
    ((1+p.PFK_Camp.*(AMP./p.PFK_Kamp))./(1+AMP./p.PFK_Kamp)).^2.*...
    ((1+((p.PFK_Cf26bp*F26BP)./(p.PFK_Kf26bp))+((p.PFK_Cf16bp.*F16BP)./(p.PFK_Kf16bp)))./(1+(F26BP./p.PFK_Kf26bp)+(F16BP./p.PFK_Kf16bp))).^2.*...
    (1+p.PFK_Catp.*(ATP./p.PFK_Katp)).^2;
v_PFK=p.PFK_ExprsCor.*(PFK_nom./PFK_denom);

%% FBA1/ALD consensus
% "D-fructose 1,6-bisphosphate" = "glyceraldehyde 3-phosphate" + "glycerone phosphate";
v_ALD=p.FBA_ExprsCor.*(((p.FBA1_kcat.*f.FBA1)./p.FBA1_Kf16bp.*(F16BP-(GLYCERAL3P.*DHAP)./p.FBA1_Keq))./...
    (1+F16BP./p.FBA1_Kf16bp+(1+GLYCERAL3P./p.FBA1_Kglyceral3p).*(1+DHAP./p.FBA1_Kdhap)-1));

%% TPI1 consensus
% "glycerone phosphate" = "glyceraldehyde 3-phosphate";  "triose-phosphate isomerase 
v_TPI1=(((p.TPI1_kcat*f.TPI1)./p.TPI1_Kdhap.*(DHAP-GLYCERAL3P./p.TPI1_Keq))./...
    (1+DHAP./p.TPI1_Kdhap+1+GLYCERAL3P./p.TPI1_Kglyceral3p-1));

%% LOWER GLYCOLYSIS

%% TDH1/GAPDH consensus
% "glyceraldehyde 3-phosphate" + NAD+ + phosphate = "3-phospho-D-glyceroyl dihydrogen phosphate" + H+ + NADH;
GLYCERATE13BP=BPG;
v_GAPDH=p.GAPDH_ExprsCor.*((((p.TDH1_kcat.*(f.TDH1+f.TDH2+f.TDH3))./(p.TDH1_Kglyceral3p.*p.TDH1_Knad.*p.TDH1_Kpi)).*(GLYCERAL3P.*NAD.*PI-(GLYCERATE13BP.*NADH)./p.TDH1_Keq))./...
    ((1+GLYCERAL3P./p.TDH1_Kglyceral3p).*(1+NAD./p.TDH1_Knad).*(1+PI./p.TDH1_Kpi)+(1+GLYCERATE13BP./p.TDH1_Kglycerate13bp).*(1+NADH./p.TDH1_Knadh)-1));

%% PGK teusink
% "3-phospho-D-glyceroyl dihydrogen phosphate" + ADP = "3-phospho-D-glyceric acid" + ATP;  
v_PGK=p.PGK_ExprsCor.*p.PGK.VmPGK.*((p.PGK.KeqPGK.*BPG.*ADP)-ATP.*P3G)./...
      (p.PGK.KmPGKATP*p.PGK.KmPGKP3G.*(1+ADP./p.PGK.KmPGKADP + ATP./p.PGK.KmPGKATP).*(1+BPG/p.PGK.KmPGKBPG+P3G/p.PGK.KmPGKP3G));
  
%% GPM1 Consensus
% "3-phospho-D-glyceric acid" = "2-phospho-D-glyceric acid"; 
p.PGM_ExprsCor=1;
v_PGM=p.PGM_ExprsCor.*((((p.GPM1_kcat.*(f.GPM1+f.GPM2+f.GPM3))./p.GPM1_K3pg).*(P3G-P2G./p.GPM1_Keq))./...
    (1+P3G./p.GPM1_K3pg+1+P2G./p.GPM1_K2pg-1));

%% ENO1
% "2-phospho-D-glyceric acid" = phosphoenolpyruvate + water;  
p.ENO_ExprsCor=1; 
v_ENO=p.ENO_ExprsCor.*((((p.ENO1_kcat.*(f.ENO1+f.ENO2))./p.ENO1_K2pg).*(P2G-PEP./p.ENO1_Keq))./...
    (1+P2G./p.ENO1_K2pg+1+PEP./p.ENO1_Kpep-1));

%% PYK1 Consensus
% ADP + H+ + phosphoenolpyruvate -> ATP + pyruvate;  
v_PYK=p.PYK_ExprsCor.*((((p.PYK1_kcat.*(f.PYK1+f.PYK2))./(p.PYK1_Kadp.*p.PYK1_Kpep).*ADP.*PEP)./...
    ((1+ADP./p.PYK1_Kadp).*(1+PEP./p.PYK1_Kpep))).*...
    ((PEP./p.PYK1_Kpep+1).^p.PYK1_hill./(p.PYK1_L.*((ATP./p.PYK1_Katp+1)./(F16BP./p.PYK1_Kf16bp+1)).^p.PYK1_hill+(PEP./p.PYK1_Kpep+1).^p.PYK1_hill)));

%% PDC1 consensus
% pyruvate + H+ -> acetaldehyde + "carbon dioxide"
v_PDC=p.PDC_ExprsCor.*((p.PDC1_kcat.*(f.PDC1).*(PYR./p.PDC1_Kpyr).^p.PDC1_hill)./...
    (1+(PYR./p.PDC1_Kpyr).^p.PDC1_hill+PI./p.PDC1_Kpi));

%% ADH teusink
% acetaldehyde + H+ + NADH = ethanol + NAD+;  
v_ADH=-p.ADH_ExprsCor.*(p.ADH.VmADH./(p.ADH.KiADHNAD.*p.ADH.KmADHETOH).*(NAD.*ETOH-NADH.*ACE./p.ADH.KeqADH)./...
    (1+NAD./p.ADH.KiADHNAD+p.ADH.KmADHNAD.*ETOH./(p.ADH.KiADHNAD*p.ADH.KmADHETOH)+p.ADH.KmADHNADH.*ACE./(p.ADH.KiADHNADH.*p.ADH.KmADHACE)...
    +NADH./p.ADH.KiADHNADH+NAD.*ETOH./(p.ADH.KiADHNAD.*p.ADH.KmADHETOH)+p.ADH.KmADHNADH.*NAD.*ACE./(p.ADH.KiADHNAD.*p.ADH.KiADHNADH.*p.ADH.KmADHACE)...
    +p.ADH.KmADHNAD.*ETOH.*NADH./(p.ADH.KiADHNAD.*p.ADH.KmADHETOH.*p.ADH.KiADHNADH)+NADH.*ACE./(p.ADH.KiADHNADH.*p.ADH.KmADHACE)+...
    NAD.*ETOH.*ACE./(p.ADH.KiADHNAD*p.ADH.KmADHETOH.*p.ADH.KiADHACE)+ETOH.*NADH.*ACE./(p.ADH.KiADHETOH*p.ADH.KiADHNADH*p.ADH.KmADHACE)));

%% DOWNSTREAM METABOLISM
v_ACE=p.VmaxACE*ACE./(ACE+p.KmACE);%canelas.mtlD.v_PDC(exp)-canelas.mtlD.v_ADH(exp);
v_PDH=p.PDH_Vmax.*(PYR.^p.PDH_n)./(p.PDH_K50.^p.PDH_n+PYR.^p.PDH_n);

%% GLYCEROL BRANCH

%% GPD1 Consensus
% "glycerone phosphate" + H+ + NADH = NAD+ + "sn-glycerol 3-phosphate" 
v_G3PDH=((((p.GPD1_kcat.*f.GPD1)./(p.GPD1_Kdhap.*p.GPD1_Knadh)).*(DHAP.*NADH-(GLYC3P.*NAD)./p.GPD1_Keq))./...
    ((1+F16BP./p.GPD1_Kf16bp+ATP./p.GPD1_Katp+ADP./p.GPD1_Kadp).*(1+DHAP./p.GPD1_Kdhap+GLYC3P./p.GPD1_Kglyc3p).*(1+NADH./p.GPD1_Knadh+NAD./p.GPD1_Knad)));

%% HOR2 consensus
% %RHR2 consensus
% "sn-glycerol 3-phosphate" + water -> glycerol + phosphate; 
v_HOR2=(((p.HOR2_kcat.*f.HOR2)./p.HOR2_Kglyc3p.*GLYC3P)./...
    ((1+PI./p.HOR2_Kpi).*(1+GLYC3P./p.HOR2_Kglyc3p)));
v_RHR2=0;

%% TREHALOSE CYCLE

%% PGM1 consensus
% "D-glucose 1-phosphate" = "D-glucose 6-phosphate";  
v_PGM1=(((p.PGM1_kcat.*(f.PGM1+f.PGM2+f.PGM3))./p.PGM1_Kg1p.*(G1P-G6P./p.PGM1_Keq))./...
    (1+G1P./p.PGM1_Kg1p+G6P./p.PGM1_Kg6p));

%% UTP + "D-glucose 1-phosphate" + H+ -> UDP-D-glucose + diphosphate
v_UGP=0;

%% TPS1     Smallbone2011
% "D-glucose 6-phosphate" + UDP-D-glucose -> "alpha,alpha-trehalose 6-phosphate" + H+ + UDP;  
v_TPS1=(F6P./(F6P+p.TPS1_KmF6P)).*(((p.TPS1_kcat.*f.TPS1)./(p.TPS1_Kg6p.*p.TPS1_Kudp_glc).*G6P.*UDP_GLC./...
    ((1+G6P./p.TPS1_Kg6p).*(1+UDP_GLC./p.TPS1_Kudp_glc).*(1+PI./p.TPS1_Kpi))));


%% TPS2     Smallbone2011
% "alpha,alpha-trehalose 6-phosphate" + water -> alpha,alpha-trehalose + phosphate; 
v_TPS2=(((p.TPS2_kcat.*f.TPS2).*T6P.*PI)./...
    ((p.TPS2_Kt6p.*p.TPS2_Kpi)+(p.TPS2_Kt6p+T6P).*PI));

%% NTH1     Smallbone2011
% alpha,alpha-trehalose + water -> 2 * D-glucose;  
v_NTH1=(((p.NTH1_kcat.*f.NTH1)./p.NTH1_Ktre.*TRE)./...
    (1+TRE./p.NTH1_Ktre));

%% COFACTOR METABOLISM

%% vac_Pi_import
v_vacuolePi=p.vacuolePi_k.*(p.vacuolePi_steadyStatePi-PI);

%% ADK1 consensus
% 2 * ADP = AMP + ATP;
v_ADK1=p.ADK1_k.*((ADP.*ADP)-(AMP.*ATP)./p.ADK1_Keq);

%% mitoNADH
v_mitoNADH=p.mitoNADHVmax.*(NADH./(NADH+p.mitoNADHKm));

%% IXP cycle
% % % % % v_Amd1=(p.Amd1_Vmax.*AMP)./(p.Amd1_K50.*(1+PI./p.Amd1_Kpi)+AMP);
v_Amd1=(p.Amd1_Vmax.*AMP)./(p.Amd1_K50.*(1+PI./p.Amd1_Kpi)./(ATP./p.Amd1_Katp + 1)+AMP);
v_Ade13_v_Ade12=IMP.*p.Ade13_Ade12_k;
v_Isn1=IMP.*p.Isn1_k; 
v_Pnp1=INO.*p.Pnp1_k;
v_Hpt1=HYP.*p.Hpt1_k;

%% TRANSPORT REACTIONS (v_{GLT} not included)

%% Etanol transport
ETOHe=f.ETOH_e;
v_ETOHtransport=p.kETOHtransport*(ETOH-ETOHe);

%% Glycerol transport
GLYCEROLe=f.GLYCEROL_e;
v_GLYCEROLtransport=p.GlycerolTransport.*(GLYCEROL-GLYCEROLe);

%% SINK REACTIONS 
% Sink reactions when biomass is added. The equations still need to change.
% Km values adjusted ad hoc.
if setup.biomass == 1
    % G6P
    km_sinkG6P      = 0.01;
    poly_sinkG6P    = 3.6854 * d.^3 -   1.4119 * d.^2 -  0.6312 * d    - 0.0043; 
    if isfield(setup, 'sinkG6P_OFF') % option to shut down sink of G6P
        if setup.sinkG6P_OFF == 0
        elseif setup.sinkG6P_OFF == 1
            poly_sinkG6P = 0;
        end
    end
    ratioG6P        = 1;
    v_sinkG6P       = ratioG6P .* poly_sinkG6P .* (G6P./(G6P + km_sinkG6P));
    % F6P
    km_sinkF6P      = 0.0001;
    poly_sinkF6P    = 519.3740 * d.^6 - 447.7990 * d.^5 + 97.2843 * d.^4 + 8.0698 * d.^3 - 4.4005 * d.^2 + 0.6254 * d - 0.0078;
    v_sinkF6P       = poly_sinkF6P * (F6P./(F6P + km_sinkF6P));
    % GAP
    km_sinkGAP      = 0.0005;
    poly_sinkGAP    = 170.8447 * d.^6 - 113.2975 * d.^5 + 2.6494 * d.^4 + 10.2461 * d.^3 - 1.8002 * d.^2 + 0.1988 * d + 0.0012;
    v_sinkGAP       = poly_sinkGAP * (GLYCERAL3P./(GLYCERAL3P + km_sinkGAP)); 
    % P3G
    km_sinkP3G      = 0.001; 
    poly_sinkP3G    = -0.2381 * d.^2 -0.0210 * d   -0.0034;
    v_sinkP3G       = poly_sinkP3G * (P3G./(P3G + km_sinkP3G));
    % PEP
    km_sinkPEP      = 0.001;
    poly_sinkPEP    = -   0.0637 * d.^2 -   0.0617 * d   -  0.0008;
    v_sinkPEP       = poly_sinkPEP * (PEP./(PEP + km_sinkPEP));
    % PYR
    km_sinkPYR      = 0.001;
    poly_sinkPYR    = - 8.4853e+03 * d.^6 + 9.4027e+03 * d.^5 - 3.8027e+03 * d.^4 + 700.5 * d.^3 - 60.26 * d.^2 + 0.711 * d - 0.0356;
    if isfield(setup, 'oxygen_OFF') % in case no oxygen is available, assumed sink of pyruvate = 0.
        if setup.oxygen_OFF == 0
        elseif setup.oxygen_OFF == 1
            poly_sinkPYR = 0;
        end
    end
    v_sinkPYR       = poly_sinkPYR * (PYR./(PYR + km_sinkPYR));
    % ACE
    km_sinkACE      = 0.0001;
    poly_sinkACE    =     118.8562 * d.^6 - 352.3943 * d.^5 + 245.6092 * d.^4 - 75.2550 * d.^3 + 11.1153 * d.^2 - 1.0379 * d + 0.0119; 
    v_sinkACE       = poly_sinkACE * (ACE./(ACE + km_sinkACE));
else
    v_sinkG6P = 0;
    v_sinkF6P = 0;
    v_sinkGAP = 0;
    v_sinkP3G = 0;
    v_sinkPEP = 0;
    v_sinkPYR = 0;
    v_sinkACE = 0;
end

% option 2 for clamping trehalose
if setup.clamp10.TRE == 2 
    if((setup.conditionsSS == 0)&&(setup.conditionsGP == 1))
        v_PGM1 = zeros(size(v_PGM1));
        v_TPS1 = zeros(size(v_TPS1));
        v_TPS2 = zeros(size(v_TPS2));
        v_NTH1 = zeros(size(v_NTH1));
        v_vacuolePi = zeros(size(v_vacuolePi));
    end
end

%% vmito and vATPase
% %% ATPase
% v_ATPase=ATP./ADP.*p.ATPaseK;
% 
% %% mito(ATP)
% v_mito=p.mitoVmax.*ADP./(ADP+p.mitoADPKm).*(PI./(PI+p.mitoPiKm));
if((setup.conditionsSS == 0)&&(setup.conditionsGP == 1)) % glucose perturbation simulations
    % before 110 mM glucose pulse
    v_ATPase = 0.9306 .* ATP ./ ADP .* p.ATPase_ratio;
    v_mito=p.mitoVmax.*ADP./(ADP+p.mitoADPKm).*(PI./(PI+p.mitoPiKm));
    % during 110 mM glucose pulse (increased ATPase activity, value fixed in 'simulateGSvanHeerden_Y3M1.m')
    if((isfield(setup,'increase_Katpase'))&&(setup.increase_Katpase == 1))
        v_ATPase = ATP./ADP.*p.ATPaseK * 3.25 * p.ATPase_ratio2;
    end
else % steady state simulations
    RQ = [1.05 1.04 1.04 1.05 1.48 1.66 2.80 4.31];
    v_mito = (-v_sinkPYR * 3 + v_PDC * 1)/RQ(experiment)*0.95;
    GAMNGAM = [0.2639    0.5139    0.9306    1.7639    2.5972    2.8056    3.0139    3.2222];
    ATP_RATIO = [2.8280    3.2825    3.2965    3.7320    3.4540    3.3985    3.0315    2.9025];
    ADP_RATIO = [0.7765    0.8835    0.8925    1.0175    0.8315    0.8115    0.4825    0.4795];
    ADPATP_RATIO = ADP_RATIO ./ ATP_RATIO;
    v_ATPase = GAMNGAM(experiment) .* ATP ./ ADP .* ADPATP_RATIO(experiment);
end


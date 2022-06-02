function v = subsystemGAPDH_PGK(t,IC,p,f,experiment,canelas,setup)
%UNTITLED12 Summary of this function goes here
%   Detailed explanation goes here

% disp(t);
if setup.parEst.costfun == 1
    % %% STAGE 1. ASSOCIATE IC -> METABOLITE CONCENTRATIONS
    BPG=IC(1);
        ATP=canelas.mtlD.ATP(experiment);
        ADP=canelas.mtlD.ADP(experiment);
        NADH=canelas.mtlD.NADH(experiment);
        NAD=canelas.mtlD.NAD(experiment);
        GLYCERAL3P=canelas.mtlD.GAP(experiment);
        P3G=canelas.mtlD.P3G(experiment);
        PI=10;
        d=canelas.mtlD.D(experiment);
        PvHoek.D_GAPDH=[0.051867	0.072614	0.09751	0.114108	0.134855	0.155602	0.178423	0.19917	0.221992	0.248963	0.275934	0.29668	0.319502	0.348548	0.375519];
        PvHoek.GAPDH=[5.981308	5.358255	4.672897	4.23676	3.862928	3.489097	3.239875	3.05296	2.990654	3.115265	3.613707	3.925234	4.423676	5.109034	6.292835];
        PvHoek.D_PGK=[0.052523	0.072215	0.095844	0.119474	0.135221	0.154929	0.178566	0.202218	0.229816	0.255468	0.287065	0.306811	0.32855	0.34636	0.364177	0.382034];
        PvHoek.PGK=[7.979653	7.444776	6.791089	6.137403	5.662165	5.245629	4.651114	4.17494	3.639128	3.340233	3.158978	3.038297	3.035725	3.211133	3.445711	3.976145];
        p.GAPDH_ExprsCor=interp1(PvHoek.D_GAPDH,PvHoek.GAPDH,d,'pchip','extrap')./interp1(PvHoek.D_GAPDH,PvHoek.GAPDH,0.1,'pchip','extrap');
        p.PGK_ExprsCor=interp1(PvHoek.D_PGK,PvHoek.PGK,d,'pchip','extrap')./interp1(PvHoek.D_PGK,PvHoek.PGK,0.1,'pchip','extrap');

    % %% STAGE 2. CALCULATE REACTION RATES
    % describe reactions here
    
    %GAPDH
    GLYCERATE13BP=BPG;
    v_TDH1=p.GAPDH_ExprsCor.*((((p.TDH1_kcat.*(f.TDH1+f.TDH2+f.TDH3))./(p.TDH1_Kglyceral3p.*p.TDH1_Knad.*p.TDH1_Kpi)).*(GLYCERAL3P.*NAD.*PI-(GLYCERATE13BP.*NADH)./p.TDH1_Keq))./...
    ((1+GLYCERAL3P./p.TDH1_Kglyceral3p).*(1+NAD./p.TDH1_Knad).*(1+PI./p.TDH1_Kpi)+(1+GLYCERATE13BP./p.TDH1_Kglycerate13bp).*(1+NADH./p.TDH1_Knadh)-1));
    v_GAPDH=v_TDH1;

    %PGK
    v_PGK=p.PGK_ExprsCor.*p.PGK.VmPGK.*((p.PGK.KeqPGK.*BPG.*ADP)-ATP.*P3G)./...
    (p.PGK.KmPGKATP*p.PGK.KmPGKP3G.*(1+ADP./p.PGK.KmPGKADP + ATP./p.PGK.KmPGKATP).*(1+BPG/p.PGK.KmPGKBPG+P3G/p.PGK.KmPGKP3G));

    % %% STAGE 3. ASSOCIATE EACH RATE TO THE RESPECTIVE MASS BALANCE
    v(1)= +v_GAPDH - v_PGK; %BPG
end
v=v';
end
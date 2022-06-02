function v = subsystemGLT_HXK(t,IC,p,f,experiment,canelas,setup)
%UNTITLED12 Summary of this function goes here
%   Detailed explanation goes here

% disp(t);
% %% STAGE 1. ASSOCIATE IC -> METABOLITE CONCENTRATIONS
GLCi=IC(1);
    f.GLCo=canelas.mtlD.Cs(experiment);
    ATP=canelas.mtlD.ATP(experiment);
    ADP=canelas.mtlD.ADP(experiment);
    T6P=canelas.mtlD.T6P(experiment);
    G6P=canelas.mtlD.G6P(experiment);
    d=canelas.mtlD.D(experiment);
    PvHoek.D_HXK=[0.048535	0.065089	0.079577	0.094065	0.112703	0.135481	0.154118	0.176879	0.199657	0.222436	0.245214	0.270084	0.303244	0.332263	0.34886	0.375839];
    PvHoek.HXK=[2.612394 2.455598 2.330085 2.204571 2.078798 1.921615 1.795843 1.576355 1.419172 1.261989 1.104806 1.009798 0.883121 0.787855 0.78682 0.816292];
    p.HXK_ExprsCor=interp1(PvHoek.D_HXK,PvHoek.HXK,d,'pchip','extrap')./interp1(PvHoek.D_HXK,PvHoek.HXK,0.1,'pchip','extrap');

% %% STAGE 2. CALCULATE REACTION RATES
% describe reactions here

% GLT
v_GLT=p.GLT.VmGLT.*(f.GLCo-GLCi./p.GLT.KeqGLT)./...
    (p.GLT.KmGLTGLCo.*(1+f.GLCo/p.GLT.KmGLTGLCo+GLCi / p.GLT.KmGLTGLCi + ...
    0.91.*f.GLCo.*GLCi./(p.GLT.KmGLTGLCi.*p.GLT.KmGLTGLCo)));
% HXK1 
v_GLK=p.HXK_ExprsCor.*(((p.HXK1_kcat.*(f.HXK1+f.HXK2))./(p.HXK1_Katp.*p.HXK1_Kglc).*(ATP.*GLCi-((ADP.*G6P)./p.HXK1_Keq)))./...
    ((1+ATP./p.HXK1_Katp+ADP./p.HXK1_Kadp).*(1+GLCi./p.HXK1_Kglc+G6P./p.HXK1_Kg6p+T6P./p.HXK1_Kt6p)));

% %% STAGE 3. ASSOCIATE EACH RATE TO THE RESPECTIVE MASS BALANCE
v(1)= +v_GLT - v_GLK; %GLCi

v=v';
end
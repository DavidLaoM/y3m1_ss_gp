function v = ODE_PDC_ADH_BRANCHES(t,IC,p,f,experiment,canelas,x)

% set state variables
PEP=canelas.mtlD.PEP(experiment);
ATP=canelas.mtlD.ATP(experiment);
ADP=canelas.mtlD.ADP(experiment);
F16BP=canelas.mtlD.FBP(experiment);
P3G=canelas.mtlD.P3G(experiment);
DHAP=canelas.mtlD.DHAP(experiment);
NAD=canelas.mtlD.NAD(experiment);
NADH=canelas.mtlD.NADH(experiment);
PI=10;

d=canelas.mtlD.D(experiment);
PvHoek_EnzymeExpressionData
%     PvHoek.D_PDC=[0.047059	0.058824	0.068627	0.086275	0.1	0.115686	0.139216	0.170588	0.201961	0.229412	0.256863	0.276471	0.288235	0.301961	0.309804	0.319608	0.32549	0.333333	0.341176	0.35098	0.358824	0.376471];
%     PvHoek.PDC=[0.792899	0.751479	0.704142	0.650888	0.60355	0.568047	0.544379	0.52071	0.514793	0.514793	0.532544	0.579882	0.650888	0.751479	0.840237	0.952663	1.04142	1.147929	1.254438	1.378698	1.443787	1.609467];
%     PvHoek.D_ADH=[0.051181	0.076772	0.098425	0.124016	0.149606	0.173228	0.204724	0.222441	0.257874	0.281496	0.311024	0.324803	0.344488	0.375984];
%     PvHoek.ADH=[9.881657	9.053254	8.343195	7.514793	6.745562	6.094675	5.088757	4.615385	3.727811	3.254438	2.781065	2.662722	2.95858	4.260355];
%     p.PDC_ExprsCor=interp1(PvHoek.D_PDC,PvHoek.PDC,d,'pchip','extrap')./interp1(PvHoek.D_PDC,PvHoek.PDC,0.1,'pchip','extrap');
%     p.ADH_ExprsCor=interp1(PvHoek.D_ADH,PvHoek.ADH,d,'pchip','extrap')./interp1(PvHoek.D_ADH,PvHoek.ADH,0.1,'pchip','extrap');

PYR=IC(1);
ACE=IC(2);
ETOH=IC(3);
GLYC3P=IC(4);
GLYCEROL=IC(5);

rateEquations_PDC_ADH_BRANCHES
v_GAPDH=canelas.mtlD.v_GAPDH(experiment);
v_PYK=canelas.mtlD.v_PYK(experiment);
v_G3PDH=canelas.mtlD.v_G3PDH(experiment);

v(1) = + v_PYK - v_PDC - v_PDH; %PYR
v(2) = + v_PDC - v_ADH - v_ACE; %ACE
v(3) = + v_ADH - v_ETOHtransport; %ETOH
v(4) = + v_G3PDH - v_HOR2; %G3P
v(5) = + v_HOR2 - v_GLYCEROLtransport; %GLYCEROL

v=v';
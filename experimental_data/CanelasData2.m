% When checking the supplmentary material from the Canelas2011 paper, most
% of the range qS = 0.0-3.0 micromol/gDW/s is covered by the hereby
% dilution ate D=0.02-0.375.
% The numbers that were on top of this and showed a different tendency in
% that work are not there.
% The 16 datapoints contain an overlap (ramp and the other experiment) that
% has been deleted here. We work with 8 datapoints.

% In the end of this file, the sink reactions are also calculated.

% summary units
    % intracellular metabolites are displayedin [mM]
    % Cs, glycerol_e as well in [mM]
    % DW in [g/L]
    % D in [h^{-1}]
    % ETOH(_e) in [mM] (instead of [mM/gdw/h], which was the reaction rate)


%% Unpolished dataset
% concentrations
canelas.mtlD_original.D=[0.020	0.021	0.050	0.051	0.100	0.101	0.200	0.201	0.300	0.301	0.325	0.326	0.350	0.351	0.375	0.376]; %1/h
canelas.mtlD_original.Cs=[0.13	0.13	0.15	0.14	0.18	0.19	0.19	0.18	0.45	0.46	0.55	0.54	1.44	1.43	3.27	3.31]; %mM
canelas.mtlD_original.Cs_sd=[0.020	0.020	0.041	0.031	0.018	0.017	0.028	0.018	0.037	0.032	0.030	0.048	0.095	0.072	0.142	0.150]; %mM
canelas.mtlD_original.DW=[3.29	3.29	3.29	3.29	3.29	3.29	3.29	3.29	3.29	3.29	3.29	3.29	3.29	3.29	3.29	3.29]; %g/L
canelas.mtlD_original.ETOH=[0.00	0.00	0.00	0.00	0.00	0.00	0.00	0.00	28.39	26.66	33.64	30.66	48.23	40.81	54.05	54.61]; %mM (mM/gdw/h)
canelas.mtlD_original.GLYCEROL_e=[0.12	0.12	0.14	0.14	0.09	0.09	0.07	0.07	0.06	0.06	0.05	0.04	0.13	0.13	0.23	0.23]; %mM

canelas.mtlD_original.G6P=[2.13	2.09	3.58	3.65	4.99	5.03	9.18	9.00	8.91	9.02	9.24	8.97	6.69	6.53	7.01	6.92]/2; % converted from micromol/gdw to mM using conversion factor 1.7
canelas.mtlD_original.G6P_sd=[0.16	0.18	0.24	0.24	0.24	0.44	0.49	0.39	0.69	0.46	0.54	0.50	0.43	0.59	0.46	0.56]/2;% converted from micromol/gdw to mM using conversion factor 1.7
canelas.mtlD_original.PEP=[2.46	2.34	2.28	2.31	2.68	2.73	2.01	1.87	0.53	0.52	0.34	0.32	0.17	0.17	0.19	0.18]/2;% converted from micromol/gdw to mM using conversion factor 1.7
canelas.mtlD_original.PEP_sd=[0.14	0.09	0.06	0.07	0.09	0.06	0.07	0.12	0.03	0.02	0.02	0.01	0.02	0.02	0.05	0.02]/2;% converted from micromol/gdw to mM using conversion factor 1.7
canelas.mtlD_original.G1P=[0.20	0.20	0.20	0.20	0.26	0.27	0.45	0.43	0.34	0.34	0.33	0.30	0.19	0.21	0.20	0.17]/2;% converted from micromol/gdw to mM using conversion factor 1.7
canelas.mtlD_original.G1P_sd=[0.03	0.02	0.01	0.03	0.00	0.00	0.03	0.00	0.02	0.03	0.02	0.02	0.02	0.02	0.02	0.05]/2;% converted from micromol/gdw to mM using conversion factor 1.7
canelas.mtlD_original.G3P=[0.09	0.08	0.10	0.11	0.06	0.07	0.06	0.08	0.08	0.07	0.09	0.08	0.13	0.14	0.22	0.22]/2;% converted from micromol/gdw to mM using conversion factor 1.7
canelas.mtlD_original.G3P_sd=[0.02	0.02	0.02	0.02	0.01	0.02	0.01	0.01	0.02	0.02	0.02	0.01	0.01	0.02	0.02	0.01]/2;% converted from micromol/gdw to mM using conversion factor 1.7
canelas.mtlD_original.PYR=[0.43	0.48	0.64	0.72	0.65	0.82	0.97	0.83	1.83	1.86	2.51	1.87	5.59	5.21	7.4	6.6]/2;% converted from micromol/gdw to mM using conversion factor 1.7
canelas.mtlD_original.PYR_sd=[0.26	0.32	0.32	0.38	0.32	0.71	0.27	0.40	0.67	0.42	0.20	0.60	0.55	1.11	1.00	1.85]/2;% converted from micromol/gdw to mM using conversion factor 1.7
canelas.mtlD_original.GAP=[0.000	0.000	0.019	0.020	0.020	0.022	0.020	0.027	0.036	0.028	0.032	0.028	0.037	0.046	0.043	0.055]/2;% converted from micromol/gdw to mM using conversion factor 1.7
canelas.mtlD_original.GAP_sd=[0.000	0.000	0.005	0.005	0.005	0.005	0.005	0.007	0.009	0.007	0.008	0.007	0.009	0.011	0.011	0.014]/2;% converted from micromol/gdw to mM using conversion factor 1.7
canelas.mtlD_original.P2G=[0.578	0.578	0.530	0.553	0.682	0.711	0.607	0.596	0.260	0.262	0.211	0.188	0.141	0.131	0.136	0.146]/2;% converted from micromol/gdw to mM using conversion factor 1.7
canelas.mtlD_original.P2G_sd=[0.002	0.006	0.002	0.006	0.000	0.008	0.001	0.017	0.031	0.006	0.007	0.008	0.013	0.010	0.013	0.000]/2;% converted from micromol/gdw to mM using conversion factor 1.7
canelas.mtlD_original.DHAP=[0.172	0.159	0.285	0.291	0.363	0.356	0.676	0.705	0.783	0.839	0.909	0.857	1.209	1.137	1.373	1.516]/2;% converted from micromol/gdw to mM using conversion factor 1.7
canelas.mtlD_original.DHAP_sd=[0.006	0.006	0.004	0.008	0.004	0.005	0.043	0.030	0.013	0.040	0.014	0.005	0.006	0.005	0.019	0.134]/2;% converted from micromol/gdw to mM using conversion factor 1.7
canelas.mtlD_original.P3G=[5.01	4.87	4.67	4.90	6.03	6.32	5.83	5.79	3.03	2.80	2.74	2.61	1.94	1.98	2.0	1.8]/2;% converted from micromol/gdw to mM using conversion factor 1.7
canelas.mtlD_original.P3G_sd=[0.032	0.097	0.043	0.020	0.078	0.172	0.014	0.018	0.044	0.036	0.001	0.037	0.080	0.117	0.025	0.083]/2;% converted from micromol/gdw to mM using conversion factor 1.7
canelas.mtlD_original.F6P=[0.515	0.497	0.922	0.947	1.320	1.365	2.368	2.322	1.874	1.806	1.770	1.702	0.948	0.915	0.876	0.859]/2;% converted from micromol/gdw to mM using conversion factor 1.7
canelas.mtlD_original.F6P_sd=[0.003	0.024	0.039	0.028	0.031	0.039	0.041	0.056	0.086	0.047	0.072	0.062	0.029	0.038	0.034	0.045]/2;% converted from micromol/gdw to mM using conversion factor 1.7
canelas.mtlD_original.FBP=[0.16	0.15	0.30	0.28	0.41	0.41	0.85	0.85	2.37	2.42	3.14	3.13	6.13	5.84	8.48	8.56]/2;% converted from micromol/gdw to mM using conversion factor 1.7
canelas.mtlD_original.FBP_sd=[0.04	0.01	0.04	0.03	0.03	0.04	0.05	0.04	0.11	0.12	0.16	0.18	0.25	0.30	0.42	0.28]/2;% converted from micromol/gdw to mM using conversion factor 1.7
canelas.mtlD_original.T6P=[0.31	0.30	0.37	0.37	0.20	0.21	0.11	0.09	0.02	0.02	0.02	0.02	0.03	0.03	0.03	0.03]/2;% converted from micromol/gdw to mM using conversion factor 1.7
canelas.mtlD_original.T6P_sd=[0.010	0.013	0.012	0.019	0.009	0.015	0.007	0.005	0.001	0.002	0.003	0.003	0.005	0.007	0.003	0.004]/2;% converted from micromol/gdw to mM using conversion factor 1.7
canelas.mtlD_original.F26BP=[0.0020	0.0007	0.0006	0.0006	0.0021	0.0020	0.0027	0.0017	0.0026	0.0048	0.0062	0.0034	0.0095	0.0079	0.0119	0.0108]/2;
canelas.mtlD_original.F26BP_sd=[0.0013	0.0003	0.0007	0.0007	0.0007	0.0007	0.0031	0.0016	0.0008	0.0024	0.0012	0.0002	0.0018	0.0026	0.0001	0.0035]/2;
canelas.mtlD_original.AMP=[0.458	0.436	0.491	0.542	0.517	0.548	0.602	0.568	0.375	0.389	0.400	0.342	0.186	0.174	0.224	0.196]/2;% converted from micromol/gdw to mM using conversion factor 1.7
canelas.mtlD_original.AMP_sd=[0.05	0.02	0.02	0.06	0.05	0.06	0.03	0.03	0.04	0.07	0.05	0.04	0.03	0.02	0.03	0.02]/2;% converted from micromol/gdw to mM using conversion factor 1.7
canelas.mtlD_original.ADP=[1.553	1.525	1.767	1.821	1.785	1.820	2.035	2.020	1.663	1.542	1.623	1.656	0.965	1.000	0.959	0.946]/2;% converted from micromol/gdw to mM using conversion factor 1.7
canelas.mtlD_original.ADP_sd=[0.03	0.02	0.04	0.03	0.05	0.08	0.04	0.02	0.03	0.06	0.04	0.10	0.04	0.02	0.06	0.03]/2;% converted from micromol/gdw to mM using conversion factor 1.7
canelas.mtlD_original.ATP=[5.656	5.535	6.565	6.648	6.593	6.909	7.464	7.362	6.908	7.069	6.797	6.471	6.063	5.812	5.805	5.894]/2;% converted from micromol/gdw to mM using conversion factor 1.7
canelas.mtlD_original.ATP_sd=[0.11	0.06	0.11	0.10	0.09	0.27	0.20	0.18	0.09	0.08	0.07	0.44	0.09	0.07	0.09	0.14]/2;% converted from micromol/gdw to mM using conversion factor 1.7
canelas.mtlD_original.NAD_NADHratio=[154.9	152.7	147.8	148.5	149.3	150.3	121.5	118.7	75.4	75.3	67.5	67.2	59.9	59.8	55.0	54.7];
sumNAD_NADH=1.59;
canelas.mtlD_original.NADH=sumNAD_NADH./(canelas.mtlD_original.NAD_NADHratio+1);
canelas.mtlD_original.NAD=sumNAD_NADH-canelas.mtlD_original.NADH;
clear sumNAD_NADH, canelas.mtlD_original.sumNAD_NADH=1.59;

% fluxes
canelas.mtlD_original.v_HK=[0.069	0.069	0.166	0.166	0.307	0.307	0.597	0.597	1.297	1.297	1.511	1.511	2.114	2.114	2.840	2.840]/2; % converted from micromol/gdw/s to mM/s using conversion factor 1.7
canelas.mtlD_original.v_PGI=[0.04	0.037495741	0.08	0.084548105	0.15	0.151074306	0.29	0.287532924	0.85	0.851711027	1.04	1.036973443	1.64	1.637588001	2.35	2.351895735]/2;
canelas.mtlD_original.v_PFK=[0.046	0.046	0.108	0.108	0.199	0.199	0.393	0.393	1.009	1.009	1.205	1.205	1.800	1.800	2.512	2.512]/2;
canelas.mtlD_original.v_FBA=[0.046	0.046	0.108	0.108	0.199	0.199	0.393	0.393	1.009	1.009	1.205	1.205	1.800	1.800	2.512	2.512]/2;
canelas.mtlD_original.v_TPI=[0.05	0.045187394	0.11	0.105928086	0.20	0.195464041	0.39	0.387108575	1.00	1.001056189	1.20	1.197344322	1.79	1.786042241	2.49	2.486176935]/2;
canelas.mtlD_original.v_GAPDH=[0.10	0.096592845	0.23	0.228620019	0.42	0.423828708	0.84	0.843576236	2.10	2.104985213	2.50	2.504578755	3.68	3.683807775	5.10	5.096761453]/2;
canelas.mtlD_original.v_PGK=[0.10	0.096592845	0.23	0.228620019	0.42	0.423828708	0.84	0.843576236	2.10	2.104985213	2.50	2.504578755	3.68	3.683807775	5.10	5.096761453]/2;
canelas.mtlD_original.v_PGM=[0.09	0.093356048	0.22	0.219630709	0.41	0.405252164	0.80	0.801872988	2.04	2.035276722	2.43	2.427884615	3.60	3.601163147	5.01	5.005924171]/2;
canelas.mtlD_original.v_ENO=[0.09	0.093356048	0.22	0.219630709	0.41	0.405252164	0.80	0.801872988	2.04	2.035276722	2.43	2.427884615	3.60	3.601163147	5.01	5.005924171]/2;
canelas.mtlD_original.v_PYK=[0.09	0.091056218	0.21	0.213394882	0.39	0.392420173	0.77	0.77260755	1.99	1.987748204	2.38	2.376373626	3.54	3.544536272	4.94	4.942733017]/2;
canelas.mtlD_original.v_G3PDH=[0.001	0.001	0.002	0.002	0.003	0.003	0.006	0.006	0.008	0.008	0.009	0.009	0.014	0.014	0.025	0.025]/2;
canelas.mtlD_original.v_GLT=[0.07	0.068551959	0.17	0.166261743	0.31	0.306550283	0.60	0.59650278	1.30	1.297000422	1.51	1.510989011	2.11	2.113559841	2.84	2.839652449]/2;
canelas.mtlD_original.v_PDC=[0.01	0.013321976	0.04	0.035050211	0.07	0.066062369	0.13	0.129499561	1.07	1.074144487	1.43	1.428571429	2.74	2.744107744	4.19	4.190363349]/2;
canelas.mtlD_original.v_ADH=[0.000	0.000	0.000	0.000	0.000	0.000	0.000	0.000	0.844	0.844	1.172	1.172	2.410	2.410	3.774	3.774]/2;

% figure(2)
% plot(canelas.mtlD.D, canelas.mtlD.Cs, '*', canelas.mtlD.D, canelas.mtlD.DW, 'o')
% A = struct2cell(canelas.mtlD)
% figure(3)
% for i = 1:55
% subplot(5,11,i)
% plot(A{1}, A{i})
% end

%% Polished dataset

% % concentrations
canelas.mtlD.D          = [canelas.mtlD_original.D(1)           canelas.mtlD_original.D(3)      canelas.mtlD_original.D(5)      canelas.mtlD_original.D(7)      canelas.mtlD_original.D(9)      canelas.mtlD_original.D(11)         canelas.mtlD_original.D(13)         canelas.mtlD_original.D(15)];
canelas.mtlD.Cs         = [canelas.mtlD_original.Cs(1)          canelas.mtlD_original.Cs(3)      canelas.mtlD_original.Cs(5)      canelas.mtlD_original.Cs(7)      canelas.mtlD_original.Cs(9)      canelas.mtlD_original.Cs(11)         canelas.mtlD_original.Cs(13)         canelas.mtlD_original.Cs(15)];
canelas.mtlD.Cs_sd      = [canelas.mtlD_original.Cs_sd(1)       canelas.mtlD_original.Cs_sd(3)      canelas.mtlD_original.Cs_sd(5)      canelas.mtlD_original.Cs_sd(7)      canelas.mtlD_original.Cs_sd(9)      canelas.mtlD_original.Cs_sd(11)         canelas.mtlD_original.Cs_sd(13)         canelas.mtlD_original.Cs_sd(15)];
canelas.mtlD.DW         = [canelas.mtlD_original.DW(1)          canelas.mtlD_original.DW(3)      canelas.mtlD_original.DW(5)      canelas.mtlD_original.DW(7)      canelas.mtlD_original.DW(9)      canelas.mtlD_original.DW(11)         canelas.mtlD_original.DW(13)         canelas.mtlD_original.DW(15)];
canelas.mtlD.ETOH       = [canelas.mtlD_original.ETOH(1)        canelas.mtlD_original.ETOH(3)      canelas.mtlD_original.ETOH(5)      canelas.mtlD_original.ETOH(7)      canelas.mtlD_original.ETOH(9)      canelas.mtlD_original.ETOH(11)         canelas.mtlD_original.ETOH(13)         canelas.mtlD_original.ETOH(15)];
canelas.mtlD.GLYCEROL_e = [canelas.mtlD_original.GLYCEROL_e(1)  canelas.mtlD_original.GLYCEROL_e(3)      canelas.mtlD_original.GLYCEROL_e(5)      canelas.mtlD_original.GLYCEROL_e(7)      canelas.mtlD_original.GLYCEROL_e(9)      canelas.mtlD_original.GLYCEROL_e(11)         canelas.mtlD_original.GLYCEROL_e(13)         canelas.mtlD_original.GLYCEROL_e(15)];

canelas.mtlD.G6P    = [canelas.mtlD_original.G6P(1)     canelas.mtlD_original.G6P(3)      canelas.mtlD_original.G6P(5)      canelas.mtlD_original.G6P(7)      canelas.mtlD_original.G6P(9)      canelas.mtlD_original.G6P(11)         canelas.mtlD_original.G6P(13)         canelas.mtlD_original.G6P(15)];
canelas.mtlD.G6P_sd    = [canelas.mtlD_original.G6P_sd(1)     canelas.mtlD_original.G6P_sd(3)      canelas.mtlD_original.G6P_sd(5)      canelas.mtlD_original.G6P_sd(7)      canelas.mtlD_original.G6P_sd(9)      canelas.mtlD_original.G6P_sd(11)         canelas.mtlD_original.G6P_sd(13)         canelas.mtlD_original.G6P_sd(15)];
canelas.mtlD.PEP    = [canelas.mtlD_original.PEP(1)     canelas.mtlD_original.PEP(3)      canelas.mtlD_original.PEP(5)      canelas.mtlD_original.PEP(7)      canelas.mtlD_original.PEP(9)      canelas.mtlD_original.PEP(11)         canelas.mtlD_original.PEP(13)         canelas.mtlD_original.PEP(15)];
canelas.mtlD.PEP_sd    = [canelas.mtlD_original.PEP_sd(1)     canelas.mtlD_original.PEP_sd(3)      canelas.mtlD_original.PEP_sd(5)      canelas.mtlD_original.PEP_sd(7)      canelas.mtlD_original.PEP_sd(9)      canelas.mtlD_original.PEP_sd(11)         canelas.mtlD_original.PEP_sd(13)         canelas.mtlD_original.PEP_sd(15)];
canelas.mtlD.G1P    = [canelas.mtlD_original.G1P(1)     canelas.mtlD_original.G1P(3)      canelas.mtlD_original.G1P(5)      canelas.mtlD_original.G1P(7)      canelas.mtlD_original.G1P(9)      canelas.mtlD_original.G1P(11)         canelas.mtlD_original.G1P(13)         canelas.mtlD_original.G1P(15)];
canelas.mtlD.G1P_sd    = [canelas.mtlD_original.G1P_sd(1)     canelas.mtlD_original.G1P_sd(3)      canelas.mtlD_original.G1P_sd(5)      canelas.mtlD_original.G1P_sd(7)      canelas.mtlD_original.G1P_sd(9)      canelas.mtlD_original.G1P_sd(11)         canelas.mtlD_original.G1P_sd(13)         canelas.mtlD_original.G1P_sd(15)];
canelas.mtlD.G3P    = [canelas.mtlD_original.G3P(1)     canelas.mtlD_original.G3P(3)      canelas.mtlD_original.G3P(5)      canelas.mtlD_original.G3P(7)      canelas.mtlD_original.G3P(9)      canelas.mtlD_original.G3P(11)         canelas.mtlD_original.G3P(13)         canelas.mtlD_original.G3P(15)];
canelas.mtlD.G3P_sd    = [canelas.mtlD_original.G3P_sd(1)     canelas.mtlD_original.G3P_sd(3)      canelas.mtlD_original.G3P_sd(5)      canelas.mtlD_original.G3P_sd(7)      canelas.mtlD_original.G3P_sd(9)      canelas.mtlD_original.G3P_sd(11)         canelas.mtlD_original.G3P_sd(13)         canelas.mtlD_original.G3P_sd(15)];
canelas.mtlD.PYR    = [canelas.mtlD_original.PYR(1)     canelas.mtlD_original.PYR(3)      canelas.mtlD_original.PYR(5)      canelas.mtlD_original.PYR(7)      canelas.mtlD_original.PYR(9)      canelas.mtlD_original.PYR(11)         canelas.mtlD_original.PYR(13)         canelas.mtlD_original.PYR(15)];
canelas.mtlD.PYR_sd    = [canelas.mtlD_original.PYR_sd(1)     canelas.mtlD_original.PYR_sd(3)      canelas.mtlD_original.PYR_sd(5)      canelas.mtlD_original.PYR_sd(7)      canelas.mtlD_original.PYR_sd(9)      canelas.mtlD_original.PYR_sd(11)         canelas.mtlD_original.PYR_sd(13)         canelas.mtlD_original.PYR_sd(15)];
canelas.mtlD.GAP    = [canelas.mtlD_original.GAP(1)     canelas.mtlD_original.GAP(3)      canelas.mtlD_original.GAP(5)      canelas.mtlD_original.GAP(7)      canelas.mtlD_original.GAP(9)      canelas.mtlD_original.GAP(11)         canelas.mtlD_original.GAP(13)         canelas.mtlD_original.GAP(15)];
canelas.mtlD.GAP_sd    = [canelas.mtlD_original.GAP_sd(1)     canelas.mtlD_original.GAP_sd(3)      canelas.mtlD_original.GAP_sd(5)      canelas.mtlD_original.GAP_sd(7)      canelas.mtlD_original.GAP_sd(9)      canelas.mtlD_original.GAP_sd(11)         canelas.mtlD_original.GAP_sd(13)         canelas.mtlD_original.GAP_sd(15)];
canelas.mtlD.P2G    = [canelas.mtlD_original.P2G(1)     canelas.mtlD_original.P2G(3)      canelas.mtlD_original.P2G(5)      canelas.mtlD_original.P2G(7)      canelas.mtlD_original.P2G(9)      canelas.mtlD_original.P2G(11)         canelas.mtlD_original.P2G(13)         canelas.mtlD_original.P2G(15)];
canelas.mtlD.P2G_sd    = [canelas.mtlD_original.P2G_sd(1)     canelas.mtlD_original.P2G_sd(3)      canelas.mtlD_original.P2G_sd(5)      canelas.mtlD_original.P2G_sd(7)      canelas.mtlD_original.P2G_sd(9)      canelas.mtlD_original.P2G_sd(11)         canelas.mtlD_original.P2G_sd(13)         canelas.mtlD_original.P2G_sd(15)];
canelas.mtlD.DHAP    = [canelas.mtlD_original.DHAP(1)     canelas.mtlD_original.DHAP(3)      canelas.mtlD_original.DHAP(5)      canelas.mtlD_original.DHAP(7)      canelas.mtlD_original.DHAP(9)      canelas.mtlD_original.DHAP(11)         canelas.mtlD_original.DHAP(13)         canelas.mtlD_original.DHAP(15)];
canelas.mtlD.DHAP_sd    = [canelas.mtlD_original.DHAP_sd(1)     canelas.mtlD_original.DHAP_sd(3)      canelas.mtlD_original.DHAP_sd(5)      canelas.mtlD_original.DHAP_sd(7)      canelas.mtlD_original.DHAP_sd(9)      canelas.mtlD_original.DHAP_sd(11)         canelas.mtlD_original.DHAP_sd(13)         canelas.mtlD_original.DHAP_sd(15)];
canelas.mtlD.P3G    = [canelas.mtlD_original.P3G(1)     canelas.mtlD_original.P3G(3)      canelas.mtlD_original.P3G(5)      canelas.mtlD_original.P3G(7)      canelas.mtlD_original.P3G(9)      canelas.mtlD_original.P3G(11)         canelas.mtlD_original.P3G(13)         canelas.mtlD_original.P3G(15)];
canelas.mtlD.P3G_sd    = [canelas.mtlD_original.P3G_sd(1)     canelas.mtlD_original.P3G_sd(3)      canelas.mtlD_original.P3G_sd(5)      canelas.mtlD_original.P3G_sd(7)      canelas.mtlD_original.P3G_sd(9)      canelas.mtlD_original.P3G_sd(11)         canelas.mtlD_original.P3G_sd(13)         canelas.mtlD_original.P3G_sd(15)];
canelas.mtlD.F6P    = [canelas.mtlD_original.F6P(1)     canelas.mtlD_original.F6P(3)      canelas.mtlD_original.F6P(5)      canelas.mtlD_original.F6P(7)      canelas.mtlD_original.F6P(9)      canelas.mtlD_original.F6P(11)         canelas.mtlD_original.F6P(13)         canelas.mtlD_original.F6P(15)];
canelas.mtlD.F6P_sd    = [canelas.mtlD_original.F6P_sd(1)     canelas.mtlD_original.F6P_sd(3)      canelas.mtlD_original.F6P_sd(5)      canelas.mtlD_original.F6P_sd(7)      canelas.mtlD_original.F6P_sd(9)      canelas.mtlD_original.F6P_sd(11)         canelas.mtlD_original.F6P_sd(13)         canelas.mtlD_original.F6P_sd(15)];
canelas.mtlD.FBP    = [canelas.mtlD_original.FBP(1)     canelas.mtlD_original.FBP(3)      canelas.mtlD_original.FBP(5)      canelas.mtlD_original.FBP(7)      canelas.mtlD_original.FBP(9)      canelas.mtlD_original.FBP(11)         canelas.mtlD_original.FBP(13)         canelas.mtlD_original.FBP(15)];
canelas.mtlD.FBP_sd    = [canelas.mtlD_original.FBP_sd(1)     canelas.mtlD_original.FBP_sd(3)      canelas.mtlD_original.FBP_sd(5)      canelas.mtlD_original.FBP_sd(7)      canelas.mtlD_original.FBP_sd(9)      canelas.mtlD_original.FBP_sd(11)         canelas.mtlD_original.FBP_sd(13)         canelas.mtlD_original.FBP_sd(15)];
canelas.mtlD.T6P    = [canelas.mtlD_original.T6P(1)     canelas.mtlD_original.T6P(3)      canelas.mtlD_original.T6P(5)      canelas.mtlD_original.T6P(7)      canelas.mtlD_original.T6P(9)      canelas.mtlD_original.T6P(11)         canelas.mtlD_original.T6P(13)         canelas.mtlD_original.T6P(15)];
canelas.mtlD.T6P_sd    = [canelas.mtlD_original.T6P_sd(1)     canelas.mtlD_original.T6P_sd(3)      canelas.mtlD_original.T6P_sd(5)      canelas.mtlD_original.T6P_sd(7)      canelas.mtlD_original.T6P_sd(9)      canelas.mtlD_original.T6P_sd(11)         canelas.mtlD_original.T6P_sd(13)         canelas.mtlD_original.T6P_sd(15)];
canelas.mtlD.F26BP    = [canelas.mtlD_original.F26BP(1)     canelas.mtlD_original.F26BP(3)      canelas.mtlD_original.F26BP(5)      canelas.mtlD_original.F26BP(7)      canelas.mtlD_original.F26BP(9)      canelas.mtlD_original.F26BP(11)         canelas.mtlD_original.F26BP(13)         canelas.mtlD_original.F26BP(15)];
canelas.mtlD.F26BP_sd    = [canelas.mtlD_original.F26BP_sd(1)     canelas.mtlD_original.F26BP_sd(3)      canelas.mtlD_original.F26BP_sd(5)      canelas.mtlD_original.F26BP_sd(7)      canelas.mtlD_original.F26BP_sd(9)      canelas.mtlD_original.F26BP_sd(11)         canelas.mtlD_original.F26BP_sd(13)         canelas.mtlD_original.F26BP_sd(15)];
canelas.mtlD.AMP    = [canelas.mtlD_original.AMP(1)     canelas.mtlD_original.AMP(3)      canelas.mtlD_original.AMP(5)      canelas.mtlD_original.AMP(7)      canelas.mtlD_original.AMP(9)      canelas.mtlD_original.AMP(11)         canelas.mtlD_original.AMP(13)         canelas.mtlD_original.AMP(15)];
canelas.mtlD.AMP_sd    = [canelas.mtlD_original.AMP_sd(1)     canelas.mtlD_original.AMP_sd(3)      canelas.mtlD_original.AMP_sd(5)      canelas.mtlD_original.AMP_sd(7)      canelas.mtlD_original.AMP_sd(9)      canelas.mtlD_original.AMP_sd(11)         canelas.mtlD_original.AMP_sd(13)         canelas.mtlD_original.AMP_sd(15)];
canelas.mtlD.ADP    = [canelas.mtlD_original.ADP(1)     canelas.mtlD_original.ADP(3)      canelas.mtlD_original.ADP(5)      canelas.mtlD_original.ADP(7)      canelas.mtlD_original.ADP(9)      canelas.mtlD_original.ADP(11)         canelas.mtlD_original.ADP(13)         canelas.mtlD_original.ADP(15)];
canelas.mtlD.ADP_sd    = [canelas.mtlD_original.ADP_sd(1)     canelas.mtlD_original.ADP_sd(3)      canelas.mtlD_original.ADP_sd(5)      canelas.mtlD_original.ADP_sd(7)      canelas.mtlD_original.ADP_sd(9)      canelas.mtlD_original.ADP_sd(11)         canelas.mtlD_original.ADP_sd(13)         canelas.mtlD_original.ADP_sd(15)];
canelas.mtlD.ATP    = [canelas.mtlD_original.ATP(1)     canelas.mtlD_original.ATP(3)      canelas.mtlD_original.ATP(5)      canelas.mtlD_original.ATP(7)      canelas.mtlD_original.ATP(9)      canelas.mtlD_original.ATP(11)         canelas.mtlD_original.ATP(13)         canelas.mtlD_original.ATP(15)];
canelas.mtlD.ATP_sd    = [canelas.mtlD_original.ATP_sd(1)     canelas.mtlD_original.ATP_sd(3)      canelas.mtlD_original.ATP_sd(5)      canelas.mtlD_original.ATP_sd(7)      canelas.mtlD_original.ATP_sd(9)      canelas.mtlD_original.ATP_sd(11)         canelas.mtlD_original.ATP_sd(13)         canelas.mtlD_original.ATP_sd(15)];
canelas.mtlD.NAD_NADHratio    = [canelas.mtlD_original.NAD_NADHratio(1)     canelas.mtlD_original.NAD_NADHratio(3)      canelas.mtlD_original.NAD_NADHratio(5)      canelas.mtlD_original.NAD_NADHratio(7)      canelas.mtlD_original.NAD_NADHratio(9)      canelas.mtlD_original.NAD_NADHratio(11)         canelas.mtlD_original.NAD_NADHratio(13)         canelas.mtlD_original.NAD_NADHratio(15)];
canelas.mtlD.sumNAD_NADH=1.59;
canelas.mtlD.NADH          = [canelas.mtlD_original.NADH(1)           canelas.mtlD_original.NADH(3)      canelas.mtlD_original.NADH(5)      canelas.mtlD_original.NADH(7)      canelas.mtlD_original.NADH(9)      canelas.mtlD_original.NADH(11)         canelas.mtlD_original.NADH(13)         canelas.mtlD_original.NADH(15)];
canelas.mtlD.NAD          = [canelas.mtlD_original.NAD(1)           canelas.mtlD_original.NAD(3)      canelas.mtlD_original.NAD(5)      canelas.mtlD_original.NAD(7)      canelas.mtlD_original.NAD(9)      canelas.mtlD_original.NAD(11)         canelas.mtlD_original.NAD(13)         canelas.mtlD_original.NAD(15)];

% fluxes
canelas.mtlD.v_HK           = [canelas.mtlD_original.v_HK(1)           canelas.mtlD_original.v_HK(3)      canelas.mtlD_original.v_HK(5)      canelas.mtlD_original.v_HK(7)      canelas.mtlD_original.v_HK(9)      canelas.mtlD_original.v_HK(11)         canelas.mtlD_original.v_HK(13)         canelas.mtlD_original.v_HK(15)];
canelas.mtlD.v_PGI          = [canelas.mtlD_original.v_PGI(1)           canelas.mtlD_original.v_PGI(3)      canelas.mtlD_original.v_PGI(5)      canelas.mtlD_original.v_PGI(7)      canelas.mtlD_original.v_PGI(9)      canelas.mtlD_original.v_PGI(11)         canelas.mtlD_original.v_PGI(13)         canelas.mtlD_original.v_PGI(15)];
canelas.mtlD.v_PFK          = [canelas.mtlD_original.v_PFK(1)           canelas.mtlD_original.v_PFK(3)      canelas.mtlD_original.v_PFK(5)      canelas.mtlD_original.v_PFK(7)      canelas.mtlD_original.v_PFK(9)      canelas.mtlD_original.v_PFK(11)         canelas.mtlD_original.v_PFK(13)         canelas.mtlD_original.v_PFK(15)];
canelas.mtlD.v_FBA          = [canelas.mtlD_original.v_FBA(1)           canelas.mtlD_original.v_FBA(3)      canelas.mtlD_original.v_FBA(5)      canelas.mtlD_original.v_FBA(7)      canelas.mtlD_original.v_FBA(9)      canelas.mtlD_original.v_FBA(11)         canelas.mtlD_original.v_FBA(13)         canelas.mtlD_original.v_FBA(15)];
canelas.mtlD.v_TPI          = [canelas.mtlD_original.v_TPI(1)           canelas.mtlD_original.v_TPI(3)      canelas.mtlD_original.v_TPI(5)      canelas.mtlD_original.v_TPI(7)      canelas.mtlD_original.v_TPI(9)      canelas.mtlD_original.v_TPI(11)         canelas.mtlD_original.v_TPI(13)         canelas.mtlD_original.v_TPI(15)];
canelas.mtlD.v_GAPDH        = [canelas.mtlD_original.v_GAPDH(1)           canelas.mtlD_original.v_GAPDH(3)      canelas.mtlD_original.v_GAPDH(5)      canelas.mtlD_original.v_GAPDH(7)      canelas.mtlD_original.v_GAPDH(9)      canelas.mtlD_original.v_GAPDH(11)         canelas.mtlD_original.v_GAPDH(13)         canelas.mtlD_original.v_GAPDH(15)];
canelas.mtlD.v_PGK          = [canelas.mtlD_original.v_PGK(1)           canelas.mtlD_original.v_PGK(3)      canelas.mtlD_original.v_PGK(5)      canelas.mtlD_original.v_PGK(7)      canelas.mtlD_original.v_PGK(9)      canelas.mtlD_original.v_PGK(11)         canelas.mtlD_original.v_PGK(13)         canelas.mtlD_original.v_PGK(15)];
canelas.mtlD.v_PGM          = [canelas.mtlD_original.v_PGM(1)           canelas.mtlD_original.v_PGM(3)      canelas.mtlD_original.v_PGM(5)      canelas.mtlD_original.v_PGM(7)      canelas.mtlD_original.v_PGM(9)      canelas.mtlD_original.v_PGM(11)         canelas.mtlD_original.v_PGM(13)         canelas.mtlD_original.v_PGM(15)];
canelas.mtlD.v_ENO          = [canelas.mtlD_original.v_ENO(1)           canelas.mtlD_original.v_ENO(3)      canelas.mtlD_original.v_ENO(5)      canelas.mtlD_original.v_ENO(7)      canelas.mtlD_original.v_ENO(9)      canelas.mtlD_original.v_ENO(11)         canelas.mtlD_original.v_ENO(13)         canelas.mtlD_original.v_ENO(15)];
canelas.mtlD.v_PYK          = [canelas.mtlD_original.v_PYK(1)           canelas.mtlD_original.v_PYK(3)      canelas.mtlD_original.v_PYK(5)      canelas.mtlD_original.v_PYK(7)      canelas.mtlD_original.v_PYK(9)      canelas.mtlD_original.v_PYK(11)         canelas.mtlD_original.v_PYK(13)         canelas.mtlD_original.v_PYK(15)];
canelas.mtlD.v_G3PDH        = [canelas.mtlD_original.v_G3PDH(1)           canelas.mtlD_original.v_G3PDH(3)      canelas.mtlD_original.v_G3PDH(5)      canelas.mtlD_original.v_G3PDH(7)      canelas.mtlD_original.v_G3PDH(9)      canelas.mtlD_original.v_G3PDH(11)         canelas.mtlD_original.v_G3PDH(13)         canelas.mtlD_original.v_G3PDH(15)];
canelas.mtlD.v_GLT          = [canelas.mtlD_original.v_GLT(1)           canelas.mtlD_original.v_GLT(3)      canelas.mtlD_original.v_GLT(5)      canelas.mtlD_original.v_GLT(7)      canelas.mtlD_original.v_GLT(9)      canelas.mtlD_original.v_GLT(11)         canelas.mtlD_original.v_GLT(13)         canelas.mtlD_original.v_GLT(15)];
canelas.mtlD.v_PDC          = [canelas.mtlD_original.v_PDC(1)           canelas.mtlD_original.v_PDC(3)      canelas.mtlD_original.v_PDC(5)      canelas.mtlD_original.v_PDC(7)      canelas.mtlD_original.v_PDC(9)      canelas.mtlD_original.v_PDC(11)         canelas.mtlD_original.v_PDC(13)         canelas.mtlD_original.v_PDC(15)];
canelas.mtlD.v_ADH          = [canelas.mtlD_original.v_ADH(1)           canelas.mtlD_original.v_ADH(3)      canelas.mtlD_original.v_ADH(5)      canelas.mtlD_original.v_ADH(7)      canelas.mtlD_original.v_ADH(9)      canelas.mtlD_original.v_ADH(11)         canelas.mtlD_original.v_ADH(13)         canelas.mtlD_original.v_ADH(15)];

% sink reactions
canelas.mtlD.sinkG6P    = canelas.mtlD.v_HK - canelas.mtlD.v_PGI;
canelas.mtlD.sinkF6P    = canelas.mtlD.v_PGI - canelas.mtlD.v_PFK;
canelas.mtlD.sinkGAP    = 2 * canelas.mtlD.v_FBA - canelas.mtlD.v_G3PDH - canelas.mtlD.v_GAPDH;
canelas.mtlD.sinkP3G    = canelas.mtlD.v_PGK - canelas.mtlD.v_PGM;
canelas.mtlD.sinkPEP    = canelas.mtlD.v_ENO - canelas.mtlD.v_PYK;
canelas.mtlD.sinkPYR    = canelas.mtlD.v_PYK - canelas.mtlD.v_PDC;
canelas.mtlD.sinkACE    = canelas.mtlD.v_PDC - canelas.mtlD.v_ADH;


% %%
% figure
% subplot(331), plot(canelas.mtlD.D, canelas.mtlD.sinkG6P, '-o')
% subplot(332), plot(canelas.mtlD.D, canelas.mtlD.sinkF6P, '-o')
% subplot(333), plot(canelas.mtlD.D, canelas.mtlD.sinkGAP, '-o')
% subplot(334), plot(canelas.mtlD.D, canelas.mtlD.sinkP3G, '-o')
% subplot(335), plot(canelas.mtlD.D, canelas.mtlD.sinkPEP, '-o')
% subplot(336), plot(canelas.mtlD.D, canelas.mtlD.sinkPYR, '-o')
% subplot(337), plot(canelas.mtlD.D, canelas.mtlD.sinkACE, '-o')

%% last edit

canelas_SS = canelas; clear canelas


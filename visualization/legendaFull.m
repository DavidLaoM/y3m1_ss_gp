function [legenda] = legendaFull
%LEGENDAFULL produces the legend for the different variables of study and
%writes them down in the variable legenda
%   Detailed explanation goes here

% legenda metabolites           %TO FILL TEMPLEG
legenda.metabolites = {'ACE, x_{1}, [mM]';
    'BPG, x_{2}, [mM]';
    'F16BP, x_{3}, [mM]';
    'F6P, x_{4}, [mM]';
    'G6P, x_{5}, [mM]';
    'GLCi, x_{6}, [mM]';
    'NAD, x_{7}, [mM]';
    'NADH, x_{8}, [mM]';
    'ATP, x_{9}, [mM]';
    'P2G, x_{10}, [mM]';
    'P3G, x_{11}, [mM]';
    'PEP, x_{12}, [mM]';
    'PYR, x_{13}, [mM]';
    'GLYCERAL3P, x_{14}, [mM]';
    'ADP, x_{15}, [mM]';
    'AMP, x_{16}, [mM]';
    'DHAP, x_{17}, [mM]';
    'GLYC3P, x_{18}, [mM]';
    'GLYCEROL, x_{19}, [mM]';
    'ETOH, x_{20}, [mM]';
    'G1P, x_{21}, [mM]';
    'UTP, x_{22}, [mM]';
    'UDP, x_{23}, [mM]';
    'UDP_GLC, x_{24}, [mM]';
    'TRE, x_{25}, [mM]';
    'T6P, x_{26}, [mM]';
    'PI, x_{27}, [mM]';
    'IMP, x_{28}, [mM]';
    'INO, x_{29}, [mM]';
    'HYP, x_{30}, [mM]';
    'ETOH_{EC}, x_{31}, [mM]';
    'GLYC_{EC}, x_{32}, [mM]';
    };

% legena reaction rates
legenda.fluxes = {'v_{GLT}, r_{1}, [mM s^{-1}]';
    'v_{GLK}, r_{2}, [mM s^{-1}]';
    'v_{PGI}, r_{3}, [mM s^{-1}]';
    'v_{PFK}, r_{4}, [mM s^{-1}]';
    'V_{ALD}, r_{5}, [mM s^{-1}]';
    'V_{G3PDH}, r_{6}, [mM s^{-1}]';
    'V_{TPI1}, r_{7}, [mM s^{-1}]';
    'V_{TDH1}, r_{8}, [mM s^{-1}]';
    'V_{PGK}, r_{9}, [mM s^{-1}]';
    'V_{PGM}, r_{10}, [mM s^{-1}]';
    'V_{ENO}, r_{11}, [mM s^{-1}]';
    'V_{PYK}, r_{12}, [mM s^{-1}]';
    'V_{PDC}, r_{13}, [mM s^{-1}]';
    'V_{ADK1}, r_{14}, [mM s^{-1}]';
    'V_{HOR2}, r_{15}, [mM s^{-1}]';
    'V_{RHR2}, r_{16}, [mM s^{-1}]';
    'V_{PGM1}, r_{17}, [mM s^{-1}]';
    'V_{UGP}, r_{18}, [mM s^{-1}]';
    'V_{TPS2}, r_{19}, [mM s^{-1}]';
    'V_{NTH1}, r_{20}, [mM s^{-1}]';
    'V_{TPS1}, r_{21}, [mM s^{-1}]';
    'V_{ACE}, r_{22}, [mM s^{-1}]';
    'V_{ETOHtransport}, r_{23}, [mM s^{-1}]';
    'V_{GLYCEROLtransport}, r_{24}, [mM s^{-1}]';
    'V_{PDH}, r_{25}, [mM s^{-1}]';
    'V_{mitoNADH}, r_{26}, [mM s^{-1}]';
    'V_{ATPase}, r_{27}, [mM s^{-1}]';
    'V_{mito}, r_{28}, [mM s^{-1}]';
    'V_{Amd1}, r_{29}, [mM s^{-1}]';
    'V_{Ade13_v_Ade12}, r_{30}, [mM s^{-1}]';
    'V_{Isn1}, r_{31}, [mM s^{-1}]';
    'V_{Pnp1}, r_{32}, [mM s^{-1}]';
    'V_{Hpt1}, r_{33}, [mM s^{-1}]';
    'V_{sinkG6P}, r_{34}, [mM s^{-1}]';
    'V_{sinkF6P}, r_{35}, [mM s^{-1}]';
    'V_{sinkGAP}, r_{36}, [mM s^{-1}]';
    'V_{sinkP3G}, r_{37}, [mM s^{-1}]';
    'V_{sinkPEP}, r_{38}, [mM s^{-1}]';
    'V_{sinkPYR}, r_{39}, [mM s^{-1}]';
    'V_{sinkACE}, r_{40}, [mM s^{-1}]';
    'V_{ADH}, r_{41}, [mM s^{-1}]';
    'V_{vacuolePi}, r_{42}, [mM s^{-1}]'};

% legenda parameters
legenda.parameters = cell(1,128);

%ADH
legenda.parameters(1) = {'px1, Keq, mM'};
legenda.parameters(2) = {'px2, KiACE, mM'};
legenda.parameters(3) = {'px3, KiETOH, mM'};
legenda.parameters(4) = {'px4, KiNAD, mM'};
legenda.parameters(5) = {'px5, KiNADH, mM'};
legenda.parameters(6) = {'px6, KmACE, mM'};
legenda.parameters(7) = {'px7, KmETOH, mM'};
legenda.parameters(8) = {'px8, KmNAD, mM'};
legenda.parameters(9) = {'px9, KmNADH, mM'};
legenda.parameters(10) = {'px10, Vm, s^{-1}'};

%FBA
legenda.parameters(11) = {'px11, kdhap, mM'};
legenda.parameters(12) = {'px12, keq, []'};
legenda.parameters(13) = {'px13, kf16bp, mM'};
legenda.parameters(14) = {'px14, kgap, mM'};
legenda.parameters(15) = {'px15, kcat, mM s^{-1}'};

%ENO
legenda.parameters(16) = {'px16, k2pg, mM'};
legenda.parameters(17) = {'px17, keq, []'};
legenda.parameters(18) = {'px18, kpep, mM'};
legenda.parameters(19) = {'px19, kcat, mM s^{-1}'};

%TDH1/GAPDH
legenda.parameters(21) = {'px21, keq, []'};
legenda.parameters(22) = {'px22, kgap, mM'};
legenda.parameters(23) = {'px23, kbpg, mM'};
legenda.parameters(24) = {'px24, knad, mM'};
legenda.parameters(25) = {'px25, knadh, mM'};
legenda.parameters(26) = {'px26, kpi, mM'};
legenda.parameters(27) = {'px27, kcat, mM s^{-1}'};

%GLK
legenda.parameters(28) = {'px_{28}, k_{ADP}, mM'};
legenda.parameters(29) = {'px_{29}, k_{ATP}, mM'};
legenda.parameters(30) = {'px_{30}, k_{eq}, []'};
legenda.parameters(31) = {'px_{31}, k_{G6P}, mM'};
legenda.parameters(32) = {'px_{32}, k_{GLC}, mM'};
legenda.parameters(33) = {'px_{33}, k_{T6P}, mM'};
legenda.parameters(34) = {'px_{34}, k_{cat}, s^{-1}'};

%GLT
legenda.parameters(35) = {'px_{35}, k_{GLCoGLCi}, mM'};
legenda.parameters(36) = {'px_{36}, r_{m}, mM'};

%PDC
legenda.parameters(39) = {'px39, kpi, mM'};
legenda.parameters(40) = {'px40, kpyr, mM'};
legenda.parameters(41) = {'px41, hill, []'};
legenda.parameters(42) = {'px42, kcat, mM s^{-1}'};

%PFK
legenda.parameters(43) = {'px_{43}, c_{iAMP}, []'};
legenda.parameters(44) = {'px_{44}, c_{ATP}, []'};
legenda.parameters(45) = {'px_{45}, c_{F16BP}, []'};
legenda.parameters(46) = {'px_{46}, c_{F26BP}, [],'};
legenda.parameters(47) = {'px_{47}, c_{ATP}, []'};
legenda.parameters(48) = {'px_{48}, k_{AMP}, mM'};
legenda.parameters(49) = {'px_{49}, k_{ATP}, mM'};
legenda.parameters(50) = {'px_{50}, k_{F16BP}, mM'};
legenda.parameters(51) = {'px_{51}, k_{F26BP}, mM'};
legenda.parameters(52) = {'px_{52}, k_{F6P}, mM'};
legenda.parameters(53) = {'px_{53}, k_{iATP}, mM'};
legenda.parameters(54) = {'px_{54}, L, []'};
legenda.parameters(55) = {'px_{55}, gR, []'};
legenda.parameters(56) = {'px_{56}, k_{cat}, s^{-1}'};
legenda.parameters(87) = {'px_{87}, F26BP, mM'};

% % PGI
% legenda.parameters(57) = {'px_{57}, k_{eq}, []'};
% legenda.parameters(58) = {'px_{58}, k_{F6P}, mM'};
% legenda.parameters(59) = {'px_{59}, k_{G6P}, mM'};
% legenda.parameters(60) = {'px_{60}, k_{cat}, s^{-1}'};

%TPI
legenda.parameters(57) = {'px57, Keq, []'};
legenda.parameters(58) = {'px58, Kf6p, mM'};
legenda.parameters(59) = {'px59, Kg6p, mM'};
legenda.parameters(60) = {'px60, kcat, mM s^{-1}'};

%PGK
legenda.parameters(61) = {'px61, keq, []'};
legenda.parameters(62) = {'px62, kadp, mM'};
legenda.parameters(63) = {'px63, katp, mM'};
legenda.parameters(64) = {'px64, kbpg, mM'};
legenda.parameters(65) = {'px65, kp3g, mM'};
legenda.parameters(66) = {'px66, rmax, mM s^{-1}'};

%PGM (GPM1)
legenda.parameters(67) = {'px67, k2pg, mM'};
legenda.parameters(68) = {'px68, k3pg, mM'};
legenda.parameters(69) = {'px69, keq, []'};
legenda.parameters(70) = {'px70, kcat, mM s^{-1}'};

%PYK
legenda.parameters(71) = {'px71, kadp, mM'};
legenda.parameters(72) = {'px72, katp, mM'};
legenda.parameters(73) = {'px73, kf16bp, mM'};
legenda.parameters(74) = {'px74, kpep, mM'};
legenda.parameters(75) = {'px75, L, []'};
legenda.parameters(76) = {'px76, hill, []'};
legenda.parameters(77) = {'px77, kcat, mM s^{-1}'};

%TPI
legenda.parameters(79) = {'px79, kdhap, mM'};
legenda.parameters(80) = {'px80, keq, []'};
legenda.parameters(81) = {'px81, kgap, mM'};
legenda.parameters(82) = {'px82, kcat, mM s^{-1}'};

%PGM1
legenda.parameters(83) = {'px83, keq, '};
legenda.parameters(84) = {'px84, kg1p, '};
legenda.parameters(85) = {'px85, kg6p, '};
legenda.parameters(86) = {'px86, kcat, '};

%GPD (G3PDH)
legenda.parameters(96)  = {'px96, kadp, mM'};
legenda.parameters(97)  = {'px97, katp, mM'};
legenda.parameters(98)  = {'px98, kdhap, mM'};
legenda.parameters(99)  = {'px99, keq, []'};
legenda.parameters(100) = {'px100, kf16bp, mM'};
legenda.parameters(101) = {'px101, kglyc3p, mM'};
legenda.parameters(102) = {'px102, knad, mM'};
legenda.parameters(103) = {'px103, knadh, mM'};
legenda.parameters(104) = {'px104, kcat, mM s^{-1}'};

%HOR2 (GlycT)
legenda.parameters(105) = {'px105, kglyc3p, '};
legenda.parameters(106) = {'px106, Kpi, '};
legenda.parameters(107) = {'px107, kcat, '};
legenda.parameters(108) = {'px108, KGlycerolTransport, '};

% % % % % Trehalose cycle % % % % %

%TPS1
legenda.parameters(124) = {'px124, k_{m,G6P}^{TPS1}, mM'};
legenda.parameters(125) = {'px125, k_{m,udpglc}^{TPS1}, mM'};
legenda.parameters(126) = {'px126, k_{cat}^{TPS1}, s^{-1}'};
legenda.parameters(127) = {'px127, k_{m,pi}^{TPS1}, '};
legenda.parameters(128) = {'px128, k_{m,F6P}^{TPS1}, '};
    
%TPS2
legenda.parameters(119) = {'px119, k^{TPS2}_{m,t6p}, mM'};
legenda.parameters(120) = {'px120, k^{TPS2}_{cat}, s^{-1}'};
legenda.parameters(121) = {'px121, k^{TPS2}_{m,pi}, '};
    
%NTH1
legenda.parameters(122) = {'px122, k^{NTH1}_{m,tre}, mM'};
legenda.parameters(123) = {'px123, k^{NTH1}_{cat}, s^{-1}'};

% % % % % Cofactors % % % % %

%mitoNADH
legenda.parameters(92) = {'px92, v^{mitoNADH}_{max}, mM'};
legenda.parameters(93) = {'px93, k^{mitoNADH}_{m,nadh}, mM'};

%vacImport
legenda.parameters(37) = {'px37, k^{vacPi}_{m}, mM^{-1} s^{-1}'};
legenda.parameters(38) = {'px38, Pi^{vacPi}_{SS}, mM'};

%mito: mitochondrial ATP synthesis
legenda.parameters(109) = {'px109, v^{mito}_{max}, '};
legenda.parameters(110) = {'px110, k^{mito}_{m,adp}, '};
legenda.parameters(111) = {'px111, k^{mito}_{m,pi}, '};

%ATPase: cytoplasmic ATPase
legenda.parameters(129) = {'px129, k^{ATPase}_{k}, '};

%ADK: ADK1 consensus
legenda.parameters(130) = {'px130, k^{ADK1}_{keq}, '};
legenda.parameters(131) = {'px131, k^{ADK1}_{k}, '};

% legenda parameter enzyme name
legenda.eName = cell(1,128);

%ADH
legenda.eName(1) = {'ADH'}; %{'px1, Keq, mM'};
legenda.eName(2) = {'ADH'}; %{'px2, KiACE, mM'};
legenda.eName(3) = {'ADH'}; %{'px3, KiETOH, mM'};
legenda.eName(4) = {'ADH'}; %{'px4, KiNAD, mM'};
legenda.eName(5) = {'ADH'}; %{'px5, KiNADH, mM'};
legenda.eName(6) = {'ADH'}; %{'px6, KmACE, mM'};
legenda.eName(7) = {'ADH'}; %{'px7, KmETOH, mM'};
legenda.eName(8) = {'ADH'}; %{'px8, KmNAD, mM'};
legenda.eName(9) = {'ADH'}; %{'px9, KmNADH, mM'};
legenda.eName(10) = {'ADH'}; %{'px10, Vm, s^{-1}'};

%FBA
legenda.eName(11) = {'FBA'}; %{'px11, kdhap, mM'};
legenda.eName(12) = {'FBA'}; %{'px12, keq, []'};
legenda.eName(13) = {'FBA'}; %{'px13, kf16bp, mM'};
legenda.eName(14) = {'FBA'}; %{'px14, kgap, mM'};
legenda.eName(15) = {'FBA'}; %{'px15, kcat, mM s^{-1}'};

%ENO
legenda.eName(16) = {'ENO'}; %{'px16, k2pg, mM'};
legenda.eName(17) = {'ENO'}; %{'px17, keq, []'};
legenda.eName(18) = {'ENO'}; %{'px18, kpep, mM'};
legenda.eName(19) = {'ENO'}; %{'px19, kcat, mM s^{-1}'};

%TDH1/GAPDH
legenda.eName(21) = {'TDH'}; %{'px21, keq, []'};
legenda.eName(22) = {'TDH'}; %{'px22, kgap, mM'};
legenda.eName(23) = {'TDH'}; %{'px23, kbpg, mM'};
legenda.eName(24) = {'TDH'}; %{'px24, knad, mM'};
legenda.eName(25) = {'TDH'}; %{'px25, knadh, mM'};
legenda.eName(26) = {'TDH'}; %{'px26, kpi, mM'};
legenda.eName(27) = {'TDH'}; %{'px27, kcat, mM s^{-1}'};

%GLK
legenda.eName(28) = {'GLK'}; %{'px_{28}, k_{ADP}, mM'};
legenda.eName(29) = {'GLK'}; %{'px_{29}, k_{ATP}, mM'};
legenda.eName(30) = {'GLK'}; %{'px_{30}, k_{eq}, []'};
legenda.eName(31) = {'GLK'}; %{'px_{31}, k_{G6P}, mM'};
legenda.eName(32) = {'GLK'}; %{'px_{32}, k_{GLC}, mM'};
legenda.eName(33) = {'GLK'}; %{'px_{33}, k_{T6P}, mM'};
legenda.eName(34) = {'GLK'}; %{'px_{34}, k_{cat}, s^{-1}'};

%GLT
legenda.eName(35) = {'GLT'}; %{'px_{35}, k_{GLCoGLCi}, mM'};
legenda.eName(36) = {'GLT'}; %{'px_{36}, r_{m}, mM'};

%PDC
legenda.eName(39) = {'PDC'}; %{'px39, kpi, mM'};
legenda.eName(40) = {'PDC'}; %{'px40, kpyr, mM'};
legenda.eName(41) = {'PDC'}; %{'px41, hill, []'};
legenda.eName(42) = {'PDC'}; %{'px42, kcat, mM s^{-1}'};

%PFK
legenda.eName(43) = {'PFK'}; %{'px_{43}, c_{iAMP}, []'};
legenda.eName(44) = {'PFK'}; %{'px_{44}, c_{ATP}, []'};
legenda.eName(45) = {'PFK'}; %{'px_{45}, c_{F16BP}, []'};
legenda.eName(46) = {'PFK'}; %{'px_{46}, c_{F26BP}, [],'};
legenda.eName(47) = {'PFK'}; %{'px_{47}, c_{ATP}, []'};
legenda.eName(48) = {'PFK'}; %{'px_{48}, k_{AMP}, mM'};
legenda.eName(49) = {'PFK'}; %{'px_{49}, k_{ATP}, mM'};
legenda.eName(50) = {'PFK'}; %{'px_{50}, k_{F16BP}, mM'};
legenda.eName(51) = {'PFK'}; %{'px_{51}, k_{F26BP}, mM'};
legenda.eName(52) = {'PFK'}; %{'px_{52}, k_{F6P}, mM'};
legenda.eName(53) = {'PFK'}; %{'px_{53}, k_{iATP}, mM'};
legenda.eName(54) = {'PFK'}; %{'px_{54}, L, []'};
legenda.eName(55) = {'PFK'}; %{'px_{55}, gR, []'};
legenda.eName(56) = {'PFK'}; %{'px_{56}, k_{cat}, s^{-1}'};
legenda.eName(87) = {'PFK'}; %{'px_{87}, F26BP, mM'};

% PGI
legenda.eName(57) = {'PGI'}; %{'px_{57}, k_{eq}, []'};
legenda.eName(58) = {'PGI'}; %{'px_{58}, k_{F6P}, mM'};
legenda.eName(59) = {'PGI'}; %{'px_{59}, k_{G6P}, mM'};
legenda.eName(60) = {'PGI'}; %{'px_{60}, k_{cat}, s^{-1}'};

%PGK
legenda.eName(61) = {'PGK'}; %{'px61, keq, []'};
legenda.eName(62) = {'PGK'}; %{'px62, kadp, mM'};
legenda.eName(63) = {'PGK'}; %{'px63, katp, mM'};
legenda.eName(64) = {'PGK'}; %{'px64, kbpg, mM'};
legenda.eName(65) = {'PGK'}; %{'px65, kp3g, mM'};
legenda.eName(66) = {'PGK'}; %{'px66, rmax, mM s^{-1}'};

%PGM (GPM1)
legenda.eName(67) = {'PGM'}; %{'px67, k2pg, mM'};
legenda.eName(68) = {'PGM'}; %{'px68, k3pg, mM'};
legenda.eName(69) = {'PGM'}; %{'px69, keq, []'};
legenda.eName(70) = {'PGM'}; %{'px70, kcat, mM s^{-1}'};

%PYK
legenda.eName(71) = {'PYK'}; %{'px71, kadp, mM'};
legenda.eName(72) = {'PYK'}; %{'px72, katp, mM'};
legenda.eName(73) = {'PYK'}; %{'px73, kf16bp, mM'};
legenda.eName(74) = {'PYK'}; %{'px74, kpep, mM'};
legenda.eName(75) = {'PYK'}; %{'px75, L, []'};
legenda.eName(76) = {'PYK'}; %{'px76, hill, []'};
legenda.eName(77) = {'PYK'}; %{'px77, kcat, mM s^{-1}'};

%TPI
legenda.eName(79) = {'TPI'}; %{'px79, kdhap, mM'};
legenda.eName(80) = {'TPI'}; %{'px80, keq, []'};
legenda.eName(81) = {'TPI'}; %{'px81, kgap, mM'};
legenda.eName(82) = {'TPI'}; %{'px82, kcat, mM s^{-1}'};

%PGM1
legenda.eName(83) = {'PGM1'}; %{'px83, keq, '};
legenda.eName(84) = {'PGM1'}; %{'px84, kg1p, '};
legenda.eName(85) = {'PGM1'}; %{'px85, kg6p, '};
legenda.eName(86) = {'PGM1'}; %{'px86, kcat, '};

%GPD (G3PDH)
legenda.eName(96) = {'GPD'}; %{'px96, kadp, mM'};
legenda.eName(97) = {'GPD'}; %{'px97, katp, mM'};
legenda.eName(98) = {'GPD'}; %{'px98, kdhap, mM'};
legenda.eName(99) = {'GPD'}; %{'px99, keq, []'};
legenda.eName(100) = {'GPD'}; %{'px100, kf16bp, mM'};
legenda.eName(101) = {'GPD'}; %{'px101, kglyc3p, mM'};
legenda.eName(102) = {'GPD'}; %{'px102, knad, mM'};
legenda.eName(103) = {'GPD'}; %{'px103, knadh, mM'};
legenda.eName(104) = {'GPD'}; %{'px104, kcat, mM s^{-1}'};

%HOR2 (GlycT)
legenda.eName(105) = {'glycT'}; %{'px105, kglyc3p, '};
legenda.eName(106) = {'glycT'}; %{'px106, Kpi, '};
legenda.eName(107) = {'glycT'}; %{'px107, kcat, '};
legenda.eName(108) = {'glycT'}; %{'px108, KGlycerolTransport, '};

% % % % % Trehalose cycle % % % % %

%TPS1
legenda.eName(124) = {'TPS1'}; %{'px124, k_{m,G6P}^{TPS1}, mM'};
legenda.eName(125) = {'TPS1'}; %{'px125, k_{m,udpglc}^{TPS1}, mM'};
legenda.eName(126) = {'TPS1'}; %{'px126, k_{cat}^{TPS1}, s^{-1}'};
legenda.eName(127) = {'TPS1'}; %{'px127, k_{m,pi}^{TPS1}, '};
legenda.eName(128) = {'TPS1'}; %{'px128, k_{m,F6P}^{TPS1}, '};
    
%TPS2
legenda.eName(119) = {'TPS2'}; %{'px119, k^{TPS2}_{m,t6p}, mM'};
legenda.eName(120) = {'TPS2'}; %{'px120, k^{TPS2}_{cat}, s^{-1}'};
legenda.eName(121) = {'TPS2'}; %{'px121, k^{TPS2}_{m,pi}, '};
    
%NTH1
legenda.eName(122) = {'NTH1'}; %{'px122, k^{NTH1}_{m,tre}, mM'};
legenda.eName(123) = {'NTH1'}; %{'px123, k^{NTH1}_{cat}, s^{-1}'};

% % % % % Cofactors % % % % %

%mitoNADH
legenda.eName(92) = {'mitoNADH'}; %{'px92, v^{mitoNADH}_{max}, mM'};
legenda.eName(93) = {'mitoNADH'}; %{'px93, k^{mitoNADH}_{m,nadh}, mM'};

%vacImport
legenda.eName(37) = {'vacImport'}; %{'px37, k^{vacPi}_{m}, mM^{-1} s^{-1}'};
legenda.eName(38) = {'vacImport'}; %{'px38, Pi^{vacPi}_{SS}, mM'};

%mito: mitochondrial ATP synthesis
legenda.eName(109) = {'mito'}; %{'px109, v^{mito}_{max}, '};
legenda.eName(110) = {'mito'}; %{'px110, k^{mito}_{m,adp}, '};
legenda.eName(111) = {'mito'}; %{'px111, k^{mito}_{m,pi}, '};

%ATPase: cytoplasmic ATPase
legenda.eName(129) = {'ATPase'}; %{'px129, k^{ATPase}_{k}, '};

%ADK: ADK1 consensus
legenda.eName(130) = {'ADK'}; %{'px130, k^{ADK1}_{keq}, '};
legenda.eName(131) = {'ADK'}; %{'px131, k^{ADK1}_{k}, '};

end


% Here the names of all parameters, metabolites and fluxes in the
% simulations.

%% Parameters
legenda.parameters = cell(1,128);

%ADH
legenda.parameters(1) = {'x1, Keq, mM'};
legenda.parameters(2) = {'x2, KiACE, mM'};
legenda.parameters(3) = {'x3, KiETOH, mM'};
legenda.parameters(4) = {'x4, KiNAD, mM'};
legenda.parameters(5) = {'x5, KiNADH, mM'};
legenda.parameters(6) = {'x6, KmACE, mM'};
legenda.parameters(7) = {'x7, KmETOH, mM'};
legenda.parameters(8) = {'x8, KmNAD, mM'};
legenda.parameters(9) = {'x9, KmNADH, mM'};
legenda.parameters(10) = {'x10, Vm, s^{-1}'};

%FBA
legenda.parameters(11) = {'x11, kdhap, mM'};
legenda.parameters(12) = {'x12, keq, []'};
legenda.parameters(13) = {'x13, kf16bp, mM'};
legenda.parameters(14) = {'x14, kgap, mM'};
legenda.parameters(15) = {'x15, kcat, mM s^{-1}'};

%ENO
legenda.parameters(16) = {'x16, k2pg, mM'};
legenda.parameters(17) = {'x17, keq, []'};
legenda.parameters(18) = {'x18, kpep, mM'};
legenda.parameters(19) = {'x19, kcat, mM s^{-1}'};

%TDH1/GAPDH
legenda.parameters(21) = {'x21, keq, []'};
legenda.parameters(22) = {'x22, kgap, mM'};
legenda.parameters(23) = {'x23, kbpg, mM'};
legenda.parameters(24) = {'x24, knad, mM'};
legenda.parameters(25) = {'x25, knadh, mM'};
legenda.parameters(26) = {'x26, kpi, mM'};
legenda.parameters(27) = {'x27, kcat, mM s^{-1}'};

%GLK
legenda.parameters(28) = {'x_{28}, k_{ADP}, mM'};
legenda.parameters(29) = {'x_{29}, k_{ATP}, mM'};
legenda.parameters(30) = {'x_{30}, k_{eq}, []'};
legenda.parameters(31) = {'x_{31}, k_{G6P}, mM'};
legenda.parameters(32) = {'x_{32}, k_{GLC}, mM'};
legenda.parameters(33) = {'x_{33}, k_{T6P}, mM'};
legenda.parameters(34) = {'x_{34}, k_{cat}, s^{-1}'};

%GLT
legenda.parameters(35) = {'x_{35}, k_{GLCoGLCi}, mM'};
legenda.parameters(36) = {'x_{36}, v_{m}, mM'};

%PDC
legenda.parameters(39) = {'x39, kpi, mM'};
legenda.parameters(40) = {'x40, kpyr, mM'};
legenda.parameters(41) = {'x41, hill, []'};
legenda.parameters(42) = {'x42, kcat, mM s^{-1}'};

%PFK
legenda.parameters(43) = {'x_{43}, c_{iAMP}, []'};
legenda.parameters(44) = {'x_{44}, c_{ATP}, []'};
legenda.parameters(45) = {'x_{45}, c_{F16BP}, []'};
legenda.parameters(46) = {'x_{46}, c_{F26BP}, [],'};
legenda.parameters(47) = {'x_{47}, c_{ATP}, []'};
legenda.parameters(48) = {'x_{48}, k_{AMP}, mM'};
legenda.parameters(49) = {'x_{49}, k_{ATP}, mM'};
legenda.parameters(50) = {'x_{50}, k_{F16BP}, mM'};
legenda.parameters(51) = {'x_{51}, k_{F26BP}, mM'};
legenda.parameters(52) = {'x_{52}, k_{F6P}, mM'};
legenda.parameters(53) = {'x_{53}, k_{iATP}, mM'};
legenda.parameters(54) = {'x_{54}, L, []'};
legenda.parameters(55) = {'x_{55}, gR, []'};
legenda.parameters(56) = {'x_{56}, k_{cat}, s^{-1}'};
legenda.parameters(87) = {'x_{87}, F26BP, mM'};

% PGI
legenda.parameters(57) = {'x_{57}, k_{eq}, []'};
legenda.parameters(58) = {'x_{58}, k_{F6P}, mM'};
legenda.parameters(59) = {'x_{59}, k_{G6P}, mM'};
legenda.parameters(60) = {'x_{60}, k_{cat}, s^{-1}'};

%TPI
legenda.parameters(57) = {'x57, Keq, []'};
legenda.parameters(58) = {'x58, Kf6p, mM'};
legenda.parameters(59) = {'x59, Kg6p, mM'};
legenda.parameters(60) = {'x60, kcat, mM s^{-1}'};

%PGK
legenda.parameters(61) = {'x61, keq, []'};
legenda.parameters(62) = {'x62, kadp, mM'};
legenda.parameters(63) = {'x63, katp, mM'};
legenda.parameters(64) = {'x64, kbpg, mM'};
legenda.parameters(65) = {'x65, kp3g, mM'};
legenda.parameters(66) = {'x66, vmax, mM s^{-1}'};

%PGM (GPM1)
legenda.parameters(67) = {'x67, k2pg, mM'};
legenda.parameters(68) = {'x68, k3pg, mM'};
legenda.parameters(69) = {'x69, keq, []'};
legenda.parameters(70) = {'x70, kcat, mM s^{-1}'};

%PYK
legenda.parameters(71) = {'x71, kadp, mM'};
legenda.parameters(72) = {'x72, katp, mM'};
legenda.parameters(73) = {'x73, kf16bp, mM'};
legenda.parameters(74) = {'x74, kpep, mM'};
legenda.parameters(75) = {'x75, L, []'};
legenda.parameters(76) = {'x76, hill, []'};
legenda.parameters(77) = {'x77, kcat, mM s^{-1}'};

%TPI
legenda.parameters(79) = {'x79, kdhap, mM'};
legenda.parameters(80) = {'x80, keq, []'};
legenda.parameters(81) = {'x81, kgap, mM'};
legenda.parameters(82) = {'x82, kcat, mM s^{-1}'};

%GPM1
legenda.parameters(83) = {'x83, keq, '};
legenda.parameters(84) = {'x84, kg1p, '};
legenda.parameters(85) = {'x85, kg6p, '};
legenda.parameters(86) = {'x86, kcat, '};

%GPD (G3PDH)
legenda.parameters(96)  = {'x96, kadp, mM'};
legenda.parameters(97)  = {'x97, katp, mM'};
legenda.parameters(98)  = {'x98, kdhap, mM'};
legenda.parameters(99)  = {'x99, keq, []'};
legenda.parameters(100) = {'x100, kf16bp, mM'};
legenda.parameters(101) = {'x101, kglyc3p, mM'};
legenda.parameters(102) = {'x102, knad, mM'};
legenda.parameters(103) = {'x103, knadh, mM'};
legenda.parameters(104) = {'x104, kcat, mM s^{-1}'};

%HOR2 (GlycT)
legenda.parameters(105) = {'x105, kglyc3p, '};
legenda.parameters(106) = {'x106, Kpi, '};
legenda.parameters(107) = {'x107, kcat, '};
legenda.parameters(108) = {'x108, KGlycerolTransport, '};



%% metabolites

legenda.metabolites = {'ACE=IC(:,1), [mM]';
    'BPG=IC(:,2), [mM]';
    'F16BP=IC(:,3), [mM]';
    'F6P=IC(:,4), [mM]';
    'G6P=IC(:,5), [mM]';
    'GLCi=IC(:,6), [mM]';
    'NAD=IC(:,7), [mM]';
    'NADH=IC(:,8), [mM]';
    'ATP=IC(:,9), [mM]';
    'P2G=IC(:,10), [mM]';
    'P3G=IC(:,11), [mM]';
    'PEP=IC(:,12), [mM]';
    'PYR=IC(:,13), [mM]';
    'GLYCERAL3P=IC(:,14), [mM]';
    'ADP=IC(:,15), [mM]';
    'AMP=IC(:,16), [mM]';
    'DHAP=IC(:,17), [mM]';
    'GLYC3P=IC(:,18), [mM]';
    'GLYCEROL=IC(:,19), [mM]';
    'ETOH=IC(:,20), [mM]';
    'G1P=IC(:,21), [mM]';
    'UTP=IC(22), [mM]';
    'UDP=IC(23), [mM]';
    'UDP_GLC=IC(:,24), [mM]';
    'TRE=IC(:,25), [mM]';
    'T6P=IC(:,26), [mM]';
    'PI=IC(:,27), [mM]';
    'IMP=IC(28), [mM]';
    'INO=IC(29), [mM]';
    'HYP=IC(30), [mM]'};

%% Fluxes

legenda.fluxes = {'v_{GLT}, v_{1}, [mM s^{-1}]';
    'v_{GLK}, v_{2}, [mM s^{-1}]';
    'v_{PGI}, v_{3}, [mM s^{-1}]';
    'v_{PFK}, v_{4}, [mM s^{-1}]';
    'V_{ALD}, v_{5}, [mM s^{-1}]';
    'V_{G3PDH}, v_{6}, [mM s^{-1}]';
    'V_{TPI1}, v_{7}, [mM s^{-1}]';
    'V_{TDH1}, v_{8}, [mM s^{-1}]';
    'V_{PGK}, v_{9}, [mM s^{-1}]';
    'V_{PGM}, v_{10}, [mM s^{-1}]';
    'V_{ENO}, v_{11}, [mM s^{-1}]';
    'V_{PYK}, v_{12}, [mM s^{-1}]';
    'V_{PDC}, v_{13}, [mM s^{-1}]';
    'V_{ADK1}, v_{14}, [mM s^{-1}]';
    'V_{HOR2}, v_{15}, [mM s^{-1}]';
    'V_{RHR2}, v_{16}, [mM s^{-1}]';
    'V_{PGM1}, v_{17}, [mM s^{-1}]';
    'V_{UGP}, v_{18}, [mM s^{-1}]';
    'V_{TPS2}, v_{19}, [mM s^{-1}]';
    'V_{NTH1}, v_{20}, [mM s^{-1}]';
    'V_{TPS1}, v_{21}, [mM s^{-1}]';
    'V_{ACE}, v_{22}, [mM s^{-1}]';
    'V_{ETOHtransport}, v_{23}, [mM s^{-1}]';
    'V_{GLYCEROLtransport}, v_{24}, [mM s^{-1}]';
    'V_{PDH}, v_{25}, [mM s^{-1}]';
    'V_{mitoNADH}, v_{26}, [mM s^{-1}]';
    'V_{ATPase}, v_{27}, [mM s^{-1}]';
    'V_{mito}, v_{28}, [mM s^{-1}]';
    'V_{Amd1}, v_{29}, [mM s^{-1}]';
    'V_{Ade13_v_Ade12}, v_{30}, [mM s^{-1}]';
    'V_{Isn1}, v_{31}, [mM s^{-1}]';
    'V_{Pnp1}, v_{32}, [mM s^{-1}]';
    'V_{Hpt1}, v_{33}, [mM s^{-1}]';
    'V_{sinkG6P}, v_{34}, [mM s^{-1}]';
    'V_{sinkF6P}, v_{35}, [mM s^{-1}]';
    'V_{sinkGAP}, v_{36}, [mM s^{-1}]';
    'V_{sinkP3G}, v_{37}, [mM s^{-1}]';
    'V_{sinkPEP}, v_{38}, [mM s^{-1}]';
    'V_{sinkPYR}, v_{39}, [mM s^{-1}]';
    'V_{sinkACE}, v_{40}, [mM s^{-1}]';
    'V_{ADH}, v_{41}, [mM s^{-1}]';
    'V_{vacuolePi}, v_{42}, [mM s^{-1}]'};

%% put everything inside the setup structure
setup.legenda = legenda;
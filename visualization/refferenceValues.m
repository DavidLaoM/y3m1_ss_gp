% referenceVals
% This script selects the reference metabolites and parameter valeus to
% then be easily called in the specific assay.

%% metabolites
InitCond_GP; 
setup.refVals.metabolites   = IC; clear f

%% parameter values
x = zeros(1,128);
setParameterStructure; 
% x = zeros(1,136);
% setup.startParamsOFF = 1;
% setParameterStructure_Y3M1; 

setup.refVals.parameters    = [...

p.ADH.KeqADH .* 10 .^ x(1); %1
p.ADH.KiADHACE .* 10 .^ x(2); %2
p.ADH.KiADHETOH .* 10 .^ x(3); %3
p.ADH.KiADHNAD .* 10 .^ x(4); %4
p.ADH.KiADHNADH .* 10 .^ x(5); %5
p.ADH.KmADHACE .* 10 .^ x(6); %6
p.ADH.KmADHETOH .* 10 .^ x(7); %7
p.ADH.KmADHNAD .* 10 .^ x(8); %8
p.ADH.KmADHNADH .* 10 .^ x(9); %9
p.ADH.VmADH .* 10 .^ x(10); %10

p.FBA1_Kdhap .* 10 .^ x(11); %11
p.FBA1_Keq .* 10 .^ x(12); %12
p.FBA1_Kf16bp .* 10 .^ x(13); %13
p.FBA1_Kglyceral3p .* 10 .^ x(14); %14
p.FBA1_kcat .* 10 .^ x(15); %15
p.ENO1_K2pg .* 10 .^ x(16); %16
p.ENO1_Keq .* 10 .^ x(17); %17
p.ENO1_Kpep .* 10 .^ x(18); %18
p.ENO1_kcat .* 10 .^ x(19); %19
0 .* 10 .^ x(1); %20

p.TDH1_Keq .* 10 .^ x(21); %21
p.TDH1_Kglyceral3p .* 10 .^ x(22); %22
p.TDH1_Kglycerate13bp .* 10 .^ x(23); %23
p.TDH1_Knad .* 10 .^ x(24); %24
p.TDH1_Knadh .* 10 .^ x(25); %25
p.TDH1_Kpi .* 10 .^ x(26); %26
p.TDH1_kcat .* 10 .^ x(27); %27
p.HXK1_Kadp .* 10 .^ x(28); %28
p.HXK1_Katp .* 10 .^ x(29); %29
p.HXK1_Keq .* 10 .^ x(30); %30

p.HXK1_Kg6p .* 10 .^ x(31); %31
p.HXK1_Kglc .* 10 .^ x(32); %32
p.HXK1_Kt6p .* 10 .^ x(33); %33
p.HXK1_kcat .* 10 .^ x(34); %34 
p.GLT.KmGLTGLCi .* 10 .^ x(35); %35 
p.GLT.VmGLT .* 10 .^ x(36); %36
0 .* 10 .^ x(37); %37
0 .* 10 .^ x(38); %38
p.PDC1_Kpi .* 10 .^ x(39); %39
p.PDC1_Kpyr .* 10 .^ x(40); %40

p.PDC1_hill .* 10 .^ x(41); %41
p.PDC1_kcat .* 10 .^ x(42); %42
p.PFK_Camp .* 10 .^ x(43); %43
p.PFK_Catp .* 10 .^ x(44); %44
p.PFK_Cf16bp .* 10 .^ x(45); %45
p.PFK_Cf26bp .* 10 .^ x(46); %46
p.PFK_Ciatp .* 10 .^ x(47); %47
p.PFK_Kamp .* 10 .^ x(48); %48
p.PFK_Katp .* 10 .^ x(49); %49
p.PFK_Kf16bp .* 10 .^ x(50); %50

p.PFK_Kf26bp .* 10 .^ x(51); %51
p.PFK_Kf6p .* 10 .^ x(52); %52
p.PFK_Kiatp .* 10 .^ x(53); %53
p.PFK_L .* 10 .^ x(54); %54
p.PFK_gR .* 10 .^ x(55); %55
p.PFK_kcat .* 10 .^ x(56); %56
p.PGI1_Keq .* 10 .^ x(57); %57
p.PGI1_Kf6p .* 10 .^ x(58); %58
p.PGI1_Kg6p .* 10 .^ x(59); %59
p.PGI1_kcat .* 10 .^ x(60); %60

p.PGK.KeqPGK .* 10 .^ x(61); %61
p.PGK.KmPGKADP .* 10 .^ x(62); %62
p.PGK.KmPGKATP .* 10 .^ x(63); %63
p.PGK.KmPGKBPG .* 10 .^ x(64); %64
p.PGK.KmPGKP3G .* 10 .^ x(65); %65
p.PGK.VmPGK .* 10 .^ x(66); %66
p.GPM1_K2pg .* 10 .^ x(67); %67
p.GPM1_K3pg .* 10 .^ x(68); %68
p.GPM1_Keq .* 10 .^ x(69); %69
p.GPM1_kcat .* 10 .^ x(70); %70

p.PYK1_Kadp .* 10 .^ x(71); %71
p.PYK1_Katp .* 10 .^ x(72); %72
p.PYK1_Kf16bp .* 10 .^ x(73); %73
p.PYK1_Kpep .* 10 .^ x(74); %74
p.PYK1_L .* 10 .^ x(75); %75
p.PYK1_hill .* 10 .^ x(76); %76
p.PYK1_kcat .* 10 .^ x(77); %77
0 .* 10 .^ x(78); %78
p.TPI1_Kdhap .* 10 .^ x(79); %79
p.TPI1_Keq .* 10 .^ x(80); %80

p.TPI1_Kglyceral3p .* 10 .^ x(81); %81
p.TPI1_kcat .* 10 .^ x(82); %82
p.PGM1_Keq .* 10 .^ x(83); %83
p.PGM1_Kg1p .* 10 .^ x(84); %84
p.PGM1_Kg6p .* 10 .^ x(85); %85
p.PGM1_kcat .* 10 .^ x(86); %86
p.PFK.F26BP .* 10 .^ x(87); %87
p.VmaxACE .* 10 .^ x(88); %88
p.PDH_Vmax .* 10 .^ x(89); %89
p.PDH_n .* 10 .^ x(90); %90

p.PDH_K50 .* 10 .^ x(91); %91
p.mitoNADHVmax .* 10 .^ x(92); %92
p.mitoNADHKm .* 10 .^ x(93); %93
p.KmACE .* 10 .^ x(94); %94
p.kETOHtransport .* 10 .^ x(95); %95
p.GPD1_Kadp .* 10 .^ x(96); %96
p.GPD1_Katp .* 10 .^ x(97); %97
p.GPD1_Kdhap .* 10 .^ x(98); %98
p.GPD1_Keq .* 10 .^ x(99); %99
p.GPD1_Kf16bp .* 10 .^ x(100); %100

p.GPD1_Kglyc3p .* 10 .^ x(101); %101
p.GPD1_Knad .* 10 .^ x(102); %102
p.GPD1_Knadh .* 10 .^ x(103); %103
p.GPD1_kcat .* 10 .^ x(104); %104
p.HOR2_Kglyc3p .* 10 .^ x(105); %105
p.HOR2_Kpi .* 10 .^ x(106); %106
p.HOR2_kcat .* 10 .^ x(107); %107
p.GlycerolTransport .* 10 .^ x(108); %108
p.mitoVmax .* 10 .^ x(109); %109
p.mitoADPKm .* 10 .^ x(110); %110

p.mitoPiKm .* 10 .^ x(111); %111
p.Amd1_Vmax .* 10 .^ x(112); %112
p.Amd1_K50 .* 10 .^ x(113); %113
p.Amd1_Kpi .* 10 .^ x(114); %114
p.Ade13_Ade12_k .* 10 .^ x(115); %115
p.Isn1_k .* 10 .^ x(116); %116
p.Pnp1_k .* 10 .^ x(117); %117
p.Hpt1_k .* 10 .^ x(118); %118
p.TPS2_Kt6p .* 10 .^ x(119); %119
p.TPS2_kcat .* 10 .^ x(120); %120

p.TPS2_Kpi .* 10 .^ x(121); %121
p.NTH1_Ktre .* 10 .^ x(122); %122
p.NTH1_kcat .* 10 .^ x(123); %123
p.TPS1_Kg6p .* 10 .^ x(124); %124
p.TPS1_Kudp_glc .* 10 .^ x(125); %125
p.TPS1_kcat .* 10 .^ x(126); %126
p.TPS1_Kpi .* 10 .^ x(127); %127
p.TPS1_KmF6P .* 10 .^ x(128); %128

    ];

clear p x IC f

% parameter number  enzyme name
% 35 to 36          GLT
% 28 to 24          HXK
% 43 to 56 + 87     PFK
% 57 to 60          PGI
% 11 to 15          FBA
% 79 to 82          TPI
% 67 to 70          GPM1 (PGM)
% 16 to 19          ENO
% 71 to 77          PYK
% 21 to 27          TDH (GAPDH)
% 61 to 66          PGK
% 96 to 104         GPD1 (G3PDH)
% 39 to 42          PDC
% 1 to 10           ADH
% 105 to 108        HOR2 + GLYCtransport
% 83 to 86          PGM1
% 88 to 91,94,95    branches (ACE, PDH, ETOHtransport)
% 92 to 93          mitoNADH
% 109 to 111        mitoATP
% 112 to 118        axp metabolism
% 119 to 128        trehalose cycle
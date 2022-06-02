% selectCaseStudy
% select the parameters, metabolites and fluxes data for the case study.

if setup.caseStudy.PGM == 1
    setup.caseStudy.parameters  = (67:70);
    setup.caseStudy.metabolites = [10,11];
    setup.caseStudy.fluxes      = 10;
elseif setup.caseStudy.PGI == 1
    setup.caseStudy.parameters  = (57:60);
    setup.caseStudy.metabolites = [5,4]; %G6P %F6P (reverted order)
    setup.caseStudy.fluxes      = 3;
elseif setup.caseStudy.ENO == 1
    setup.caseStudy.parameters  = (16:19);
    setup.caseStudy.metabolites = [10,12]; %P2G PEP
    setup.caseStudy.fluxes      = 11;
elseif setup.caseStudy.PYK == 1
    setup.caseStudy.parameters  = (71:77);
%     setup.caseStudy.metabolites = [15,12,9,3]; %ADP,PEP,ATP,F16BP %use in case of sPSA
    setup.caseStudy.metabolites = [12,13]; %PEP,PYR
    setup.caseStudy.fluxes      = 12;
elseif setup.caseStudy.PDC == 1
    setup.caseStudy.parameters  = (39:42);
%     setup.caseStudy.metabolites = [13,27]; %PYR, PI %use in case of sPSA
    setup.caseStudy.metabolites = [13,1]; %PYR, ACE
    setup.caseStudy.fluxes      = 13;
elseif setup.caseStudy.GLT == 1
    setup.caseStudy.parameters  = (35:36);
    setup.caseStudy.metabolites = 6; %GLCi sPSA
    setup.caseStudy.fluxes      = 1;
elseif setup.caseStudy.HXK == 1
    setup.caseStudy.parameters  = (28:34);
    setup.caseStudy.metabolites = [9,6,15,5,26]; % ATP,GLCi,ADP,G6P,T6P sPSA
    setup.caseStudy.fluxes      = 2;
elseif setup.caseStudy.GLT_HXK == 1
    setup.caseStudy.parameters  = [(28:36),38];
    setup.caseStudy.metabolites = [9,6,15,5,26]; % ATP,GLCi,ADP,G6P,T6P sPSA
    setup.caseStudy.fluxes      = 1; %wrong at 1, but kept there for functionality of PLA code. %[1,2]; %GLT,HXK
elseif setup.caseStudy.PFK == 1
    setup.caseStudy.parameters  = (43:56);
    setup.caseStudy.metabolites = [4,9,16,3]; % F6P,ATP,AMP,F16BP sPSA
    setup.caseStudy.fluxes      = 4;
elseif setup.caseStudy.FBA == 1
    setup.caseStudy.parameters  = (11:15);
    setup.caseStudy.metabolites = [3,14,17]; % F16BP,GAP,DHAP sPSA
    setup.caseStudy.fluxes      = 5;    
elseif setup.caseStudy.TPI == 1
    setup.caseStudy.parameters  = (79:82);
    setup.caseStudy.metabolites = [17,14]; %DHAP,GAP
    setup.caseStudy.fluxes      = 7;
elseif setup.caseStudy.GPD == 1
    setup.caseStudy.parameters  = (96:104);
    setup.caseStudy.metabolites = [17,8,18,7,3,9,15]; %DHAP,NADH,G3P,NAD,FBP,ATP,ADP
    setup.caseStudy.fluxes      = 6;
elseif setup.caseStudy.TDH1 == 1
    setup.caseStudy.parameters  = (21:27);
    setup.caseStudy.metabolites = [14,7,27,2,8]; %GLYCERAL3P,NAD,PI,GLYCERATE13BP(BPG),NADH
    setup.caseStudy.fluxes      = 8;
elseif setup.caseStudy.PGK == 1
    setup.caseStudy.parameters  = (61:66);
    setup.caseStudy.metabolites = [2,15,9,11]; %BPG,ADP,ATP,P3G
    setup.caseStudy.fluxes      = 9;
elseif setup.caseStudy.GAPDH_PGK == 1
    setup.caseStudy.parameters  = [(21:27),(61:66)];
    setup.caseStudy.metabolites = [14,7,27,2,8,2,15,9,11]; %added from TDH1 and PGK
    setup.caseStudy.fluxes      = 8; %9 also from previous cases
elseif setup.caseStudy.ADH == 1  % % understand the negative sign in the front of its equation!!
    setup.caseStudy.parameters  = (1:10);
    setup.caseStudy.metabolites = [7,20,8,1]; %NAD,ETOH,NADH,ACE
    setup.caseStudy.fluxes      = 41;
elseif setup.caseStudy.HOR2 == 1
    setup.caseStudy.parameters  = (105:107);
    setup.caseStudy.metabolites = [18,27]; %G3P,PI
    setup.caseStudy.fluxes      = 15;
elseif setup.caseStudy.GlycT == 1
    setup.caseStudy.parameters  = (108);
    setup.caseStudy.metabolites = [19]; %glyc
    setup.caseStudy.fluxes      = 24;
% elseif setup.caseStudy.HOR2_glycT == 1
%     setup.caseStudy.parameters  = (105:108);
%     setup.caseStudy.metabolites = [18,19,27]; %G3P,glyc,PI
%     setup.caseStudy.fluxes      = 15;
elseif setup.caseStudy.PGM1 == 1
    setup.caseStudy.parameters  = (83:86);
    setup.caseStudy.metabolites = [21,5]; %G1P,G6P
    setup.caseStudy.fluxes      = 17;
elseif setup.caseStudy.PDC_ADH_branches == 1
    setup.caseStudy.parameters  = [(88:95),(1:10),(39:42),(105:108)]; %boundaries, %ADH, %PDC, %HOR2+glycerol transport
    setup.caseStudy.metabolites = 1; % just a bypass. Not to be run for sensitivity analysis
    setup.caseStudy.fluxes      = 13; % just a bypass. Not to be run for sensitivity analysis
elseif setup.caseStudy.PDC_ADH == 1
    setup.caseStudy.parameters  = [(1:10),(39:42)]; %ADH, %PDC
    setup.caseStudy.metabolites = 1; % just a bypass. Not to be run for sensitivity analysis
    setup.caseStudy.fluxes      = 13; % just a bypass. Not to be run for sensitivity analysis
elseif setup.caseStudy.branches == 1
    setup.caseStudy.parameters  = [(88:95),(105:108)]; %boundaries, %HOR2+glycerol transport
    setup.caseStudy.metabolites = 1; % just a bypass. Not to be run for sensitivity analysis
    setup.caseStudy.fluxes      = 13; % just a bypass. Not to be run for sensitivity analysis
elseif setup.caseStudy.GLT_HXK_old == 1
    setup.caseStudy.parameters  = [(28:36)];
    setup.caseStudy.metabolites = [9,6,15,5,26]; % ATP,GLCi,ADP,G6P,T6P sPSA
    setup.caseStudy.fluxes      = 1; %wrong at 1, but kept there for functionality of PLA code. %[1,2]; %GLT,HXK
% elseif setup.caseStudy.branches == 1
end



% parameter number  enzyme name
% 28 to 34          HXK
% 35 to 36          GLT
% 43 to 56 + 87     PFK
% 57 to 60          PGI
% 11 to 15          FBA
% 79 to 82          TPI
% 67 to 70          GPM1 (PGM)      %
% 16 to 19          ENO             %
% 71 to 77          PYK             %
% 21 to 27          TDH (GAPDH)
% 61 to 66          PGK             %
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
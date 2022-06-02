% initial metaboltie concentrations have been adjusted to the canelas
% dataset, dilution rate 0.1h-1, when possible. Old setup in the end of
% this script, under memoryDump.
% all concentrations, both extracellular and intracellular, are in mM
% (respective volume of bioreactor or cytosol, respectively).
ACE=0.04; % not adjusted to canelas
BPG=0; % not adjusted to canelas
F16BP=canelas_SS.mtlD.FBP(experiment); %0.18;
F6P=canelas_SS.mtlD.F6P(experiment);
G6P=canelas_SS.mtlD.G6P(experiment); %1.7;
GLCi=0.2;
% NAD=1.49;
NAD=canelas_SS.mtlD.NAD(experiment); % Adjusted to the canelas dataset (calculated concentrations) on 2020-02-06
NADH=canelas_SS.mtlD.NADH(experiment);
ATP=canelas_SS.mtlD.ATP(experiment); %4.43./2;
P2G=canelas_SS.mtlD.P2G(experiment); %0.1;
P3G=canelas_SS.mtlD.P3G(experiment); %3.1;
PEP=canelas_SS.mtlD.PEP(experiment); %1.4;
PYR=canelas_SS.mtlD.PYR(experiment); %1;
GLYCERAL3P=canelas_SS.mtlD.GAP(experiment); %0.007; %glyceraldehyde 3-phosphate [cytoplasm] (mmol/L)
ADP=canelas_SS.mtlD.ADP(experiment); %1.06./2;
AMP=canelas_SS.mtlD.AMP(experiment); %0.0693./2;
DHAP=canelas_SS.mtlD.DHAP(experiment); %1; %glycerone phosphate [cytoplasm] (mmol/L)
GLYC3P=canelas_SS.mtlD.G3P(experiment); %0.2; %sn-glycerol 3-phosphate [cytoplasm] (mmol/L)
GLYCEROL=0.1; %glycorol [cytoplasm] (mmol/L)
ETOH=10;%25;
G1P=canelas_SS.mtlD.G1P(experiment); %0.11; %D-glucose 1-phosphate [cytoplasm] (mmol/L)
UTP=0.649; %UTP [cytoplasm] (mmol/L)
UDP=0.281; %UDP [cytoplasm] (mmol/L)
UDP_GLC=0.07; %UPD-D-glucose [cytoplasm] (mmol/L)
% TRE=canelas_SS.mtlD.TRE(experiment); %56; %alpha,alpha-trehalose [cytoplasm] (mmol/L)
TRE=data.metabolites(2,16);
T6P=canelas_SS.mtlD.T6P(experiment); %0.11; %alpha,aplha-trehalose 6-phosphate [cytoplasm] (mmol/L)
PI=10;
IMP=0.1;%0.5; %IMP
INO=0.1;%0.5; %inosine
HYP=1.5; %1 %hypoxyanthine
ETOHec=canelas_SS.mtlD.ECetoh(experiment); %0; % external ethanol product concentration [EC] (mmol/L)
GLYCec=canelas_SS.mtlD.ECglycerol(experiment); %0.0904; % external glycerol by product concentration [EC] (mmol/L)

IC=[ACE;
BPG;
F16BP;
F6P;
G6P;
GLCi;
NAD;
NADH;
ATP;
P2G;
P3G;
PEP;
PYR;
GLYCERAL3P;
ADP;
AMP;
DHAP;
GLYC3P;
GLYCEROL;
ETOH;
G1P
UTP;
UDP;
UDP_GLC;
TRE;
T6P;
PI;
IMP;
INO;
HYP;
ETOHec;
GLYCec];

% clear the variable names
clear ACE BPG F16BP F6P G6P GLCi NAD NADH ATP P2G P3G PEP PYR GLYCERAL3P ADP
clear AMP DHAP GLYC3P GLYCEROL ETOH G1P UTP UDP UDP_GLC TRE T6P PI IMP INO HYP

%fixed states
f.CO2=1;
f.GLCo=0.1;%110;
f.GLY=0.15;
f.Glyc=0;
f.SUCC=0;
f.Trh=0;
f.X=0;
f.GLYCEROL_e=0; %glycerol_out [extracellular] (mmol/L)
f.ETOH_e=0; %ethanol_out [extracellular] (mmol/L)
f.TRE_e=0; %alpha,alpha-trehalose_out [extracellular] (mmol/L)

% %% memoryDump
% ACE=0.04;
% BPG=0;
% F16BP=0.18;
% F6P=0.46;
% G6P=1.7;
% GLCi=0.2;
% % NAD=1.49;
% NAD=1.58; % Adjusted to the canelas dataset (calculated concentrations) on 2020-02-06
% NADH=0.01;
% ATP=4.43./2;
% P2G=0.1;
% P3G=3.1;
% PEP=1.4;
% PYR=1;
% GLYCERAL3P=0.007; %glyceraldehyde 3-phosphate [cytoplasm] (mmol/L)
% ADP=1.06./2;
% AMP=0.0693./2;
% DHAP=1; %glycerone phosphate [cytoplasm] (mmol/L)
% GLYC3P=0.2; %sn-glycerol 3-phosphate [cytoplasm] (mmol/L)
% GLYCEROL=0.1; %glycorol [cytoplasm] (mmol/L)
% ETOH=10;%25;
% G1P=0.11; %D-glucose 1-phosphate [cytoplasm] (mmol/L)
% UTP=0.649; %UTP [cytoplasm] (mmol/L)
% UDP=0.281; %UDP [cytoplasm] (mmol/L)
% UDP_GLC=0.07; %UPD-D-glucose [cytoplasm] (mmol/L)
% TRE=56; %alpha,alpha-trehalose [cytoplasm] (mmol/L)
% T6P=0.11; %alpha,aplha-trehalose 6-phosphate [cytoplasm] (mmol/L)
% PI=10;
% IMP=0.1;%0.5; %IMP
% INO=0.1;%0.5; %inosine
% HYP=1.5; %1 %hypoxyanthine
% ETOHec=0; % external ethanol product concentration [EC] (mmol/L)
% GLYCec=0.0904; % external glycerol by product concentration [EC] (mmol/L)
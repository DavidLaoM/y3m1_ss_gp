ACE=0.04;
BPG=0;
F16BP=0.18;
F6P=0.46;
G6P=1.7;
GLCi=0.2;
NAD=1.49;
NADH=0.01;
ATP=4.43./2;
P2G=0.1;
P3G=3.1;
PEP=1.4;
PYR=1;
GLYCERAL3P=0.007; %glyceraldehyde 3-phosphate [cytoplasm] (mmol/L)
ADP=1.06./2;
AMP=0.0693./2;
DHAP=1; %glycerone phosphate [cytoplasm] (mmol/L)
GLYC3P=0.2; %sn-glycerol 3-phosphate [cytoplasm] (mmol/L)
GLYCEROL=0.1; %glycorol [cytoplasm] (mmol/L)
ETOH=10;%25;
G1P=0.11; %D-glucose 1-phosphate [cytoplasm] (mmol/L)
UTP=0.649; %UTP [cytoplasm] (mmol/L)
UDP=0.281; %UDP [cytoplasm] (mmol/L)
UDP_GLC=0.07; %UPD-D-glucose [cytoplasm] (mmol/L)
TRE=56; %alpha,alpha-trehalose [cytoplasm] (mmol/L)
T6P=0.11; %alpha,aplha-trehalose 6-phosphate [cytoplasm] (mmol/L)
PI=10;
IMP=0.1;%0.5; %IMP
INO=0.1;%0.5; %inosine
HYP=1.5; %1 %hypoxyanthine

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
HYP];

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
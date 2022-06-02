% CLEAN TO GET LITERATURE PARAMETERS HERE. DO ON THE PROCESS


%% GLT
%GLT Teusink
p.GLT.KeqGLT=1;
p.GLT.KmGLTGLCi=10.^x(35).*1.1918;
p.GLT.KmGLTGLCo=10.^x(35).*1.1918;
p.GLT.VmGLT=.2*10.^x(36).*3.67; %--> convert to seconds 

%% GLK 
% HXK1 Concensus
p.HXK1_Kadp=10.^x(28).*0.23; %(UNIT!)
p.HXK1_Katp=10.^x(29).*0.15; %(UNIT!)
p.HXK1_Keq=10.^x(30).*3800; %(UNIT!)
p.HXK1_Kg6p=10.^x(31).*30; %(UNIT!)
p.HXK1_Kglc=10.^x(32).*0.08; %(UNIT!)
p.HXK1_Kt6p=10.^x(33).*0.2; %(UNIT!)
p.HXK1_kcat=.2*10.^x(34).*4.75; %(UNIT!) 

%% ATP + "D-fructose 6-phosphate" = ADP + "D-fructose 1,6-bisphosphate" + H+
% PFK Consensus
p.PFK_Camp=10.^x(43).*0.0845; %(1) %CiAMP
p.PFK_Catp=10.^x(44).*3; %(1)
p.PFK_Cf16bp=10.^x(45).*0.397; %(1) CiF16BP
p.PFK_Cf26bp=10.^x(46).*0.0174; %(1) CiF26BP
p.PFK_Ciatp=10.^x(47).*100; %(1)
p.PFK_Kamp=10.^x(48).*0.095; %(mmol/L)  %corrected value, new value from Teusink model
p.PFK_Katp=10.^x(49).*0.71;%0.71; %(mmol/L) %changed from 0.068 (Manchester kinetics)
p.PFK_Kf16bp=10.^x(50).*0.111; %(mmol/L) %corrected value, new value from Teusink model
p.PFK_Kf26bp=10.^x(51).*0.000682; %(mmol/L) %corrected value, new value from Teusink model
p.PFK_Kf6p=10.^x(52).*0.1; %(mmol/L)  %changed from 0.229 (Manchester kinetics)
p.PFK_Kiatp=10.^x(53).*0.65; %(mmol/L)
p.PFK_L=10.^x(54).*0.66; %(UNIT!)
p.PFK_gR=10.^x(55).*5.12; %(1)
p.PFK_kcat=.2*10.^x(56).*3.55; %(UNIT!)
p.PFK.F26BP=10.^x(87).*0.001; %adjusted from 0.02 (Teusink model), new value based on Canelas data

%% D-glucose 6-phosphate = D-fructose 6-phosphate
% PGI1 consensus
p.PGI1_Keq=10.^x(57).*0.2586; % []
p.PGI1_Kf6p=10.^x(58).*0.307; %mM
p.PGI1_Kg6p=10.^x(59).*1.0257; % mM
p.PGI1_kcat=10.^x(60).*13.4667; %mM/s

%% "D-fructose 1,6-bisphosphate" = "glyceraldehyde 3-phosphate" + "glycerone phosphate";
%FBA1 concensus
p.FBA1_Kdhap=10.^x(11).*2.4; %mM
p.FBA1_Keq=10.^x(12).*0.069; % []
p.FBA1_Kf16bp=10.^x(13).*0.451; % mM
p.FBA1_Kglyceral3p=10.^x(14).*2; % mM
p.FBA1_kcat=10.^x(15).*3.15; % mM/s

%% "glycerone phosphate" = "glyceraldehyde 3-phosphate";  "triose-phosphate isomerase 
%TPI1 consensus
p.TPI1_Kdhap=10.^x(79).*6.45; % mM
p.TPI1_Keq=10.^x(80).*0.0391; %[]
p.TPI1_Kglyceral3p=10.^x(81).*5.25; % mM
p.TPI1_kcat=10.^x(82).*78.396; %mM/s

%% "3-phospho-D-glyceric acid" = "2-phospho-D-glyceric acid"; 
%GPM1 Consensus
% p.GPM1_K2pg=10.^x(67).*0.782; % mM          % Teusink (van Eunen=0.08) (Smallbone=1.4)
% p.GPM1_K3pg=10.^x(68).*1.2; % mM            % (Teusink=0.1) vanEunen = 1.2
% p.GPM1_Keq=10.^x(69).*0.11567; % []         % Teusink & vanEunen 0.19
% p.GPM1_kcat=10.^x(70).*14.2667; % mM/s      % vanEunen=856/60
p.GPM1_K2pg=10.^x(67)   .*0.08;     % mM        % Teusink (van Eunen=0.08) (Smallbone=1.4)
p.GPM1_K3pg=10.^x(68)   .*1.2;      % mM        % (Teusink=0.1) vanEunen = 1.2
p.GPM1_Keq=10.^x(69)    .*0.19;     % []        % Teusink & vanEunen 0.19
p.GPM1_kcat=10.^x(70)   .*14.2667;  % mM/s      % vanEunen=856/60

%% "2-phospho-D-glyceric acid" = phosphoenolpyruvate + water;  
%ENO1 Concensus
p.ENO1_K2pg=10.^x(16).*0.043; % mM
p.ENO1_Keq=4.011*10.^x(17); %[]
p.ENO1_Kpep=10.^x(18).*0.5; % mM
p.ENO1_kcat=10.^x(19).*5.95; % mM/s

%% ADP + H+ + phosphoenolpyruvate -> ATP + pyruvate;  
%PYK1 Consensus
p.PYK1_Kadp=10.^x(71).*0.243; % mM
p.PYK1_Katp=10.^x(72).*9.3; % mM
p.PYK1_Kf16bp=10.^x(73).*0.2; % mM
p.PYK1_Kpep=10.^x(74).*0.281; % mM
p.PYK1_L=10.^x(75).*60000; % []
p.PYK1_hill=10.^x(76).*4; % []
p.PYK1_kcat=10.^x(77).*9.3167; % mM/s

%% "glyceraldehyde 3-phosphate" + NAD+ + phosphate = "3-phospho-D-glyceroyl dihydrogen phosphate" + H+ + NADH;
%TDH1 consensus
p.TDH1_Keq=10.^x(21).*0.0056; % []
p.TDH1_Kglyceral3p=10.^x(22).*0.459; % mM
p.TDH1_Kglycerate13bp=10.^x(23).*0.908; % mM
p.TDH1_Knad=10.^x(24).*2.92; %mM
p.TDH1_Knadh=10.^x(25).*0.022; % mM
p.TDH1_Kpi=10.^x(26).*1.5; % mM
p.TDH1_kcat=10.^x(27).*30.98; % mM/s

%% "3-phospho-D-glyceroyl dihydrogen phosphate" + ADP = "3-phospho-D-glyceric acid" + ATP;  
%PGK Teusink
p.PGK.KeqPGK=10.^x(61).*3200; % []
p.PGK.KmPGKADP=10.^x(62).*0.2; % mM
p.PGK.KmPGKATP=10.^x(63).*0.3; % mM
p.PGK.KmPGKBPG=10.^x(64).*0.003; % mM
p.PGK.KmPGKP3G=10.^x(65).*0.53; % mM
p.PGK.VmPGK=10.^x(66).*44.5; % mM/s

%% "glycerone phosphate" + H+ + NADH = NAD+ + "sn-glycerol 3-phosphate" 
%GPD1 consensus
p.GPD1_Kadp=2*10.^x(96); %(UNIT!)
p.GPD1_Katp=0.73*10.^x(97); %(UNIT!)
p.GPD1_Kdhap=0.54*10.^x(98); %(UNIT!)
p.GPD1_Keq=10000*10.^x(99); %(UNIT!)
p.GPD1_Kf16bp=4.8*10.^x(100); %(UNIT!)
p.GPD1_Kglyc3p=1.2*10.^x(101); %(UNIT!)
p.GPD1_Knad=0.93*10.^x(102); %(UNIT!)
p.GPD1_Knadh=0.023*10.^x(103); %(UNIT!)
p.GPD1_kcat=10.^x(104).*1.169; %(UNIT!)

%% pyruvate + H+ -> acetaldehyde + "carbon dioxide" 
%PDC1 consensus
p.PDC1_Kpi=10.^x(39).*14.7; %(UNIT!)
p.PDC1_Kpyr=10.^x(40).*8.5; %(UNIT!)
p.PDC1_hill=10.^x(41).*1.9; %(UNIT!)
p.PDC1_kcat=10.^x(42).*12.2; %(UNIT!)

%% acetaldehyde + H+ + NADH = ethanol + NAD+;  
%ADH Teusink
p.ADH.KeqADH=10.^x(1).*6.9e-5;
p.ADH.KiADHACE=10.^x(2).*1.1;
p.ADH.KiADHETOH=10.^x(3).*90;
p.ADH.KiADHNAD=10.^x(4).*0.92;
p.ADH.KiADHNADH=10.^x(5).*0.031;
p.ADH.KmADHACE=10.^x(6).*1.11;
p.ADH.KmADHETOH=10.^x(7).*17;
p.ADH.KmADHNAD=10.^x(8).*0.17;
p.ADH.KmADHNADH=10.^x(9).*0.11;
p.ADH.VmADH=10.^x(10).*810./60; %--> convert to seconds

%% "sn-glycerol 3-phosphate" + water -> glycerol + phosphate; 
%HOR2 consensus
p.HOR2_Kglyc3p=5.99.*10.^x(105); %(UNIT!)
p.HOR2_Kpi=1.*10.^x(106); %(UNIT!)
p.HOR2_kcat=10.^x(107).*0.5; %(UNIT!)

%% "D-glucose 1-phosphate" = "D-glucose 6-phosphate";  
%PGM1 consensus
p.PGM1_Keq=10.^x(83).*9.07*3; %(UNIT!)
% % % % p.PGM1_Keq=10.^x(83).*1/6; %(UNIT!) % % % changed update
p.PGM1_Kg1p=10.^x(84).*0.023; %(UNIT!)
p.PGM1_Kg6p=10.^x(85).*0.05; %(UNIT!)
p.PGM1_kcat=10.^x(86).*222; %(UNIT!)
% % % % p.PGM1_kcat=10.^x(86).*100; %(UNIT!) % % % changed update

%% Branches
p.GlycerolTransport=0.1.*10.^x(108);

p.VmaxACE=0.2*10.^x(88);
p.KmACE=0.1*10.^x(94);

p.kETOHtransport=0.1*10.^x(95);

p.PDH_Vmax=0.5*10.^x(89);
p.PDH_n=3*10.^x(90);
p.PDH_K50=0.5*10.^x(91);

p.mitoNADHVmax=1.0*10.^x(92);
p.mitoNADHKm=4*0.030*10.^x(93);

% mitochondrial ATP synthesis
p.mitoVmax=1*10.^x(109); 
p.mitoADPKm=.5*10.^x(110); 
p.mitoPiKm=1*10.^x(111);

% cytoplasmic ATPase
p.ATPaseK=.12;

% AXP metabolism
p.Amd1_Vmax=4*10.^x(112); 
p.Amd1_K50=0.3*10.^x(113);
p.Amd1_Kpi=0.35*10.^x(114);
p.Ade13_Ade12_k=0.05*10.^x(115);
p.Isn1_k=.1*10.^x(116);
p.Pnp1_k=0.03*10.^x(117);
p.Hpt1_k=0.02*10.^x(118);

% vacuole Pi
p.vacuolePi_k=0.0017;
p.vacuolePi_steadyStatePi=10;

%% 2 * ADP = AMP + ATP;  
%ADK1 consensus
p.ADK1_Keq=0.45; %0.275%(UNIT!)
p.ADK1_k=100;% %(UNIT!)

%% trehalose pathway
% "alpha,alpha-trehalose 6-phosphate" + water -> alpha,alpha-trehalose + phosphate;  
p.TPS2_Kt6p=5*10.^x(119); %(UNIT!)
% % % % p.TPS2_Kt6p=0.5*10.^x(119); % % % change update
p.TPS2_kcat=48.7642*10.^x(120); %(UNIT!)
% % % % p.TPS2_kcat=81.45*10.^x(120); % % % change update
p.TPS2_Kpi=0.6*10.^x(121);
% % % % p.TPS2_Kpi=1*10.^x(121); % % % change update

% alpha,alpha-trehalose + water -> 2 * D-glucose;  
p.NTH1_Ktre=2.99*10.^x(122); %(UNIT!)
p.NTH1_kcat=10^(0.54).*.9*10.^x(123); %(UNIT!) 
% % % % p.NTH1_kcat=100*10.^x(123); % % % change update

% "D-glucose 6-phosphate" + UDP-D-glucose -> "alpha,alpha-trehalose 6-phosphate" + H+ + UDP;  
p.TPS1_Kg6p=3.8*10.^x(124); %(UNIT!)
p.TPS1_Kudp_glc=0.886*10.^x(125); %(UNIT!)
p.TPS1_kcat=8.0594e+03*10.^x(126); 
% % % % p.TPS1_kcat=1e+03*10.^x(126); % % % change update
p.TPS1_Kpi=2*1.5*10.^x(127);
% % % % p.TPS1_Kpi=1*10.^x(127); % % % change update
p.TPS1_KmF6P=0.5*1.4*10.^x(128);
% % % % p.TPS1_KmF6P=1*10.^x(128); % % % change update

%% Enzyme concentrations
f.GLK1=0; %glucokinase (GLK1, YCL040W) [cytoplasm] (mmol/L)
f.HXK1=1; %hexokinase 1 (HXK1, YFR053C) [cytoplasm] (mmol/L)
f.HXK2=0; %hexokinase 2 (HXK2, YGL253W) [cytoplasm] (mmol/L)
f.PGI1=1; %glucose-6-phosphate isomerase (PGI1, YBR196C) [cytoplasm] (mmol/L)
f.PFK=1; %PFK1:PFK2 (YGR240C:YMR205C) [cytoplasm] (mmol/L)
f.FBA1=1; %fructose-bisphosphate aldolase (FBA1, YKL060C) [cytoplasm] (mmol/L)
f.GPD1=1; %glycerol-3-phosphate dehydrogenase (GPD1, YDL022W) [cytoplasm] (mmol/L)
f.GPD2=0; %glycerol-3-phosphate dehydrogenase (GPD2, YOL059W) [cytoplasm] (mmol/L)
f.TDH1=1; %glyceraldehyde-3-phosphate dehydrogenase (TDH1, YJL052W) [cytoplasm] (mmol/L)
f.TDH2=0; %glyceraldehyde-3-phosphate dehydrogenase (TDH2, YJR009C) [cytoplasm] (mmol/L)
f.TDH3=0; %glyceraldehyde-3-phosphate dehydrogenase (TDH3, YGR192C) [cytoplasm] (mmol/L)
f.PGK1=0.132; %phosphoglycerate kinase (PGK1, YCR012W) [cytoplasm] (mmol/L)
f.GPM1=1; %phosphoglycerate mutase 1 (GPM1, YKL152C) [cytoplasm] (mmol/L)
f.GPM2=0; %phosphoglycerate mutase 2 (GPM2,  YDL021W) [cytoplasm] (mmol/L)
f.GPM3=0; %phosphoglycerate mutase 3 (GPM3, YOL056W) [cytoplasm] (mmol/L)
f.ENO1=1; %enolase (ENO1, YGR254W) [cytoplasm] (mmol/L)
f.ENO2=0; %enolase (ENO1, YHR174W) [cytoplasm] (mmol/L)
f.PYK1=1; %pyruvate kinase (PYK1, CDC19, YAL038W) [cytoplasm] (mmol/L)
f.PYK2=0; %pyruvate kinase (PYK1, CDC19, YOR347C) [cytoplasm] (mmol/L) 
f.PDC1=0.529; %pyruvate decarboxylase (PDC1, YLR044C) [cytoplasm] (mmol/L) 
f.PDC5=0.00592; %pyruvate decarboxylase (PDC5,  YLR134W) [cytoplasm] (mmol/L) 
f.PDC6=0.00338; %pyruvate decarboxylase (PDC1, YGR087C) [cytoplasm] (mmol/L) 
f.ADH1=0.0933; %alchol dehydrogenase 1 (ADH1, YOL086C) [cytoplasm] (mmol/L) 
f.ADH2=0; %alchol dehydrogenase 2 (ADH2, YMR303C) [cytoplasm] (mmol/L) 
f.ADH3=0.00192; %alchol dehydrogenase 3 (ADH3, YMR083W)) [cytoplasm] (mmol/L) 
f.ADH4=0.0359; %alchol dehydrogenase 4 (ADH4, YGL256W) [cytoplasm] (mmol/L) 
f.ADH5=0.00229; %alchol dehydrogenase 5 (ADH5, YBR145W) [cytoplasm] (mmol/L) 
f.ADH6=0.0171; %alchol dehydrogenase 6 (ADH6, YMR318C) [cytoplasm] (mmol/L) 
f.ADH7=0.0; %alchol dehydrogenase 7 (ADH7, YCR105W) [cytoplasm] (mmol/L) 
f.TPI1=1; %triose-phosphate isomerase (TPI1, YDR050C) [cytoplasm] (mmol/L)
f.HOR2=1; %glycerol-3-phosphatase (HOR2, YER062C) [cytoplasm] (mmol/L)
f.RHR2=0; %glycerol-3-phosphatase  (RHR2, YIL053W)[cytoplasm] (mmol/L)
f.PGM1=1;%0 %phosphoglucomutase 1 (PGM1, YKL127W) [cytoplasm] (mmol/L)
f.PGM2=0; %phosphoglucomutase 2 (PGM2, YMR105C) [cytoplasm] (mmol/L)
f.PGM3=0; %phosphoglucomutase 3 (PGM3, YMR278W) [cytoplasm] (mmol/L)
f.UGP1=0.00031; %UTP-glucose 1-phosphate uridylyltransferase (UGP1, YKL035W) [cytoplasm] (mmol/L)
f.TPS2=0.00133; %alpha,alpha-trehalose phosphatase (TPS2, YDR074W) [cytoplasm] (mmol/L)
f.NTH1=0.00196; %alpha,alpha-trehalase (NTH1, YDR001C) [cytoplasm] (mmol/L)
f.TPS1=0.00145; %TPS1:TPS3:TSL1 (YBR126C:YMR261C:YML100W) [cytoplasm] (mmol/L)

% parameter number  enzyme name
% 35 to 36          GLT
% 28 to 34          HXK
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
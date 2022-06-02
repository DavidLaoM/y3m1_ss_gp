% recall metabolite concentrations
IC=Y();

ACE         =IC(:,1);
BPG         =IC(:,2);
F16BP       =IC(:,3);
F6P         =IC(:,4);
G6P         =IC(:,5);
GLCi        =IC(:,6);
NAD         =IC(:,7);
NADH        =IC(:,8);
ATP         =IC(:,9);
P2G         =IC(:,10);
P3G         =IC(:,11);
PEP         =IC(:,12);
PYR         =IC(:,13);
GLYCERAL3P  =IC(:,14);
ADP         =IC(:,15);
AMP         =IC(:,16);
DHAP        =IC(:,17);
GLYC3P      =IC(:,18);
GLYCEROL    =IC(:,19);
ETOH        =IC(:,20);
G1P         =IC(:,21);
UTP         =0;
UDP         =0;
UDP_GLC=    IC(:,24);
TRE=        IC(:,25);
T6P=        IC(:,26);
PI=         IC(:,27);
IMP=        IC(:,28);
INO=        IC(:,29);
HYP=        IC(:,30);

% calculate reaction rates
rateEquations_Y3M1;

% re-structure results
Q = size(v_GLT);
v = zeros(Q(1),42);
clear Q
    v(:,1)  = v_GLT;
    v(:,2)  = v_GLK;
    v(:,3)  = v_PGI;
    v(:,4)  = v_PFK;
    v(:,5)  = v_ALD;
    v(:,6)  = v_G3PDH;
    v(:,7)  = v_TPI1;
    v(:,8)  = v_GAPDH;
    v(:,9)  = v_PGK;
    v(:,10)  = v_PGM;
    v(:,11)  = v_ENO;
    v(:,12)  = v_PYK;
    v(:,13)  = v_PDC;
    v(:,14)  = v_ADK1;
    v(:,15)  = v_HOR2;
    v(:,16)  = v_RHR2;
    v(:,17)  = -v_PGM1;
    v(:,18)  = v_UGP;
    v(:,19)  = v_TPS2;
    v(:,20)  = v_NTH1;
    v(:,21)  = v_TPS1;
    v(:,22)  = v_ACE;
    v(:,23)  = v_ETOHtransport;
    v(:,24)  = v_GLYCEROLtransport;
    v(:,25)  = v_PDH;
    v(:,26)  = v_mitoNADH;
    v(:,27)  = v_ATPase;
    v(:,28)  = v_mito;
    v(:,29)  = v_Amd1;
    v(:,30)  = v_Ade13_v_Ade12;
    v(:,31)  = v_Isn1;
    v(:,32)  = v_Pnp1;
    v(:,33)  = v_Hpt1;
    v(:,34)  = -v_sinkG6P;
    v(:,35)  = v_sinkF6P;
    v(:,36)  = v_sinkGAP;
    v(:,37)  = -v_sinkP3G;
    v(:,38)  = -v_sinkPEP;
    v(:,39)  = -v_sinkPYR;
    v(:,40)  = -v_sinkACE;
    v(:,41)  = v_ADH;
    v(:,42)  = v_vacuolePi;
    
    
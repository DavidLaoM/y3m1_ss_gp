% commented from function to script
% function [q,Tss,Yss,Vss,Tgs,Ygs,Vgs,canelas,dprofile,...
% exp_G6P,exp_F6P,exp_FBP,exp_GAP,exp_DHAP,exp_G1P,...
% exp_T6P,exp_G3P,exp_P3G,exp_P2G,exp_PEP,exp_PYR,...
% exp_ETOH,exp_AMP,exp_ADP,exp_ATP,exp_NADH,exp_NAD,...
% exp_PI,exp_BPG,exp_GLCi,exp_ACE,exp_GLYC,exp_TRE,...
% exp_IMP,exp_INO,exp_HYP,exp_v_GLT,exp_v_HK,...
% exp_v_PGI,exp_v_PFK,exp_v_FBA,exp_v_TPI,exp_v_G3PDH,...
% exp_v_GAPDH,exp_v_PGK,exp_v_PGM,exp_v_ENO,exp_v_PYK,...
% exp_v_PDC,exp_v_ADH,exp_v_PDH,exp_v_ACE,exp_v_ET_tr,...
% exp_v_HOR2,exp_v_GH_tr,exp_v_PGM1,exp_v_TPS1,exp_v_TPS2,...
% exp_v_NTH1,exp_v_sinkG6P,exp_v_sinkF6P,exp_v_sinkGAP,...
% exp_v_sinkP3G,exp_v_sinkPEP,exp_v_sinkPYR,...
% exp_v_sinkACE,exp_v_mitoNADH,exp_v_mito,...
% exp_v_vacPi,exp_v_ATPase,exp_v_ADK1,exp_v_Amd1,...
% exp_v_Ade1213,exp_v_Isn1,exp_v_Pnp1,exp_v_Hpt1] = ...
%     recallSimsData_Y3M1(simResults, canelas_SS, npSets)

% create emptyCell and fill all the variables that need it
q = 8;
emptyCell = cell(npSets,1);
Tss = emptyCell;
Yss = emptyCell;
Vss = emptyCell;
Tgs = emptyCell;
Ygs = emptyCell;
Vgs = emptyCell;
canelas     = emptyCell;
dprofile    = emptyCell; 
T.UG        = emptyCell;
T.newss     = emptyCell;
Y.newss     = emptyCell;
V.newss     = emptyCell;
%metabolites
exp_G6P     = emptyCell;
exp_F6P     = emptyCell;
exp_FBP     = emptyCell;
exp_GAP     = emptyCell;
exp_DHAP    = emptyCell;
exp_G1P     = emptyCell;
exp_T6P     = emptyCell;
exp_G3P     = emptyCell;
exp_P3G     = emptyCell;
exp_P2G     = emptyCell;
exp_PEP     = emptyCell;
exp_PYR     = emptyCell;
exp_ETOH    = emptyCell;
exp_AMP     = emptyCell;
exp_ADP     = emptyCell;
exp_ATP     = emptyCell;
exp_NADH    = emptyCell;
exp_NAD     = emptyCell;
exp_PI      = emptyCell;
exp_BPG     = emptyCell;
exp_GLCi    = emptyCell;
exp_ACE     = emptyCell;
exp_GLYC    = emptyCell;
exp_TRE     = emptyCell;
exp_IMP     = emptyCell;
exp_INO     = emptyCell;
exp_HYP     = emptyCell;
exp_ECetoh  = emptyCell;
exp_ECglycerol  = emptyCell;
% fluxes
exp_v_GLT       = emptyCell;
exp_v_HK        = emptyCell;
exp_v_PGI       = emptyCell;
exp_v_PFK       = emptyCell;
exp_v_FBA       = emptyCell;
exp_v_TPI       = emptyCell;
exp_v_G3PDH     = emptyCell;
exp_v_GAPDH     = emptyCell;
exp_v_PGK       = emptyCell;
exp_v_PGM       = emptyCell;
exp_v_ENO       = emptyCell;
exp_v_PYK       = emptyCell;
exp_v_PDC       = emptyCell;
exp_v_ADH       = emptyCell;
exp_v_PDH       = emptyCell;
exp_v_ACE       = emptyCell;
exp_v_ET_tr     = emptyCell;
exp_v_HOR2      = emptyCell;
exp_v_GH_tr     = emptyCell;
exp_v_PGM1      = emptyCell;
exp_v_TPS1      = emptyCell;
exp_v_TPS2      = emptyCell;
exp_v_NTH1      = emptyCell;  
exp_v_sinkG6P   = emptyCell;
exp_v_sinkF6P   = emptyCell;
exp_v_sinkGAP   = emptyCell;
exp_v_sinkP3G   = emptyCell;
exp_v_sinkPEP   = emptyCell;
exp_v_sinkPYR   = emptyCell;
exp_v_sinkACE   = emptyCell;
exp_v_mitoNADH  = emptyCell;
exp_v_mito              = emptyCell;
exp_v_vacPi             = emptyCell;
exp_v_ATPase            = emptyCell;
exp_v_ADK1              = emptyCell;
exp_v_Amd1              = emptyCell;
exp_v_Ade1213           = emptyCell;
exp_v_Isn1              = emptyCell;
exp_v_Pnp1              = emptyCell;
exp_v_Hpt1              = emptyCell;

% recall data for each simulation
for j = 1:npSets
    % first recall
    if setup.runGSvanHeerden == 1
        Tgs{j} = simResults{j}.gs.T;
        Ygs{j} = simResults{j}.gs.Y;
        Vgs{j} = simResults{j}.gs.V;
    end
    if setup.runSScanelas == 1
        Tss{j} = simResults{j}.ss.T;
        Yss{j} = simResults{j}.ss.Y;
        Vss{j} = simResults{j}.ss.V;
        canelas{j} = canelas_SS;
        dprofile{j} = canelas_SS.mtlD.D;  
        % canelas data needs to be recalled
        T.UG = Tss{j}.ic2ss; %[0:1:3000]';
        T.newss = Tss{j}.ss2can; 
        Y.newss = Yss{j}.ss2can;
        V.newss = Vss{j}.ss2can;
        %metabolites
        exp_G6P{j}     = zeros(8,1);
        exp_F6P{j}     = zeros(8,1);
        exp_FBP{j}     = zeros(8,1);
        exp_GAP{j}     = zeros(8,1);
        exp_DHAP{j}    = zeros(8,1);
        exp_G1P{j}     = zeros(8,1);
        exp_T6P{j}     = zeros(8,1);
        exp_G3P{j}     = zeros(8,1);
        exp_P3G{j}     = zeros(8,1);
        exp_P2G{j}     = zeros(8,1);
        exp_PEP{j}    = zeros(8,1);
        exp_PYR{j}     = zeros(8,1);
        exp_ETOH{j}    = zeros(8,1);
        exp_AMP{j}     = zeros(8,1);
        exp_ADP{j}     = zeros(8,1);
        exp_ATP{j}     = zeros(8,1);
        exp_NADH{j}    = zeros(8,1);
        exp_NAD{j}     = zeros(8,1);
        exp_PI{j}      = zeros(8,1);
        exp_BPG{j}    = zeros(8,1);
        exp_GLCi{j}    = zeros(8,1);
        exp_ACE{j}     = zeros(8,1);
        exp_GLYC{j}    = zeros(8,1);
        exp_TRE{j}     = zeros(8,1);
        exp_IMP{j}     = zeros(8,1);
        exp_INO{j}     = zeros(8,1);
        exp_HYP{j}     = zeros(8,1);
        exp_ECetoh{j}     = zeros(8,1);
        exp_ECglycerol{j}     = zeros(8,1);
        % fluxes
        exp_v_GLT{j}       = zeros(8,1);
        exp_v_HK{j}        = zeros(8,1);
        exp_v_PGI{j}      = zeros(8,1);
        exp_v_PFK{j}       = zeros(8,1);
        exp_v_FBA{j}       = zeros(8,1);
        exp_v_TPI{j}       = zeros(8,1);
        exp_v_G3PDH{j}     = zeros(8,1);
        exp_v_GAPDH{j}     = zeros(8,1);
        exp_v_PGK{j}       = zeros(8,1);
        exp_v_PGM{j}       = zeros(8,1);
        exp_v_ENO{j}       = zeros(8,1);
        exp_v_PYK{j}       = zeros(8,1);
        exp_v_PDC{j}       = zeros(8,1);
        exp_v_ADH{j}       = zeros(8,1);
        exp_v_PDH{j}       = zeros(8,1);
        exp_v_ACE{j}       = zeros(8,1);
        exp_v_ET_tr{j}     = zeros(8,1);
        exp_v_HOR2{j}      = zeros(8,1);
        exp_v_GH_tr{j}     = zeros(8,1);
        exp_v_PGM1{j}      = zeros(8,1);
        exp_v_TPS1{j}      = zeros(8,1);
        exp_v_TPS2{j}      = zeros(8,1);
        exp_v_NTH1{j}      = zeros(8,1);  
        exp_v_sinkG6P{j}   = zeros(8,1);
        exp_v_sinkF6P{j}   = zeros(8,1);
        exp_v_sinkGAP{j}   = zeros(8,1);
        exp_v_sinkP3G{j}   = zeros(8,1);
        exp_v_sinkPEP{j}   = zeros(8,1);
        exp_v_sinkPYR{j}   = zeros(8,1);
        exp_v_sinkACE{j}   = zeros(8,1);
        exp_v_mitoNADH{j}  = zeros(8,1);
        exp_v_mito{j}              = zeros(8,1);
        exp_v_vacPi{j}             = zeros(8,1);
        exp_v_ATPase{j}            = zeros(8,1);
        exp_v_ADK1{j}              = zeros(8,1);
        exp_v_Amd1{j}              = zeros(8,1);
        exp_v_Ade1213{j}          = zeros(8,1);
        exp_v_Isn1{j}              = zeros(8,1);
        exp_v_Pnp1{j}              = zeros(8,1);
        exp_v_Hpt1{j}              = zeros(8,1);

        % recall the data
        for i = 1:8      
            % Metabolites
            exp_G6P{j}(i)  = Y.newss{i}(end,5);
            exp_F6P{j}(i)   = Y.newss{i}(end,4);
            exp_FBP{j}(i)  = Y.newss{i}(end,3);
            exp_GAP{j}(i)  = Y.newss{i}(end,14);
            exp_DHAP{j}(i) = Y.newss{i}(end,17);
            exp_G1P{j}(i)  = Y.newss{i}(end,21);
            exp_T6P{j}(i)  = Y.newss{i}(end,26);
            exp_G3P{j}(i)  = Y.newss{i}(end,18);
            exp_P3G{j}(i)  = Y.newss{i}(end,11);
            exp_P2G{j}(i)  = Y.newss{i}(end,10);
            exp_PEP{j}(i)  = Y.newss{i}(end,12);
            exp_PYR{j}(i)  = Y.newss{i}(end,13);
            exp_ETOH{j}(i) = Y.newss{i}(end,20);
            exp_AMP{j}(i)  = Y.newss{i}(end,16);
            exp_ADP{j}(i)  = Y.newss{i}(end,15);
            exp_ATP{j}(i)  = Y.newss{i}(end,9);
            exp_NADH{j}(i) = Y.newss{i}(end,8);
            exp_NAD{j}(i)  = Y.newss{i}(end,7);
            exp_PI{j}(i)   = Y.newss{i}(end,27); 
            exp_BPG{j}(i)  = Y.newss{i}(end,2);
            exp_GLCi{j}(i) = Y.newss{i}(end,6);
            exp_ACE{j}(i)  = Y.newss{i}(end,1);
            exp_GLYC{j}(i) = Y.newss{i}(end,19);
            exp_TRE{j}(i)  = Y.newss{i}(end,25);
            exp_IMP{j}(i)  = Y.newss{i}(end,28);
            exp_INO{j}(i)  = Y.newss{i}(end,29);
            exp_HYP{j}(i)  = Y.newss{i}(end,30);
            exp_ECetoh{j}(i)     = Y.newss{i}(end,31);
            exp_ECglycerol{j}(i)     = Y.newss{i}(end,32);
            % Reaction rates
            exp_v_GLT{j}(i)    = V.newss{i}(end,1);
            exp_v_HK{j}(i)     = V.newss{i}(end,2);
            exp_v_PGI{j}(i)    = V.newss{i}(end,3);
            exp_v_PFK{j}(i)    = V.newss{i}(end,4);
            exp_v_FBA{j}(i)    = V.newss{i}(end,5);
            exp_v_TPI{j}(i)    = V.newss{i}(end,7);
            exp_v_G3PDH{j}(i)  = V.newss{i}(end,6);
            exp_v_GAPDH{j}(i)  = V.newss{i}(end,8);
            exp_v_PGK{j}(i)    = V.newss{i}(end,9);
            exp_v_PGM{j}(i)    = V.newss{i}(end,10);
            exp_v_ENO{j}(i)    = V.newss{i}(end,11);
            exp_v_PYK{j}(i)    = V.newss{i}(end,12);
            exp_v_PDC{j}(i)    = V.newss{i}(end,13);
            exp_v_ADH{j}(i)    = V.newss{i}(end,41);
            exp_v_PDH{j}(i)    = V.newss{i}(end,25);
            exp_v_ACE{j}(i)    = V.newss{i}(end,22);
            exp_v_ET_tr{j}(i)  = V.newss{i}(end,23);
            exp_v_HOR2{j}(i)   = V.newss{i}(end,15);
            exp_v_GH_tr{j}(i)  = V.newss{i}(end,24);
            exp_v_PGM1{j}(i)   = V.newss{i}(end,17);
            exp_v_TPS1{j}(i)   = V.newss{i}(end,21);
            exp_v_TPS2{j}(i)   = V.newss{i}(end,21);
            exp_v_NTH1{j}(i)   = V.newss{i}(end,20);
            exp_v_sinkG6P{j}(i)    = V.newss{i}(end,34);
            exp_v_sinkF6P{j}(i)    = V.newss{i}(end,35);
            exp_v_sinkGAP{j}(i)    = V.newss{i}(end,36);
            exp_v_sinkP3G{j}(i)    = V.newss{i}(end,37);
            exp_v_sinkPEP{j}(i)    = V.newss{i}(end,38);
            exp_v_sinkPYR{j}(i)    = V.newss{i}(end,39);
            exp_v_sinkACE{j}(i)    = V.newss{i}(end,40);
            exp_v_mitoNADH{j}(i)   = V.newss{i}(end,26);
            exp_v_mito{j}(i)       = V.newss{i}(end,28);
            exp_v_vacPi{j}(i)      = V.newss{i}(end,42);
            exp_v_ATPase{j}(i)     = V.newss{i}(end,27);
            exp_v_ADK1{j}(i)       = V.newss{i}(end,14);
            exp_v_Amd1{j}(i)       = V.newss{i}(end,29);
            exp_v_Ade1213{j}(i)    = V.newss{i}(end,30);
            exp_v_Isn1{j}(i)       = V.newss{i}(end,31);
            exp_v_Pnp1{j}(i)       = V.newss{i}(end,32);
            exp_v_Hpt1{j}(i)       = V.newss{i}(end,33);
        end
    end
end 
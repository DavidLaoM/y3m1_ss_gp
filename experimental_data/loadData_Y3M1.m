% LOADDATA
% loads canelas and vanHeerden data
CanelasData2; canelas.mtlD = canelas_SS.mtlD_original; %adjustment in hte middle of the versions
    canelas.mtlD.sumAXP     = canelas.mtlD.AMP + canelas.mtlD.ADP + canelas.mtlD.ATP;
    canelas.mtlD.sumNADX    = canelas.mtlD.NAD + canelas.mtlD.NADH;
    % addedData from canelas2011 (supplementary)
%         % growth rate in their units ()
%         d
%         % substrate uptake rate in their units ()
%         d
    canelas_SS.mtlD.TRE = 1/1.7*[230.45 166.36 111.36 18.18 3.18 3.18 3.18 3.18];
    canelas_SS.mtlD.qCO2 = 1/1.7*[0.266990291262136 0.4894822006472497 0.8454692556634313 1.5485436893203888 2.883495145631069 3.221682847896441  3.9959546925566354 5.135113268608415];
    canelas_SS.mtlD.qO2 = 1/1.7*[0.23264401772525822 0.4420236336779908 0.7644017725258492 1.4158050221565732 1.927621861152142 1.8977104874446087 1.3925406203840474 1.1698670605612997];
    canelas_SS.mtlD.qAcet = 1/1.7*[0 0 0 0.0013215822670930022 0.006435574102534064 0.032913440158206936 0.11959147491948519 0.18707008545879839];
    canelas_SS.mtlD.qEtoh = 1/1.7*[0.012376268598816687 0.00036916069126924356 0.006726232934139986 0.007333158138429496 1.0133773823378505 1.3568094505136967 2.436147714332195 3.797531003241106];
    canelas_SS.mtlD.qGlyc = 1/1.7*[0.0003384615384615375 0.00043076923076922763 0.000676923076923075 0.0011692307692307662 0.0018461538461538446 0.0018461538461538446 0.006984615384615384 0.01716923076923077];
    canelas_SS.mtlD.ECgluc = [0.1084337349397595 0.1506024096385543 0.20481927710843362 0.24698795180722932 0.46987951807228967 0.572289156626506 1.4578313253012052 3.289156626506024]; %mM
    canelas_SS.mtlD.ECbiom = [3.578710644677661 3.7188905547226385 3.7683658170914542 3.7518740629685157 3.001499250374813 2.795352323838081 1.8470764617691149 1.4430284857571216]; % g/L
    canelas_SS.mtlD.ECetoh = [0 0 0 0 13.247699142117824 19.69061109765581 41.230783693864815 54.94268569389031]; %mM
    canelas_SS.mtlD.ECco2 = [0.004010490985893234 0.006583376946371615 0.010832136770488061 0.020167861738467893 0.03105182719992366 0.03194615541748488 0.02498775107317895 0.0251528394697251]; %bar
    canelas_SS.mtlD.ECo2 = [0.2716435401441852 0.26846462112911645 0.2634758393485212 0.2548449105457242 0.2510350507697298 0.2518745113983388 0.2620336890378341 0.2636175431725813]; %bar
    canelas_SS.mtlD.ECacet = [0.0914634146341462 0.10518292682926855 0.10060975609756095 0.10060975609756095 1.3170731707317072 1.4725609756097562 2.236280487804878 2.4466463414634143]; %mM
    canelas_SS.mtlD.ECacetaldehyde = [0 0 0 0 0.32472060213091813 0.8354772271658586 2.382203594744902 1.7429385264666992]; %mM
%     canelas_SS.mtlD.ECglycerol = [0.10253719674208356 0.09397662903307935 0.09035910170654188 0.08446888114716039 0.0425323440163195 0.045200893860103675 0.13011775348752214 0.2298592409977443]; %mM
    canelas_SS.mtlD.ECglycerol = canelas_SS.mtlD.GLYCEROL_e;

%         %%
%         % interpolation of the changing dAXPdD (really lousy methof now.
%         % After testing this pproach is useful and valid, best to 
%         % interpolate and then recall like Joep did with dAXPdt. 
%         AXP=sum([canelas_SS.mtlD.ATP; canelas_SS.mtlD.ADP; canelas_SS.mtlD.AMP]);
%         D=canelas_SS.mtlD.ATP;
%         temp=diff(AXP);
%         dAXPdD = [...
%             -temp(1)-temp(2),...
%             -temp(2),...
%             0,... % 0.1 h^{-1} is the reference growth rate.
%             +temp(3),...
%             +temp(3)+temp(4),...
%             +temp(3)+temp(4)+temp(5),...
%             +temp(3)+temp(4)+temp(5)+temp(6),...
%             +temp(3)+temp(4)+temp(5)+temp(6)+temp(7),...
%             ];     
% %         time=[0:1:340];
% %         AXP_intp=interp1(data.time_nucleotides,AXP,time,'linear','extrap');
% % %         AXP_intp=interp1(data.time_nucleotides,AXP,time,'makima','extrap');
% % %         'linear' (default) | 'nearest' | 'next' | 'previous' | 'spline' | 'pchip' | 'cubic' | 'makima'
% %         dAXPdt=diff(AXP_intp);
%         canelas_SS.mtlD.AXP = AXP;
%         canelas_SS.mtlD.D = D;
%         canelas_SS.mtlD.dAXPdD = dAXPdD;
% 
%         clear AXP D temp dAXPdD    
%         %%   
    
    
load('vHeerden_trehalose_data_micromolgdw.mat'); 
    data.nucleotides=data.nucleotides/2;
    data.metabolites=data.metabolites/2;
    data.fluxes=data.fluxes/2;
        AXP=sum(data.nucleotides(:,2:4)'); 
        time=[0:1:340];
        AXP_intp=interp1(data.time_nucleotides,AXP,time,'linear','extrap');
%         AXP_intp=interp1(data.time_nucleotides,AXP,time,'makima','extrap');
%         'linear' (default) | 'nearest' | 'next' | 'previous' | 'spline' | 'pchip' | 'cubic' | 'makima'
        dAXPdt=diff(AXP_intp);
    data.time_dAXPdt=time(2:end);
    data.AXP=AXP;
    data.dAXPdt=dAXPdt;
    clear AXP time AXP_intp dAXPdt
load('TUDdata.mat')
    dataset.FF01.time_mets = [0;5;10;15;20;25;30;60;90;120;150;180;220;250;300;350;400];
    dataset.FF01.timeECgluc= [0;5;11;15;20;   30;60;90;    150;180;220;250;300;350;400];
    dataset.FF03.time_mets = [0;5;10;20;30;40;60;90;120;150;200;250;300;400;550;700;800;900;1000;1200;1400;1600;1700;1803];
    dataset.FF04.time_mets = [0;5;10;15;20;30;60;90;120;150;180;220;250;300;350;398];

% %% visualizing interpolation
% figure,
% subplot(3,3,1)
% plot(data.time_nucleotides, data.nucleotides(:,2),'.')
% title(data.legenda_nucleotides{2})
% subplot(3,3,2)
% plot(data.time_nucleotides, data.nucleotides(:,3),'.')
% title(data.legenda_nucleotides{3})
% subplot(3,3,3)
% plot(data.time_nucleotides, data.nucleotides(:,4),'.')
% title(data.legenda_nucleotides{4})
% subplot(3,3,4)
% plot(data.time_nucleotides, AXP,'.')
% title('AXP')
% subplot(3,3,5)
% plot(time, AXP_intp,'.')
% title('AXP_{interp}')
% subplot(3,3,6)
% plot(data.time_dAXPdt, dAXPdt,'.')
% title('dAXPdt')
% disp('done here')
    
%% restructuring
% to put the data into an array that can be automatically recalled when
% plotting
canelas_SS.ordered.metabolites.onoff = zeros(1,32);
    canelas_SS.ordered.metabolites.onoff(5) = 1;% G6P
    canelas_SS.ordered.metabolites.onoff(12) = 1;% PEP
    canelas_SS.ordered.metabolites.onoff(21) = 1;% G1P
    canelas_SS.ordered.metabolites.onoff(18) = 1;% G3P
    canelas_SS.ordered.metabolites.onoff(13) = 1;% PYR
    canelas_SS.ordered.metabolites.onoff(14) = 1;% GAP
    canelas_SS.ordered.metabolites.onoff(10) = 1;% P2G
    canelas_SS.ordered.metabolites.onoff(17) = 1;% DHAP
    canelas_SS.ordered.metabolites.onoff(11) = 1;% P3G
    canelas_SS.ordered.metabolites.onoff(4) = 1;% F6P
    canelas_SS.ordered.metabolites.onoff(3) = 1;% FBP
    canelas_SS.ordered.metabolites.onoff(26) = 1;% T6P
    canelas_SS.ordered.metabolites.onoff(16) = 1;% AMP
    canelas_SS.ordered.metabolites.onoff(15) = 1;% ADP
    canelas_SS.ordered.metabolites.onoff(9) = 1;% ATP

    canelas_SS.ordered.metabolites.onoff(8) = 1;% NADH
    canelas_SS.ordered.metabolites.onoff(7) = 1;% NAD
    canelas_SS.ordered.metabolites.onoff(25) = 1;% TRE
%     canelas_SS.ordered.metabolites.onoff(32) = 1;% EC_glycerol
%     canelas_SS.ordered.metabolites.onoff(31) = 1;% EC_ETOH
    canelas_SS.ordered.metabolites.onoff(19) = 1;% glycerol
    canelas_SS.ordered.metabolites.onoff(20) = 1;% ETOH
canelas_SS.ordered.metabolites.data = cell(1,32);
    canelas_SS.ordered.metabolites.data{5} = canelas_SS.mtlD.G6P;
    canelas_SS.ordered.metabolites.data{12} = canelas_SS.mtlD.PEP;
    canelas_SS.ordered.metabolites.data{21} = canelas_SS.mtlD.G1P;
    canelas_SS.ordered.metabolites.data{18} = canelas_SS.mtlD.G3P;
    canelas_SS.ordered.metabolites.data{13} = canelas_SS.mtlD.PYR;
    canelas_SS.ordered.metabolites.data{14} = canelas_SS.mtlD.GAP;
    canelas_SS.ordered.metabolites.data{10} = canelas_SS.mtlD.P2G;
    canelas_SS.ordered.metabolites.data{17} = canelas_SS.mtlD.DHAP;
    canelas_SS.ordered.metabolites.data{11} = canelas_SS.mtlD.P3G;
    canelas_SS.ordered.metabolites.data{4} = canelas_SS.mtlD.F6P;
    canelas_SS.ordered.metabolites.data{3} = canelas_SS.mtlD.FBP;
    canelas_SS.ordered.metabolites.data{26} = canelas_SS.mtlD.T6P;
    canelas_SS.ordered.metabolites.data{16} = canelas_SS.mtlD.AMP;
    canelas_SS.ordered.metabolites.data{15} = canelas_SS.mtlD.ADP;
    canelas_SS.ordered.metabolites.data{9} = canelas_SS.mtlD.ATP;

    canelas_SS.ordered.metabolites.data{8} = canelas_SS.mtlD.NADH;
    canelas_SS.ordered.metabolites.data{7} = canelas_SS.mtlD.NAD;
    canelas_SS.ordered.metabolites.data{25} = canelas_SS.mtlD.TRE;
%     canelas_SS.ordered.metabolites.data{32} = canelas_SS.mtlD.ECglycerol;
%     canelas_SS.ordered.metabolites.data{31} = canelas_SS.mtlD.ETOH;
    canelas_SS.ordered.metabolites.data{19} = canelas_SS.mtlD.ECglycerol;
    canelas_SS.ordered.metabolites.data{20} = canelas_SS.mtlD.ETOH;
canelas_SS.ordered.metabolites.std = cell(1,32);
    canelas_SS.ordered.metabolites.std{5} = canelas_SS.mtlD.G6P_sd;
    canelas_SS.ordered.metabolites.std{12} = canelas_SS.mtlD.PEP_sd;
    canelas_SS.ordered.metabolites.std{21} = canelas_SS.mtlD.G1P_sd;
    canelas_SS.ordered.metabolites.std{18} = canelas_SS.mtlD.G3P_sd;
    canelas_SS.ordered.metabolites.std{13} = canelas_SS.mtlD.PYR_sd;
    canelas_SS.ordered.metabolites.std{14} = canelas_SS.mtlD.GAP_sd;
    canelas_SS.ordered.metabolites.std{10} = canelas_SS.mtlD.P2G_sd;
    canelas_SS.ordered.metabolites.std{17} = canelas_SS.mtlD.DHAP_sd;
    canelas_SS.ordered.metabolites.std{11} = canelas_SS.mtlD.P3G_sd;
    canelas_SS.ordered.metabolites.std{4} = canelas_SS.mtlD.F6P_sd;
    canelas_SS.ordered.metabolites.std{3} = canelas_SS.mtlD.FBP_sd;
    canelas_SS.ordered.metabolites.std{26} = canelas_SS.mtlD.T6P_sd;
    canelas_SS.ordered.metabolites.std{16} = canelas_SS.mtlD.AMP_sd;
    canelas_SS.ordered.metabolites.std{15} = canelas_SS.mtlD.ADP_sd;
    canelas_SS.ordered.metabolites.std{9} = canelas_SS.mtlD.ATP_sd;

    canelas_SS.ordered.metabolites.std{8} = zeros(size(canelas_SS.mtlD.NADH));
    canelas_SS.ordered.metabolites.std{7} = zeros(size(canelas_SS.mtlD.NAD));
    canelas_SS.ordered.metabolites.std{25} = zeros(size(canelas_SS.mtlD.TRE));
%     canelas_SS.ordered.metabolites.std{32} = zeros(size(canelas_SS.mtlD.ECglycerol));
%     canelas_SS.ordered.metabolites.std{31} = zeros(size(canelas_SS.mtlD.ETOH));
    canelas_SS.ordered.metabolites.std{19} = zeros(size(canelas_SS.mtlD.ECglycerol));
    canelas_SS.ordered.metabolites.std{20} = zeros(size(canelas_SS.mtlD.ETOH));
canelas_SS.ordered.metabolites.dprof = cell(1,32);
    canelas_SS.ordered.dprof{5} = canelas_SS.mtlD.D;% G6P
    canelas_SS.ordered.dprof{12} = canelas_SS.mtlD.D;% PEP
    canelas_SS.ordered.dprof{21} = canelas_SS.mtlD.D;% G1P
    canelas_SS.ordered.dprof{18} = canelas_SS.mtlD.D;% G3P
    canelas_SS.ordered.dprof{13} = canelas_SS.mtlD.D;% PYR
    canelas_SS.ordered.dprof{14} = canelas_SS.mtlD.D;% GAP
    canelas_SS.ordered.dprof{10} = canelas_SS.mtlD.D;% P2G
    canelas_SS.ordered.dprof{17} = canelas_SS.mtlD.D;% DHAP
    canelas_SS.ordered.dprof{11} = canelas_SS.mtlD.D;% P3G
    canelas_SS.ordered.dprof{4} = canelas_SS.mtlD.D;% F6P
    canelas_SS.ordered.dprof{3} = canelas_SS.mtlD.D;% FBP
    canelas_SS.ordered.dprof{26} = canelas_SS.mtlD.D;% T6P
    canelas_SS.ordered.dprof{16} = canelas_SS.mtlD.D;% AMP
    canelas_SS.ordered.dprof{15} = canelas_SS.mtlD.D;% ADP
    canelas_SS.ordered.dprof{9} = canelas_SS.mtlD.D;% ATP

    canelas_SS.ordered.dprof{8} = canelas_SS.mtlD.D;% NADH
    canelas_SS.ordered.dprof{7} = canelas_SS.mtlD.D;% NAD
    canelas_SS.ordered.dprof{25} = canelas_SS.mtlD.D;% TRE
%     canelas_SS.ordered.dprof{32} = canelas_SS.mtlD.D;% ECglycerol
%     canelas_SS.ordered.dprof{31} = canelas_SS.mtlD.D;% ETOH
    canelas_SS.ordered.dprof{19} = canelas_SS.mtlD.D;% ECglycerol
    canelas_SS.ordered.dprof{20} = canelas_SS.mtlD.D;% ETOH
canelas_SS.ordered.fluxes.onoff = zeros(1,42);
    canelas_SS.ordered.fluxes.onoff(2) = 1;% v_HK
    canelas_SS.ordered.fluxes.onoff(3) = 1;% v_PGI
    canelas_SS.ordered.fluxes.onoff(4) = 1;% v_PFK
    canelas_SS.ordered.fluxes.onoff(5) = 1;% v_FBA
    canelas_SS.ordered.fluxes.onoff(7) = 1;% v_TPI
    canelas_SS.ordered.fluxes.onoff(8) = 1;% v_GAPDH
    canelas_SS.ordered.fluxes.onoff(9) = 1;% v_PGK
    canelas_SS.ordered.fluxes.onoff(10) = 1;% v_PGM
    canelas_SS.ordered.fluxes.onoff(11) = 1;% v_ENO
    canelas_SS.ordered.fluxes.onoff(12) = 1;% v_PYK
    canelas_SS.ordered.fluxes.onoff(6) = 1;% v_G3PDH
    canelas_SS.ordered.fluxes.onoff(1) = 1;% v_GLT
    canelas_SS.ordered.fluxes.onoff(13) = 1;% v_PDC
    canelas_SS.ordered.fluxes.onoff(41) = 1;% v_ADH
    canelas_SS.ordered.fluxes.onoff(34) = 1;% sinkG6P
    canelas_SS.ordered.fluxes.onoff(35) = 1;% sinkF6P
    canelas_SS.ordered.fluxes.onoff(36) = 1;% sinkGAP
    canelas_SS.ordered.fluxes.onoff(37) = 1;% sinkP3G
    canelas_SS.ordered.fluxes.onoff(38) = 1;% sinkPEP
    canelas_SS.ordered.fluxes.onoff(39) = 1;% sinkPYR
    canelas_SS.ordered.fluxes.onoff(40) = 1;% sinkACE
canelas_SS.ordered.fluxes.data = cell(1,42);
    canelas_SS.ordered.fluxes.data{2} = canelas_SS.mtlD.v_HK;
    canelas_SS.ordered.fluxes.data{3} = canelas_SS.mtlD.v_PGI;
    canelas_SS.ordered.fluxes.data{4} = canelas_SS.mtlD.v_PFK;
    canelas_SS.ordered.fluxes.data{5} = canelas_SS.mtlD.v_FBA;
    canelas_SS.ordered.fluxes.data{7} = canelas_SS.mtlD.v_TPI;
    canelas_SS.ordered.fluxes.data{8} = canelas_SS.mtlD.v_GAPDH;
    canelas_SS.ordered.fluxes.data{9} = canelas_SS.mtlD.v_PGK;
    canelas_SS.ordered.fluxes.data{10} = canelas_SS.mtlD.v_PGM;
    canelas_SS.ordered.fluxes.data{11} = canelas_SS.mtlD.v_ENO;
    canelas_SS.ordered.fluxes.data{12} = canelas_SS.mtlD.v_PYK;
    canelas_SS.ordered.fluxes.data{6} = canelas_SS.mtlD.v_G3PDH;
    canelas_SS.ordered.fluxes.data{1} = canelas_SS.mtlD.v_GLT;
    canelas_SS.ordered.fluxes.data{13} = canelas_SS.mtlD.v_PDC;
    canelas_SS.ordered.fluxes.data{41} = canelas_SS.mtlD.v_ADH;
    canelas_SS.ordered.fluxes.data{34} = canelas_SS.mtlD.sinkG6P;
    canelas_SS.ordered.fluxes.data{35} = -canelas_SS.mtlD.sinkF6P;
    canelas_SS.ordered.fluxes.data{36} = -canelas_SS.mtlD.sinkGAP;
    canelas_SS.ordered.fluxes.data{37} = canelas_SS.mtlD.sinkP3G;
    canelas_SS.ordered.fluxes.data{38} = canelas_SS.mtlD.sinkPEP;
    canelas_SS.ordered.fluxes.data{39} = canelas_SS.mtlD.sinkPYR;
    canelas_SS.ordered.fluxes.data{40} = canelas_SS.mtlD.sinkACE;


data.ordered.metabolites.onoff = zeros(1,32);
    data.ordered.metabolites.onoff(5) = 1;%2); %G6P
    data.ordered.metabolites.onoff(4) = 1;%3); %F6P
    data.ordered.metabolites.onoff(3) = 1;%4); %FBP
    data.ordered.metabolites.onoff(14) = 1;%5); %GAP
    data.ordered.metabolites.onoff(11) = 1;%6); %P3G
    data.ordered.metabolites.onoff(12) = 1;%7); %PEP
    data.ordered.metabolites.onoff(21) = 1;%13); %G1P
    data.ordered.metabolites.onoff(24) = 1;%14); %UDPglc
    data.ordered.metabolites.onoff(26) = 1;%15); %T6P
    data.ordered.metabolites.onoff(25) = 1;%16); %TRE
    data.ordered.metabolites.onoff(16) = 1;%2); %AMP
    data.ordered.metabolites.onoff(15) = 1;%3); %ADP
    data.ordered.metabolites.onoff(9) = 1;%4); %ATP
data.ordered.metabolites.data = cell(1,32);
    data.ordered.metabolites.data{5} = data.metabolites(:,2); %G6P
    data.ordered.metabolites.data{4} = data.metabolites(:,3); %F6P
    data.ordered.metabolites.data{3} = data.metabolites(:,4); %FBP
    data.ordered.metabolites.data{14} = data.metabolites(:,5); %GAP
    data.ordered.metabolites.data{11} = data.metabolites(:,6); %P3G
    data.ordered.metabolites.data{12} = data.metabolites(:,7); %PEP
    data.ordered.metabolites.data{21} = data.metabolites(:,13); %G1P
    data.ordered.metabolites.data{24} = data.metabolites(:,14); %UDPglc
    data.ordered.metabolites.data{26} = data.metabolites(:,15); %T6P
    data.ordered.metabolites.data{25} = data.metabolites(:,16); %TRE
    data.ordered.metabolites.data{16} = data.nucleotides(:,2); %AMP
    data.ordered.metabolites.data{15} = data.nucleotides(:,3); %ADP
    data.ordered.metabolites.data{9} = data.nucleotides(:,4); %ATP
data.ordered.metabolites.time = cell(1,32);
    data.ordered.metabolites.time{5} = data.time_metabolites; %G6P
    data.ordered.metabolites.time{4} = data.time_metabolites; %F6P
    data.ordered.metabolites.time{3} = data.time_metabolites; %FBP
    data.ordered.metabolites.time{14} = data.time_metabolites; %GAP
    data.ordered.metabolites.time{11} = data.time_metabolites; %P3G
    data.ordered.metabolites.time{12} = data.time_metabolites; %PEP
    data.ordered.metabolites.time{21} = data.time_metabolites; %G1P
    data.ordered.metabolites.time{24} = data.time_metabolites; %UDPglc
    data.ordered.metabolites.time{26} = data.time_metabolites; %T6P
    data.ordered.metabolites.time{25} = data.time_metabolites; %TRE
    data.ordered.metabolites.time{16} = data.time_nucleotides; %AMP
    data.ordered.metabolites.time{15} = data.time_nucleotides; %ADP
    data.ordered.metabolites.time{9} = data.time_nucleotides; %ATP
data.ordered.fluxes.onoff = zeros(1,42);
    data.ordered.fluxes.onoff(1) = 1; %GLT
    data.ordered.fluxes.onoff(2) = 1; %HXK
    data.ordered.fluxes.onoff(3) = 1; %PGI
    data.ordered.fluxes.onoff(4) = 1; %PFK
    data.ordered.fluxes.onoff(5) = 1; %ALD
    data.ordered.fluxes.onoff(17) = 1; %PGM1
    data.ordered.fluxes.onoff(21) = 1; %TPS1
    data.ordered.fluxes.onoff(19) = 1; %TPS2
    data.ordered.fluxes.onoff(20) = 1; %NTH1
data.ordered.fluxes.data = cell(1,42);
    data.ordered.fluxes.data{1} = data.fluxes(:,1); %GLT
    data.ordered.fluxes.data{2} = data.fluxes(:,1); %HXK
    data.ordered.fluxes.data{3} = data.fluxes(:,2) - data.fluxes(:,3); %PGI
    data.ordered.fluxes.data{4} = data.fluxes(:,4); %PFK
    data.ordered.fluxes.data{5} = data.fluxes(:,5); %ALD
    data.ordered.fluxes.data{17} = data.fluxes(:,6) - data.fluxes(:,7); %PGM1
    data.ordered.fluxes.data{21} = data.fluxes(:,9); %TPS1
    data.ordered.fluxes.data{19} = data.fluxes(:,10); %TPS2
    data.ordered.fluxes.data{20} = data.fluxes(:,11); %NTH1
data.ordered.fluxes.time = cell(1,42);
    data.ordered.fluxes.time{1} = data.time_fluxes; %GLT
    data.ordered.fluxes.time{2} = data.time_fluxes; %HXK
    data.ordered.fluxes.time{3} = data.time_fluxes; %PGI
    data.ordered.fluxes.time{4} = data.time_fluxes; %PFK
    data.ordered.fluxes.time{5} = data.time_fluxes; %ALD
    data.ordered.fluxes.time{17} = data.time_fluxes; %PGM1
    data.ordered.fluxes.time{21} = data.time_fluxes; %TPS1
    data.ordered.fluxes.time{19} = data.time_fluxes; %TPS2
    data.ordered.fluxes.time{20} = data.time_fluxes; %NTH1
   

%% memoryDump
% figure,
% subplot(2,2,1), plot(data.time_nucleotides, AXP), title('AXP')
% subplot(2,2,2), plot(time, AXP_intp), title('AXP_{intp}')
% subplot(2,2,3), plot(data.time_dAXPdt, data.dAXPdt), title('dAXPdt = diff(AXP_{intp})')
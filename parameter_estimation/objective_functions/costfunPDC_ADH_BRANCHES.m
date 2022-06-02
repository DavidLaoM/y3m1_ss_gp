function [error]=costfunPDC_ADH_BRANCHES(x_temp,canelas_SS,setup,x)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
% disp(x_temp);

%adjust parameters
x(setup.caseStudy.parameters) = x_temp;

%set parameter structure
setParameterStructure_Y3M1;

% steady state simulation at D=0.1 (Canelas exp #4)
darray = canelas_SS.mtlD.D;
for experiment=1:8
    
    d=canelas_SS.mtlD.D(experiment);
    PvHoek.D_PDC=[0.047059	0.058824	0.068627	0.086275	0.1	0.115686	0.139216	0.170588	0.201961	0.229412	0.256863	0.276471	0.288235	0.301961	0.309804	0.319608	0.32549	0.333333	0.341176	0.35098	0.358824	0.376471];
    PvHoek.PDC=[0.792899	0.751479	0.704142	0.650888	0.60355	0.568047	0.544379	0.52071	0.514793	0.514793	0.532544	0.579882	0.650888	0.751479	0.840237	0.952663	1.04142	1.147929	1.254438	1.378698	1.443787	1.609467];
    PvHoek.D_ADH=[0.051181	0.076772	0.098425	0.124016	0.149606	0.173228	0.204724	0.222441	0.257874	0.281496	0.311024	0.324803	0.344488	0.375984];
    PvHoek.ADH=[9.881657	9.053254	8.343195	7.514793	6.745562	6.094675	5.088757	4.615385	3.727811	3.254438	2.781065	2.662722	2.95858	4.260355];
    p.PDC_ExprsCor=interp1(PvHoek.D_PDC,PvHoek.PDC,d,'pchip','extrap')./interp1(PvHoek.D_PDC,PvHoek.PDC,0.1,'pchip','extrap');
    p.ADH_ExprsCor=interp1(PvHoek.D_ADH,PvHoek.ADH,d,'pchip','extrap')./interp1(PvHoek.D_ADH,PvHoek.ADH,0.1,'pchip','extrap');

    % set initial conditions (InitCond_PDC_ADH_BRANCHES)
    IC=[canelas_SS.mtlD.PYR(experiment);        %PYR;
        0.1;                                    %ACE;
        canelas_SS.mtlD.ETOH(experiment);       %ETOH;
        canelas_SS.mtlD.G3P(experiment);        %G3P;
        canelas_SS.mtlD.GLYCEROL_e(experiment)];   %GLYCEROL];
    f.GLYCEROL_e=canelas_SS.mtlD.GLYCEROL_e(experiment); %glycerol_out [extracellular] (mmol/L)
    f.ETOH_e=canelas_SS.mtlD.ETOH(experiment); %ethanol_out [extracellular] (mmol/L)
    
    tspan=[0 500];
    options=odeset('RelTol',1e-4,'RelTol',1e-4);
    [T,Y] = ode15s(@ODE_PDC_ADH_BRANCHES,tspan,IC,options,p,f,experiment,canelas_SS,x);

    PYRss(experiment)=Y(end,1);
    ACEss(experiment)=Y(end,2);
    ETOHss(experiment)=Y(end,3);
    GLYC3Pss(experiment)=Y(end,4);
    GLYCEROLss(experiment)=Y(end,5);

    PEP=canelas_SS.mtlD.PEP(experiment);
    ATP=canelas_SS.mtlD.ATP(experiment);
    ADP=canelas_SS.mtlD.ADP(experiment);
    F16BP=canelas_SS.mtlD.FBP(experiment);
    P3G=canelas_SS.mtlD.P3G(experiment);
    DHAP=canelas_SS.mtlD.DHAP(experiment);
    NAD=canelas_SS.mtlD.NAD(experiment);
    NADH=canelas_SS.mtlD.NADH(experiment);
    PI=10;

    PYR=Y(:,1);
    ACE=Y(:,2);
    ETOH=Y(:,3);
    GLYC3P=Y(:,4);
    GLYCEROL=Y(:,5);
    
    rateEquations_PDC_ADH_BRANCHES
        v_ADHss(experiment)=v_ADH(end);
        v_PDCss(experiment)=v_PDC(end);
        v_ACEss(experiment)=v_ACE(end);
        v_PDHss(experiment)=v_PDH(end);
        

%     if setup.parEst.costfun2 == 2
%         figure(5)
%         subplot(2,2,1)
%         plot(T,v_ADH)
%         title('v_{ADH}')
%         subplot(2,2,2)
%         plot(T,v_ETOHtransport)
%         title('v_{ETOHtransport}')
%         subplot(2,2,3)
%         plot(T,ETOH)
%         title('ETOH_{int}')
%         subplot(2,2,4)
%         plot(T,ETOHe*ones(size(ETOH)))
%         title('ETOH_{ect}')
%     end
end


%plotting options
if setup.parEst.drawnowCF == 1
    
    figure(1000)
    
    %G3P
    subplot(2,3,1),plot(GLYC3Pss,canelas_SS.mtlD.G3P,'*')
    title('G3P')
    hold on
    line([0 1.2*max(canelas_SS.mtlD.G3P)],[0 1.2*max(canelas_SS.mtlD.G3P)])
    
    %PYR
    subplot(2,3,2),plot(PYRss,canelas_SS.mtlD.PYR,'*')
    title('PYR')
    hold on
    line([0 1.2*max(canelas_SS.mtlD.PYR)],[0 1.2*max(canelas_SS.mtlD.PYR)])
    
    %vPDH
    subplot(2,3,3),plot(v_PDHss,canelas_SS.mtlD.v_PYK-canelas_SS.mtlD.v_PDC,'*')
    title('vPDH')
    hold on
    line([0 max(canelas_SS.mtlD.v_PYK-canelas_SS.mtlD.v_PDC)*1.2],[0 max(canelas_SS.mtlD.v_PYK-canelas_SS.mtlD.v_PDC)*1.2])
    
    %vPDC
    subplot(2,3,4),plot(v_PDCss,canelas_SS.mtlD.v_PDC,'*')
    title('vPDC')
    hold on
    line([0 max(canelas_SS.mtlD.v_PDC)*1.2],[0 max(canelas_SS.mtlD.v_PDC)*1.2])
    
    %vACE
    subplot(2,3,5),plot(v_ACEss,canelas_SS.mtlD.v_PDC-canelas_SS.mtlD.v_ADH,'*')
    title('vACE')
    hold on
    line([0 max(canelas_SS.mtlD.v_PDC-canelas_SS.mtlD.v_ADH)*1.2],[0 max(canelas_SS.mtlD.v_PDC-canelas_SS.mtlD.v_ADH)*1.2])
    
    %vADH
    subplot(2,3,6),plot(v_ADHss,canelas_SS.mtlD.v_ADH,'*')
    title('vADH')
    hold on
    line([0 max(canelas_SS.mtlD.v_ADH)*1.2],[0 max(canelas_SS.mtlD.v_ADH)*1.2])
    
%     drawnow();
end
% for i=[1,2,5,6,7,8]
% subplot(2,4,i), hold off
% end

if setup.parEst.costfun == 1
    errorPYR    = PYRss-canelas_SS.mtlD.PYR;
    errorv_ADH  = v_ADHss-canelas_SS.mtlD.v_ADH;
    errorv_PDC  = v_PDCss-canelas_SS.mtlD.v_PDC;
    errorv_ACE  = v_ACEss-(canelas_SS.mtlD.v_PDC-canelas_SS.mtlD.v_ADH);
    errorv_PDH  = v_PDHss-(canelas_SS.mtlD.v_PYK-canelas_SS.mtlD.v_PDC);
    errorG3P    = GLYC3Pss-canelas_SS.mtlD.G3P;

%     errorPYR    =(PYRss-canelas_SS.mtlD.PYR)./canelas_SS.mtlD.PYR;
%     errorv_ADH  =(v_ADHss-canelas_SS.mtlD.v_ADH)./canelas_SS.mtlD.v_ADH);
% %     errorv_ADH1 =(v_ADHss(10:end)-canelas_SS.mtlD.v_ADH(10:end))./canelas_SS.mtlD.v_ADH(10:end);
% %     errorv_ADH2 =v_ADHss(1:9)-canelas_SS.mtlD.v_ADH(1:9);
%     errorv_PDC  =(v_PDCss-canelas_SS.mtlD.v_PDC)./canelas_SS.mtlD.v_PDC;
%     errorv_ACE  =(v_ACEss-(canelas_SS.mtlD.v_PDC-canelas_SS.mtlD.v_ADH))./v_ACEss;
%     errorv_PDH  =(v_PDHss-(canelas_SS.mtlD.v_PYK-canelas_SS.mtlD.v_PDC))./(canelas_SS.mtlD.v_PYK-canelas_SS.mtlD.v_PDC);
%     errorG3P    =(GLYC3Pss-canelas_SS.mtlD.G3P)./canelas_SS.mtlD.G3P;

    lambda      = setup.parEst.lambda;
    error       =[errorPYR';
                   errorv_ADH';
%                    errorv_ADH1';
%                    errorv_ADH2';
                   errorv_PDC';
                   errorv_ACE';
                   errorv_PDH';
                   errorG3P';
                   lambda.*x_temp'];

end
end
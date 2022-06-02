function [error]=costfunGAPDH_PGK(x_temp,canelas_SS,setup,x)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
% disp(x_temp);

%adjust parameters
x(setup.caseStudy.parameters) = x_temp;

%set parameter structure
setParameterStructure_Y3M1;

%select experimental data
v_GAPDH_data    = canelas_SS.mtlD.v_GAPDH;
v_PGK_data      = canelas_SS.mtlD.v_PGK;
GLYCERAL3P_data = canelas_SS.mtlD.GAP;
P3G_data        = canelas_SS.mtlD.P3G;
darray = canelas_SS.mtlD.D;

%simulated data
if setup.parEst.costfun == 1
    for experiment=1:8
        IC(1)=0.001; % initial condution BPG

        d=canelas_SS.mtlD.D(experiment);
        PvHoek.D_GAPDH=[0.051867	0.072614	0.09751	0.114108	0.134855	0.155602	0.178423	0.19917	0.221992	0.248963	0.275934	0.29668	0.319502	0.348548	0.375519];
        PvHoek.GAPDH=[5.981308	5.358255	4.672897	4.23676	3.862928	3.489097	3.239875	3.05296	2.990654	3.115265	3.613707	3.925234	4.423676	5.109034	6.292835];
        PvHoek.D_PGK=[0.052523	0.072215	0.095844	0.119474	0.135221	0.154929	0.178566	0.202218	0.229816	0.255468	0.287065	0.306811	0.32855	0.34636	0.364177	0.382034];
        PvHoek.PGK=[7.979653	7.444776	6.791089	6.137403	5.662165	5.245629	4.651114	4.17494	3.639128	3.340233	3.158978	3.038297	3.035725	3.211133	3.445711	3.976145];
        p.GAPDH_ExprsCor=interp1(PvHoek.D_GAPDH,PvHoek.GAPDH,d,'pchip','extrap')./interp1(PvHoek.D_GAPDH,PvHoek.GAPDH,0.1,'pchip','extrap');
        p.PGK_ExprsCor=interp1(PvHoek.D_PGK,PvHoek.PGK,d,'pchip','extrap')./interp1(PvHoek.D_PGK,PvHoek.PGK,0.1,'pchip','extrap');

        options=[];
        tspan=[0 500];
        options=[];
        [T,Y] = ode15s(@subsystemGAPDH_PGK,tspan,IC,options,p,f,experiment,canelas_SS,setup);

        BPG=Y(end,1);
        ATP=canelas_SS.mtlD.ATP(experiment);
        ADP=canelas_SS.mtlD.ADP(experiment);
        NADH=canelas_SS.mtlD.NADH(experiment);
        NAD=canelas_SS.mtlD.NAD(experiment);
        GLYCERAL3P=canelas_SS.mtlD.GAP(experiment);
        P3G=canelas_SS.mtlD.P3G(experiment);
        PI=10;

        GLYCERATE13BP=BPG;
        v_TDH1=p.GAPDH_ExprsCor.*((((p.TDH1_kcat.*(f.TDH1+f.TDH2+f.TDH3))./(p.TDH1_Kglyceral3p.*p.TDH1_Knad.*p.TDH1_Kpi)).*(GLYCERAL3P.*NAD.*PI-(GLYCERATE13BP.*NADH)./p.TDH1_Keq))./...
        ((1+GLYCERAL3P./p.TDH1_Kglyceral3p).*(1+NAD./p.TDH1_Knad).*(1+PI./p.TDH1_Kpi)+(1+GLYCERATE13BP./p.TDH1_Kglycerate13bp).*(1+NADH./p.TDH1_Knadh)-1));
        v_GAPDH_predicted(experiment)=v_TDH1;

        v_PGK_predicted(experiment)=p.PGK_ExprsCor.*p.PGK.VmPGK.*((p.PGK.KeqPGK.*BPG.*ADP)-ATP.*P3G)./...
         (p.PGK.KmPGKATP*p.PGK.KmPGKP3G.*(1+ADP./p.PGK.KmPGKADP + ATP./p.PGK.KmPGKATP).*(1+BPG/p.PGK.KmPGKBPG+P3G/p.PGK.KmPGKP3G));

    end
end
 
%plotting options
if setup.parEst.drawnowCF == 1

    figure(1000)
    
    subplot(2,2,1),
    plot(v_GAPDH_predicted,v_GAPDH_data,'r*')
    line([0 max(v_GAPDH_data)*1.2],[0 max(v_GAPDH_data)*1.2])
    xlabel('predicted flux (mM/s)')
    ylabel('measured flux (mM/s)')
    title('GAPDH')
    
    subplot(2,2,3),
    plot(v_PGK_predicted,v_PGK_data,'r*')
    line([0 max(v_PGK_data)*1.2],[0 max(v_PGK_data)*1.2])
    xlabel('predicted flux (mM/s)')
    ylabel('measured flux (mM/s)')
    title('PGK')

    subplot(2,2,2) 
    bar(x_temp)
    title('parameter estimates')
    
    delete(findall(gcf,'type','annotation'))
    dim = [0.65 .35 0 .1];
    
    str1 = ['  p_{1}=  ', num2str(x_temp(1)),'  p_{2}=  ', num2str(x_temp(2))];
    str2 = ['  p_{3}=  ', num2str(x_temp(3)),'  p_{4}=  ', num2str(x_temp(4))];
    str3 = ['  p_{5}=  ', num2str(x_temp(5)),'  p_{6}=  ', num2str(x_temp(6))];
    str4 = ['  p_{7}=  ', num2str(x_temp(7)),'  p_{8}=  ', num2str(x_temp(8))];
    str5 = ['  p_{9}=  ', num2str(x_temp(9)),'  p_{10}=  ', num2str(x_temp(10))];
    str6 = ['  p_{11}=  ', num2str(x_temp(11)),'  p_{12}=  ', num2str(x_temp(12))];
    str7 = ['  p_{13}=  ', num2str(x_temp(13))];
    str = {str1, str2, str3, str4, str5, str6, str7};
    annotation('textbox',dim,'string',str,'FitBoxToText','on')
 

    drawnow()

end

% different types of cost function
if setup.parEst.costfun == 1
    lambda      = setup.parEst.lambda;
    errorGAPDH  = v_GAPDH_data-v_GAPDH_predicted;
    errorPGK    = v_PGK_data-v_PGK_predicted;
    error       = [errorGAPDH,errorPGK,lambda.*x_temp(1:13)];
end

end
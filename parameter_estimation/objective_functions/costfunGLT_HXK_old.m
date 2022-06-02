function [error]=costfunGLT_HXK(x_temp,canelas_SS,setup,x)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
% disp(x_temp);

%adjust parameters
x(setup.caseStudy.parameters) = x_temp;

%set parameter structure
setParameterStructure;

%select experimental data
v_GLT_data  = canelas_SS.mtlD.v_GLT;
v_HXK_data  = canelas_SS.mtlD.v_HK;
v_GLK_data = v_HXK_data;
G6P_data    = canelas_SS.mtlD.G6P;
darray = canelas_SS.mtlD.D;

%simulated data
for experiment=1:8
    IC(1)=0.1; %GLCi

    options=[];
    tspan=[0 500];
    options=[];

    [T,Y] = ode15s(@subsystemGLT_HXK,tspan,IC,options,p,f,experiment,canelas_SS,setup);

    ATP     = canelas_SS.mtlD.ATP(experiment);
    ADP     = canelas_SS.mtlD.ADP(experiment);
    T6P     = canelas_SS.mtlD.T6P(experiment);
    G6P     = canelas_SS.mtlD.G6P(experiment);
    f.GLCo  = canelas_SS.mtlD.Cs(experiment);
    GLCi    = Y(end,1);
    d       = canelas_SS.mtlD.D(experiment);
    PvHoek.D_HXK    = [0.048535	0.065089	0.079577	0.094065	0.112703	0.135481	0.154118	0.176879	0.199657	0.222436	0.245214	0.270084	0.303244	0.332263	0.34886	0.375839];
    PvHoek.HXK      = [2.612394 2.455598 2.330085 2.204571 2.078798 1.921615 1.795843 1.576355 1.419172 1.261989 1.104806 1.009798 0.883121 0.787855 0.78682 0.816292];
    p.HXK_ExprsCor  = interp1(PvHoek.D_HXK,PvHoek.HXK,d,'pchip','extrap')./interp1(PvHoek.D_HXK,PvHoek.HXK,0.1,'pchip','extrap');

    v_GLT_predicted(experiment)=p.GLT.VmGLT.*(f.GLCo-GLCi./p.GLT.KeqGLT)./...
        (p.GLT.KmGLTGLCo.*(1+f.GLCo/p.GLT.KmGLTGLCo+GLCi / p.GLT.KmGLTGLCi + ...
        0.91.*f.GLCo.*GLCi./(p.GLT.KmGLTGLCi.*p.GLT.KmGLTGLCo)));

    v_GLK_predicted(experiment)=p.HXK_ExprsCor.*(((p.HXK1_kcat.*...
        (f.HXK1+f.HXK2))./(p.HXK1_Katp.*p.HXK1_Kglc).*(ATP.*...
        GLCi-((ADP.*G6P)./p.HXK1_Keq)))./((1+ATP./p.HXK1_Katp+ADP./...
        p.HXK1_Kadp).*(1+GLCi./p.HXK1_Kglc+G6P./p.HXK1_Kg6p+T6P./p.HXK1_Kt6p)));
end
 
%plotting options
if setup.parEst.drawnowCF == 1
    
    figure(1000)
    
    subplot(2,2,1),
    plot(v_GLT_predicted,v_GLT_data,'r*')
    line([0 max(v_GLT_data)*1.2],[0 max(v_GLT_data)*1.2])
    xlabel('predicted flux (mM/s)')
    ylabel('measured flux (mM/s)')
    title('GLT')
    
    subplot(2,2,3),
    plot(v_GLK_predicted,v_HXK_data,'r*')
    line([0 max(v_HXK_data)*1.2],[0 max(v_HXK_data)*1.2])
    xlabel('predicted flux (mM/s)')
    ylabel('measured flux (mM/s)')
    title('HXK')

%     if setup.parEst.costfun2 == 1
%         % glc_i        
% %         subplot(3,2,3)
%         
%         % g6p
%         subplot(3,2,4)
%         plot(G6P_predicted,G6P,'r*')
%         line([0 max(G6P)*1.2],[0 max(G6P)*1.2])
%         xlabel('predicted concentration [mM]')
%         ylabel('measured concentration [mM]')
%         title('G6P')    
%     end
    subplot(2,2,2) 
    % parameters
    bar(x_temp)
    title('parameter estimates')
    
    delete(findall(gcf,'type','annotation'))

%     dim = [.3 .5 0 0.1];
%     str = ['  p_{1}=  ', num2str(x_temp(1)), ',  p_{2}=  ', num2str(x_temp(2)), ',  p_{3}=  ', num2str(x_temp(3)), ',  p_{4}=  ', num2str(x_temp(4)), ',  p_{5}=  ', num2str(x_temp(5)), ',  p_{6}=  ', num2str(x_temp(6)), ',  p_{7}=  ', num2str(x_temp(7)), ',  p_{8}=  ', num2str(x_temp(8)), ',  p_{9}=  ', num2str(x_temp(9))];
    
    dim = [0.65 .35 0 .1];
    
%     str = ['  p_{1}=  ', num2str(x_temp(1)), ',  p_{2}=  ', num2str(x_temp(2)), ',  p_{3}=  ', num2str(x_temp(3)), ',  p_{4}=  ', num2str(x_temp(4))];
    str1 = ['  p_{1}=  ', num2str(x_temp(1)),'  p_{2}=  ', num2str(x_temp(2))];
    str2 = ['  p_{3}=  ', num2str(x_temp(3)),'  p_{4}=  ', num2str(x_temp(4))];
    str3 = ['  p_{5}=  ', num2str(x_temp(5)),'  p_{6}=  ', num2str(x_temp(6))];
    str4 = ['  p_{7}=  ', num2str(x_temp(7)),'  p_{8}=  ', num2str(x_temp(8))];
    str5 = ['  p_{9}=  ', num2str(x_temp(9))];
    str = {str1, str2, str3, str4, str5};    annotation('textbox',dim,'string',str,'FitBoxToText','on')
    
%     drawnow()
    
    
    figure(1001)
    
    subplot(2,2,1),
    plot(darray, v_GLT_data, 'r*')
    hold on
    plot(darray, v_GLT_predicted, 'k-')
    hold off
    xlabel('dilution rate [h^{-1}]')
    ylabel('GLT flux [mM s^{-1}]')
    legend('experimental data','simulated','location','northwest')
     
    subplot(2,2,3),
    plot(darray, v_GLK_data, 'r*')
    hold on
    plot(darray, v_GLK_predicted, 'k-')
    hold off
    xlabel('dilution rate [h^{-1}]')
    ylabel('HXK flux [mM s^{-1}]')
    legend('experimental data','simulated','location','northwest')

%     if setup.parEst.costfun2 == 1
%         % glc_i        
%         subplot(3,2,3)
%         plot(darray, GLCi_predicted, 'k-')
%         xlabel('dilution rate [h^{-1}]')
%         ylabel('G6P concentration [mM]')
%         
%         % g6p
%         subplot(3,2,4)
%         plot(darray, G6P, 'r*')
%         hold on
%         plot(darray, G6P_predicted, 'k-')
%         hold off
%         xlabel('dilution rate [h^{-1}]')
%         ylabel('G6P concentration [mM]')
%               
%     end
    subplot(2,2,2) 
    % parameters
    bar(x_temp)
    title('parameter estimates')
    
    delete(findall(gcf,'type','annotation'))
%     dim = [.3 .5 0 0.1];
%     str = ['  p_{1}=  ', num2str(x_temp(1)), ',  p_{2}=  ', num2str(x_temp(2)), ',  p_{3}=  ', num2str(x_temp(3)), ',  p_{4}=  ', num2str(x_temp(4)), ',  p_{5}=  ', num2str(x_temp(5)), ',  p_{6}=  ', num2str(x_temp(6)), ',  p_{7}=  ', num2str(x_temp(7)), ',  p_{8}=  ', num2str(x_temp(8)), ',  p_{9}=  ', num2str(x_temp(9))];
    
    dim = [0.65 .35 0 .1];
    
%     str = ['  p_{1}=  ', num2str(x_temp(1)), ',  p_{2}=  ', num2str(x_temp(2)), ',  p_{3}=  ', num2str(x_temp(3)), ',  p_{4}=  ', num2str(x_temp(4))];
    str1 = ['  p_{1}=  ', num2str(x_temp(1)),'  p_{2}=  ', num2str(x_temp(2))];
    str2 = ['  p_{3}=  ', num2str(x_temp(3)),'  p_{4}=  ', num2str(x_temp(4))];
    str3 = ['  p_{5}=  ', num2str(x_temp(5)),'  p_{6}=  ', num2str(x_temp(6))];
    str4 = ['  p_{7}=  ', num2str(x_temp(7)),'  p_{8}=  ', num2str(x_temp(8))];
    str5 = ['  p_{9}=  ', num2str(x_temp(9))];
    str = {str1, str2, str3, str4, str5};

    annotation('textbox',dim,'string',str,'FitBoxToText','on')
%     
%     drawnow()
end

% different types of cost function
lambda = setup.parEst.lambda;
errorGLT=v_GLT_predicted-v_GLT_data;
errorHXK=v_GLK_predicted-v_HXK_data;
error=[errorGLT,errorHXK,lambda.*x_temp(1:end)];

end
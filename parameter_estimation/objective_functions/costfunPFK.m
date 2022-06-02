function [error]=costfunPFK(x_temp,canelas_SS,setup,x)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

% adjust parameters
x(setup.caseStudy.parameters) = x_temp;

% set parameter structure
setParameterStructure_Y3M1;
if(isfield(setup,'adj_lit')&&(setup.adj_lit == 1))
%     p.PFK_kcat=.2*10.^x(56).*3.55; %(UNIT!)
    p.PFK_kcat=10.^x(56).*3.55; %(UNIT!)
end

% select data
F6P=canelas_SS.mtlD.F6P;
ATP=canelas_SS.mtlD.ATP;
F16BP=canelas_SS.mtlD.FBP;
AMP=canelas_SS.mtlD.AMP;
v_PFK_data=canelas_SS.mtlD.v_PFK;
lambda = setup.parEst.lambda;
darray = canelas_SS.mtlD.D;
if length(lambda) == 1
    d=canelas_SS.mtlD.D;
else
    d=canelas_SS.mtlD.D(d);
end
PvHoek_EnzymeExpressionData

if setup.parEst.costfun == 1
    % simulate SS reaction rates
    F26BP=canelas_SS.mtlD.F26BP;
    PFK_nom=(p.PFK_kcat.*f.PFK.*p.PFK_gR.*(F6P./p.PFK_Kf6p).*(ATP./p.PFK_Katp).*(1+(F6P./p.PFK_Kf6p)+(ATP./p.PFK_Katp)+p.PFK_gR.*((F6P./p.PFK_Kf6p).*(ATP./p.PFK_Katp))));
    PFK_denom=(1+F6P./p.PFK_Kf6p+ATP./p.PFK_Katp+(p.PFK_gR.*(F6P./p.PFK_Kf6p).*(ATP./p.PFK_Katp))).^2+...
        p.PFK_L.*...
        ((1+p.PFK_Ciatp.*(ATP./p.PFK_Kiatp))./(1+ATP./p.PFK_Kiatp)).^2.*...
        ((1+p.PFK_Camp.*(AMP./p.PFK_Kamp))./(1+AMP./p.PFK_Kamp)).^2.*...
        ((1+((p.PFK_Cf26bp*F26BP)./(p.PFK_Kf26bp))+((p.PFK_Cf16bp.*F16BP)./(p.PFK_Kf16bp)))./(1+(F26BP./p.PFK_Kf26bp)+(F16BP./p.PFK_Kf16bp))).^2.*...
        (1+p.PFK_Catp.*(ATP./p.PFK_Katp)).^2;
    v_PFK_predicted=p.PFK_ExprsCor.*(PFK_nom./PFK_denom);
end

% plotting options
if setup.parEst.drawnowCF == 1

    figure(1000)
    
    subplot(2,2,1)
    plot(v_PFK_predicted,v_PFK_data,'r*')
    line([0 max(v_PFK_data)*1.2],[0 max(v_PFK_data)*1.2])
    xlabel('predicted flux (mM/s)')
    ylabel('measured flux (mM/s)')
    title('PFK')
    
    subplot(2,2,3)
    plot(darray, v_PFK_data, 'r*')
    hold on
    plot(darray, v_PFK_predicted, 'k-')
    hold off
    xlabel('dilution rate [h^{-1}]')
    ylabel('PFK flux [mM s^{-1}]')
    legend('experimental data','simulated','location','northwest')
    
    subplot(2,2,2)
    bar(x_temp)
    title('parameter estimates')
    
    delete(findall(gcf,'type','annotation'))
%     dim = [.3 .5 0 0.1];
    dim = [0.55 .35 0 .1];
%     str = ['  p_{1}=  ', num2str(x_temp(1)), ',  p_{2}=  ', num2str(x_temp(2)), ',  p_{3}=  ', num2str(x_temp(3)), ',  p_{4}=  ', num2str(x_temp(4)), ',  p_{5}=  ', num2str(x_temp(5)), ',  p_{6}=  ', num2str(x_temp(6)), ',  p_{7}=  ', num2str(x_temp(7)), ',  p_{8}=  ', num2str(x_temp(8)), ',  p_{9}=  ', num2str(x_temp(9))];
    str1 = ['  p_{1}=  ', num2str(x_temp(1)),'  p_{2}=  ', num2str(x_temp(2))];
    str2 = ['  p_{3}=  ', num2str(x_temp(3)),'  p_{4}=  ', num2str(x_temp(4))];
    str3 = ['  p_{5}=  ', num2str(x_temp(5)),'  p_{6}=  ', num2str(x_temp(6))];
    str4 = ['  p_{7}=  ', num2str(x_temp(7)),'  p_{8}=  ', num2str(x_temp(8))];
    str5 = ['  p_{9}=  ', num2str(x_temp(9)),'  p_{10}=  ', num2str(x_temp(10))];
    str6 = ['  p_{11}=  ', num2str(x_temp(11)),'  p_{12}=  ', num2str(x_temp(12))];
    str7 = ['  p_{13}=  ', num2str(x_temp(13)),'  p_{14}=  ', num2str(x_temp(14))];
%     str = {str1, str2, str3, str4, str5, str6, str7, str8, str9, str10, str11, str12, str13, str14}; 
    str = {str1, str2, str3, str4, str5, str6, str7};    
    annotation('textbox',dim,'string',str,'FitBoxToText','on')
    
% % % %     
%     
%     figure(1000)
%     
%     subplot(2,2,1)
%     plot(v_PFK_predicted,v_PFK_data,'r*')
%     line([0 max(v_PFK_data)*1.2],[0 max(v_PFK_data)*1.2])
%     xlabel('predicted flux (mM/s)')
%     ylabel('measured flux (mM/s)')
%     title('PFK') 
%     
%     if setup.parEst.costfun2 == 1
%         subplot(2,2,3)
%         plot(F6P_predicted,F6P,'r*')
%         line([0 max(F6P)*1.2],[0 max(F6P)*1.2])
%         xlabel('predicted concentration [mM]')
%         ylabel('measured concentration [mM]')
%         title('F6P') 
% 
%         subplot(2,2,4)
%         plot(F16BP_predicted,F16BP,'r*')
%         line([0 max(F16BP)*1.2],[0 max(F16BP)*1.2])
%         xlabel('predicted concentration [mM]')
%         ylabel('measured concentration [mM]')
%         title('F16BP')
%     end
%     subplot(2,2,2)
%     bar(x_temp)
%     title('parameter estimates')
%     
%     delete(findall(gcf,'type','annotation'))
%     dim = [.3 .5 0 0.1];
%     str = ['  p_{1}=  ', num2str(x_temp(1)), ',  p_{2}=  ', num2str(x_temp(2)), ',  p_{3}=  ', num2str(x_temp(3)), ',  p_{4}=  ', num2str(x_temp(4)), ',  p_{5}=  ', num2str(x_temp(5)), ',  p_{6}=  ', num2str(x_temp(6)), ',  p_{7}=  ', num2str(x_temp(7)), ',  p_{8}=  ', num2str(x_temp(8)), ',  p_{9}=  ', num2str(x_temp(9))];
%     annotation('textbox',dim,'string',str,'FitBoxToText','on')
%     
%     drawnow()
%     
%     figure(1001)
%     
%     subplot(2,2,1)
%     plot(darray, v_PFK_data, 'r*')
%     hold on
%     plot(darray, v_PFK_predicted, 'k-')
%     hold off
%     xlabel('dilution rate [h^{-1}]')
%     ylabel('PFK flux [mM s^{-1}]')
%     legend('experimental data','simulated','location','northwest')
%     
%     if setup.parEst.costfun2 == 1
%         subplot(2,2,3)
%         plot(darray, F6P, 'r*')
%         hold on
%         plot(darray, F6P_predicted, 'k-')
%         hold off
%         xlabel('dilution rate [h^{-1}]')
%         ylabel('F6P concentration [mM]')
% 
%         subplot(2,2,4)
%         plot(darray, F16BP, 'r*')
%         hold on
%         plot(darray, F16BP_predicted, 'k-')
%         hold off
%         xlabel('dilution rate [h^{-1}]')
%         ylabel('F16BP concentration [mM]')
%     end
%     subplot(2,2,2)
%     bar(x_temp)
%     title('parameter estimates')
%     
%     delete(findall(gcf,'type','annotation'))
%     dim = [.3 .5 0 0.1];
%     str = ['  p_{1}=  ', num2str(x_temp(1)), ',  p_{2}=  ', num2str(x_temp(2)), ',  p_{3}=  ', num2str(x_temp(3)), ',  p_{4}=  ', num2str(x_temp(4)), ',  p_{5}=  ', num2str(x_temp(5)), ',  p_{6}=  ', num2str(x_temp(6)), ',  p_{7}=  ', num2str(x_temp(7)), ',  p_{8}=  ', num2str(x_temp(8)), ',  p_{9}=  ', num2str(x_temp(9))];
%     annotation('textbox',dim,'string',str,'FitBoxToText','on')
%     
%     drawnow() 
% % % % 
end

% different types of cost function
if setup.parEst.costfun == 1
    error = [(v_PFK_data-v_PFK_predicted),lambda.*x_temp([1:14])];
end

end
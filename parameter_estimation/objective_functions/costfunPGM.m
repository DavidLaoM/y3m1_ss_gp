function [error]=costfunPGM(x_temp,canelas_SS,setup,x)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

% adjust parameters
x(setup.caseStudy.parameters) = x_temp;

% set parameter structure
setParameterStructure_Y3M1;

% select data
P3G=canelas_SS.mtlD.P3G;
P2G=canelas_SS.mtlD.P2G;
v_PGM_data=canelas_SS.mtlD.v_PGM;
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
    v_PGM_predicted=p.PGM_ExprsCor.*((((p.GPM1_kcat.*(f.GPM1+f.GPM2+f.GPM3))./p.GPM1_K3pg).*(P3G-P2G./p.GPM1_Keq))./...
        (1+P3G./p.GPM1_K3pg+1+P2G./p.GPM1_K2pg-1));
end

% plotting options
if setup.parEst.drawnowCF == 1
    figure(1000)
    
    subplot(2,2,1)
    plot(v_PGM_predicted,v_PGM_data,'r*')
    line([0 max(v_PGM_data)*1.2],[0 max(v_PGM_data)*1.2])
    xlabel('predicted flux (mM/s)')
    ylabel('measured flux (mM/s)')
    title('PGM') 
    
    subplot(2,2,3)
    plot(darray, v_PGM_data, 'r*')
    hold on
    plot(darray, v_PGM_predicted, 'k-')
    hold off
    xlabel('dilution rate [h^{-1}]')
    ylabel('PGM flux [mM s^{-1}]')
    legend('experimental data','simulated','location','northwest')

    subplot(2,2,2)
    bar(x_temp)
    title('parameter estimates')
    
    delete(findall(gcf,'type','annotation'))
%     dim = [.3 .5 0 0.1];
    dim = [0.65 .35 0 .1];
    
%     str = ['  p_{1}=  ', num2str(x_temp(1)), ',  p_{2}=  ', num2str(x_temp(2)), ',  p_{3}=  ', num2str(x_temp(3)), ',  p_{4}=  ', num2str(x_temp(4))];
    str1 = ['  p_{1}=  ', num2str(x_temp(1))];
    str2 = ['  p_{2}=  ', num2str(x_temp(2))];
    str3 = ['  p_{3}=  ', num2str(x_temp(3))];
    str4 = ['  p_{4}=  ', num2str(x_temp(4))];
    str = {str1, str2, str3, str4};
    
    annotation('textbox',dim,'string',str,'FitBoxToText','on')
    
    drawnow()
end
% if setup.parEst.drawnow == 1
%     figure(1000)
%     
%     subplot(2,2,1)
%     plot(v_PGM_predicted,v_PGM_data,'r*')
%     line([0 max(v_PGM_data)*1.2],[0 max(v_PGM_data)*1.2])
%     xlabel('predicted flux (mM/s)')
%     ylabel('measured flux (mM/s)')
%     title('PGM') 
%     
%     if setup.parEst.costfun2 == 1
%     subplot(2,2,3)
%     plot(P2G_predicted,P2G,'r*')
%     line([0 max(P2G)*1.2],[0 max(P2G)*1.2])
%     xlabel('predicted concentration [mM]')
%     ylabel('measured concentration [mM]')
%     title('P2G') 
%     
%     subplot(2,2,4)
%     plot(P3G_predicted,P3G,'r*')
%     line([0 max(P3G)*1.2],[0 max(P3G)*1.2])
%     xlabel('predicted concentration [mM]')
%     ylabel('measured concentration [mM]')
%     title('P3G')
%     end
%     subplot(2,2,2)
%     bar(x_temp)
%     title('parameter estimates')
%     
%     delete(findall(gcf,'type','annotation'))
%     dim = [.3 .5 0 0.1];
%     str = ['  p_{1}=  ', num2str(x_temp(1)), ',  p_{2}=  ', num2str(x_temp(2)), ',  p_{3}=  ', num2str(x_temp(3)), ',  p_{4}=  ', num2str(x_temp(4))];
%     annotation('textbox',dim,'string',str,'FitBoxToText','on')
%     
%     drawnow()
%     
%     figure(1001)
%     
%     subplot(2,2,1)
%     plot(darray, v_PGM_data, 'r*')
%     hold on
%     plot(darray, v_PGM_predicted, 'k-')
%     hold off
%     xlabel('dilution rate [h^{-1}]')
%     ylabel('PGM flux [mM s^{-1}]')
%     legend('experimental data','simulated','location','northwest')
%     
%     if setup.parEst.costfun2 == 1
%         subplot(2,2,3)
%         plot(darray, P3G, 'r*')
%         hold on
%         plot(darray, P3G_predicted, 'k-')
%         hold off
%         xlabel('dilution rate [h^{-1}]')
%         ylabel('P3G concentration [mM]')
% %         legend('experimental data','simulated','location','northwest')
% 
%         subplot(2,2,4)
%         plot(darray, P2G, 'r*')
%         hold on
%         plot(darray, P2G_predicted, 'k-')
%         hold off
%         xlabel('dilution rate [h^{-1}]')
%         ylabel('P2G concentration [mM]')
% %         legend('experimental data','simulated','location','northwest')
%     end
%     subplot(2,2,2)
%     bar(x_temp)
%     title('parameter estimates')
%     
%     delete(findall(gcf,'type','annotation'))
%     dim = [.3 .5 0 0.1];
%     str = ['  p_{1}=  ', num2str(x_temp(1)), ',  p_{2}=  ', num2str(x_temp(2)), ',  p_{3}=  ', num2str(x_temp(3)), ',  p_{4}=  ', num2str(x_temp(4))];
%     annotation('textbox',dim,'string',str,'FitBoxToText','on')
%     
%     drawnow()   
% % % % %     
% end

% different types of cost function
if setup.parEst.costfun == 1
    error = [(v_PGM_data-v_PGM_predicted),lambda.*x_temp([1:4])];
end

end
function [error]=costfunENO(x_temp,canelas_SS,setup,x)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

% adjust parameters
x(setup.caseStudy.parameters) = x_temp;

% set parameter structure
setParameterStructure_Y3M1;
if(isfield(setup,'adj_lit')&&(setup.adj_lit == 1))
%     p.ENO1_Keq=4.011*10.^x(17); %[]
    p.ENO1_Keq=6.7*10.^x(17); %[]
end

% select data
PEP=canelas_SS.mtlD.PEP;
P2G=canelas_SS.mtlD.P2G;
v_ENO_data=canelas_SS.mtlD.v_ENO;
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
    v_ENO_predicted=p.ENO_ExprsCor.*((((p.ENO1_kcat.*(f.ENO1))./p.ENO1_K2pg).*(P2G-PEP./p.ENO1_Keq))./...
        (1+P2G./p.ENO1_K2pg+1+PEP./p.ENO1_Kpep-1));
end

% plotting options
if setup.parEst.drawnowCF == 1
    figure(1000)
    
    subplot(2,2,1)
    plot(v_ENO_predicted,v_ENO_data,'r*')
    line([0 max(v_ENO_data)*1.2],[0 max(v_ENO_data)*1.2])
    xlabel('predicted flux (mM/s)')
    ylabel('measured flux (mM/s)')
    title('ENO') 
    
    subplot(2,2,3)
    plot(darray, v_ENO_data, 'r*')
    hold on
    plot(darray, v_ENO_predicted, 'k-')
    hold off
    xlabel('dilution rate [h^{-1}]')
    ylabel('ENO flux [mM s^{-1}]')
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
    
% %     
% %     figure(1001)
% %     
% %     subplot(2,2,1)
% %     plot(darray, v_ENO_data, 'r*')
% %     hold on
% %     plot(darray, v_ENO_predicted, 'k-')
% %     hold off
% %     xlabel('dilution rate [h^{-1}]')
% %     ylabel('ENO flux [mM s^{-1}]')
% %     legend('experimental data','simulated','location','northwest')
% %     
% % % %
% %     subplot(2,2,2)
% %     bar(x_temp)
% %     title('parameter estimates')
% %     
% %     delete(findall(gcf,'type','annotation'))
% %     dim = [.3 .5 0 0.1];
% %     str = ['  p_{1}=  ', num2str(x_temp(1)), ',  p_{2}=  ', num2str(x_temp(2)), ',  p_{3}=  ', num2str(x_temp(3)), ',  p_{4}=  ', num2str(x_temp(4))];
% %     annotation('textbox',dim,'string',str,'FitBoxToText','on')
% %     
% %     drawnow()   
    
%     figure(1000)
%     
%     subplot(2,2,1)
%     plot(v_ENO_predicted,v_ENO_data,'r*')
%     line([0 max(v_ENO_data)*1.2],[0 max(v_ENO_data)*1.2])
%     xlabel('predicted flux (mM/s)')
%     ylabel('measured flux (mM/s)')
%     title('ENO') 
%     
%     if(setup.parEst.costfun2 == 1 || setup.parEst.costfun3 == 1)
%         subplot(2,2,3)
%         plot(P2G_predicted,P2G,'r*')
%         line([0 max(P2G)*1.2],[0 max(P2G)*1.2])
%         xlabel('predicted concentration [mM]')
%         ylabel('measured concentration [mM]')
%         title('P2G') 
% 
%         subplot(2,2,4)
%         plot(PEP_predicted,PEP,'r*')
%         line([0 max(PEP)*1.2],[0 max(PEP)*1.2])
%         xlabel('predicted concentration [mM]')
%         ylabel('measured concentration [mM]')
%         title('PEP')
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
%     plot(darray, v_ENO_data, 'r*')
%     hold on
%     plot(darray, v_ENO_predicted, 'k-')
%     hold off
%     xlabel('dilution rate [h^{-1}]')
%     ylabel('ENO flux [mM s^{-1}]')
%     legend('experimental data','simulated','location','northwest')
%     
%     if(setup.parEst.costfun2 == 1 || setup.parEst.costfun3 == 1)
%         subplot(2,2,3)
%         plot(darray, P2G, 'r*')
%         hold on
%         plot(darray, P2G_predicted, 'k-')
%         hold off
%         xlabel('dilution rate [h^{-1}]')
%         ylabel('P2G concentration [mM]')
% %         legend('experimental data','simulated','location','northwest')
% 
%         subplot(2,2,4)
%         plot(darray, PEP, 'r*')
%         hold on
%         plot(darray, PEP_predicted, 'k-')
%         hold off
%         xlabel('dilution rate [h^{-1}]')
%         ylabel('PEP concentration [mM]')
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
% %     persistent piter
% %     persistent pparam
% %     if sum(x_temp) == 0
% %         piter   = 1;
% %         pparam  = x_temp;
% %     else
% %         piter = [piter; piter + 1];
% %         pparam = [pparam; x_temp];
% %     end
% %     
% %     subplot(2,2,2)
% %     plot(piter(end), pparam(end,:),'-o')
% %     hold on
%     
%     
%     drawnow()    
    
end

% different types of cost function
if setup.parEst.costfun == 1
    error = [(v_ENO_data-v_ENO_predicted),lambda.*x_temp([1:end])];
end

end
function [error] = costfunPGI(x_temp,canelas_SS,setup,x)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

% adjust parameters
x(setup.indexParameters.PGI) = x_temp;

% set parameter structure
setParameterStructure_Y3M1;
if(isfield(setup,'adj_lit')&&(setup.adj_lit == 1))
%     p.PGI1_Keq=10.^x(57).*0.2586; % []
    p.PGI1_Keq=10.^x(57).*0.314; % []
end

% select data
F6P=canelas_SS.mtlD.F6P;
G6P=canelas_SS.mtlD.G6P;
v_PGI_data=canelas_SS.mtlD.v_PGI;
lambda = setup.parEst.lambda;
darray = canelas_SS.mtlD.D;
if length(lambda) == 1
    d=canelas_SS.mtlD.D;
else
    d=canelas_SS.mtlD.D(d);
end
PvHoek_EnzymeExpressionData;

if setup.parEst.costfun == 1
    v_PGI_predicted=p.PGI_ExprsCor.*((((p.PGI1_kcat*f.PGI1)./p.PGI1_Kg6p).*(G6P-(F6P./p.PGI1_Keq)))./...
        (1+G6P./p.PGI1_Kg6p+1+F6P./p.PGI1_Kf6p-1));
end

% plotting options
if setup.parEst.drawnowCF == 1
    
    figure(1000)
    
    subplot(2,2,1)
    plot(v_PGI_predicted,v_PGI_data,'r*')
    line([0 max(v_PGI_data)*1.2],[0 max(v_PGI_data)*1.2])
    xlabel('predicted flux (mM/s)')
    ylabel('measured flux (mM/s)')
    title('PGI') 
     
    subplot(2,2,3)
    plot(darray, v_PGI_data, 'r*')
    hold on
    plot(darray, v_PGI_predicted, 'k-')
    hold off
    xlabel('dilution rate [h^{-1}]')
    ylabel('PGI flux [mM s^{-1}]')
    legend('experimental data','simulated','location','northwest')
    
    subplot(2,2,2)
    bar(x_temp)
    title('parameter estimates')
    
    delete(findall(gcf,'type','annotation'))
%     dim = [.3 .5 0.1 0.1];
%     dim = [.3 .5 0 0.1];
%     str = ['  p_{1}=  ', num2str(x_temp(1)), ',  p_{2}=  ', num2str(x_temp(2)), ',  p_{3}=  ', num2str(x_temp(3)), ',  p_{4}=  ', num2str(x_temp(4))];
    dim = [0.65 .35 0 .1];
    
%     str = ['  p_{1}=  ', num2str(x_temp(1)), ',  p_{2}=  ', num2str(x_temp(2)), ',  p_{3}=  ', num2str(x_temp(3)), ',  p_{4}=  ', num2str(x_temp(4))];
    str1 = ['  p_{1}=  ', num2str(x_temp(1))];
    str2 = ['  p_{2}=  ', num2str(x_temp(2))];
    str3 = ['  p_{3}=  ', num2str(x_temp(3))];
    str4 = ['  p_{4}=  ', num2str(x_temp(4))];
    str = {str1, str2, str3, str4};

    annotation('textbox',dim,'string',str,'FitBoxToText','on')
    
    drawnow()


    
 
    
    
    
% %     figure(1000)
% %     
% %     subplot(2,2,1)
% %     plot(v_PGI_predicted,v_PGI_data,'r*')
% %     line([0 max(v_PGI_data)*1.2],[0 max(v_PGI_data)*1.2])
% %     xlabel('predicted flux (mM/s)')
% %     ylabel('measured flux (mM/s)')
% %     title('PGI') 
% %     
% %     if setup.parEst.costfun2 == 1
% %         subplot(2,2,3)
% %         plot(G6P_predicted,G6P,'r*')
% %         line([0 max(G6P)*1.2],[0 max(G6P)*1.2])
% %         xlabel('predicted concentration [mM]')
% %         ylabel('measured concentration [mM]')
% %         title('G6P') 
% % 
% %         subplot(2,2,4)
% %         plot(F6P_predicted,F6P,'r*')
% %         line([0 max(F6P)*1.2],[0 max(F6P)*1.2])
% %         xlabel('predicted concentration [mM]')
% %         ylabel('measured concentration [mM]')
% %         title('F6P')
% %         %     drawnow
% %         
% %         % insert textbox with parameter values
% %         
% %     end
% %     subplot(2,2,2)
% %     bar(x_temp)
% %     title('parameter estimates')
% %     
% %     delete(findall(gcf,'type','annotation'))
% % %     dim = [.3 .5 0.1 0.1];
% %     dim = [.3 .5 0 0.1];
% %     str = ['  p_{1}=  ', num2str(x_temp(1)), ',  p_{2}=  ', num2str(x_temp(2)), ',  p_{3}=  ', num2str(x_temp(3)), ',  p_{4}=  ', num2str(x_temp(4))];
% %     annotation('textbox',dim,'string',str,'FitBoxToText','on')
% %     
% %     drawnow()
% %     
% %     figure(1001)
% %     
% %     subplot(2,2,1)
% %     plot(darray, v_PGI_data, 'r*')
% %     hold on
% %     plot(darray, v_PGI_predicted, 'k-')
% %     hold off
% %     xlabel('dilution rate [h^{-1}]')
% %     ylabel('PGI flux [mM s^{-1}]')
% %     legend('experimental data','simulated','location','northwest')
% %     
% %     if setup.parEst.costfun2 == 1
% %         subplot(2,2,3)
% %         plot(darray, G6P, 'r*')
% %         hold on
% %         plot(darray, G6P_predicted, 'k-')
% %         hold off
% %         xlabel('dilution rate [h^{-1}]')
% %         ylabel('G6P concentration [mM]')
% % %         legend('experimental data','simulated','location','northwest')
% % 
% %         subplot(2,2,4)
% %         plot(darray, F6P, 'r*')
% %         hold on
% %         plot(darray, F6P_predicted, 'k-')
% %         hold off
% %         xlabel('dilution rate [h^{-1}]')
% %         ylabel('F6P concentration [mM]')
% % %         legend('experimental data','simulated','location','northwest')
% %     end
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
% % % % % % %     
end

% different types of cost function
if setup.parEst.costfun == 1
    error = [(v_PGI_data-v_PGI_predicted),lambda.*x_temp([1:end])];
end


end
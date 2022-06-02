function [error]=costfunTPI(x_temp,canelas_SS,setup,x)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

% adjust parameters
x(setup.caseStudy.parameters) = x_temp;

% set parameter structure
setParameterStructure_Y3M1;
if(isfield(setup,'adj_lit')&&(setup.adj_lit == 1))
%     p.TPI1_Keq=10.^x(80).*0.0391; %[]
%     p.TPI1_kcat=10.^x(82).*78.396; %mM/s
    p.TPI1_Keq=10.^x(80).*0.045; %[]
    p.TPI1_kcat=10.^x(82).*37.74; %mM/s
end

% select data
GLYCERAL3P = canelas_SS.mtlD.GAP;
DHAP = canelas_SS.mtlD.DHAP;
v_TPI_data = canelas_SS.mtlD.v_TPI;
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
    v_TPI_predicted=p.TPI_ExprsCor.*(((p.TPI1_kcat*f.TPI1)./p.TPI1_Kdhap.*(DHAP-GLYCERAL3P./p.TPI1_Keq))./...
        (1+DHAP./p.TPI1_Kdhap+1+GLYCERAL3P./p.TPI1_Kglyceral3p-1));
end

% plotting options
if setup.parEst.drawnowCF == 1
    figure(1000)
    
    subplot(2,2,1)
    plot(v_TPI_predicted,v_TPI_data,'r*')
    line([0 max(v_TPI_data)*1.2],[0 max(v_TPI_data)*1.2])
    xlabel('predicted flux (mM/s)')
    ylabel('measured flux (mM/s)')
    title('TPI') 
    
    subplot(2,2,3)
    plot(darray, v_TPI_data, 'r*')
    hold on
    plot(darray, v_TPI_predicted, 'k-')
    hold off
    xlabel('dilution rate [h^{-1}]')
    ylabel('TPI flux [mM s^{-1}]')
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
%     plot(v_TPI_predicted,v_TPI_data,'r*')
%     line([0 max(v_TPI_data)*1.2],[0 max(v_TPI_data)*1.2])
%     xlabel('predicted flux (mM/s)')
%     ylabel('measured flux (mM/s)')
%     title('TPI')
%     
% %     if setup.parEst.costfun2 == 1
% %     subplot(2,2,3)
% %     plot(P2G_predicted,P2G,'r*')
% %     line([0 max(P2G)*1.2],[0 max(P2G)*1.2])
% %     xlabel('predicted concentration [mM]')
% %     ylabel('measured concentration [mM]')
% %     title('P2G') 
% %     
% %     subplot(2,2,4)
% %     plot(PEP_predicted,PEP,'r*')
% %     line([0 max(PEP)*1.2],[0 max(PEP)*1.2])
% %     xlabel('predicted concentration [mM]')
% %     ylabel('measured concentration [mM]')
% %     title('PEP')
% %     end
% end

% different types of cost function
if setup.parEst.costfun == 1
    error = [(v_TPI_data-v_TPI_predicted),lambda.*x_temp([1:4])];
end

end
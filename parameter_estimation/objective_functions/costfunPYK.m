function [error]=costfunPYK(x_temp,canelas_SS,setup,x)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
plength = length(setup.caseStudy.parameters);

% adjust parameters
x(setup.caseStudy.parameters) = x_temp;

% set parameter structure
setParameterStructure_Y3M1;

% select data
ADP=canelas_SS.mtlD.ADP;
ATP=canelas_SS.mtlD.ATP;
F16BP=canelas_SS.mtlD.FBP;
PEP=canelas_SS.mtlD.PEP;
PYR=canelas_SS.mtlD.PYR;
v_PYK_data=canelas_SS.mtlD.v_PYK;
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
    v_PYK_predicted=p.PYK_ExprsCor.*((((p.PYK1_kcat.*(f.PYK1+f.PYK2))./(p.PYK1_Kadp.*p.PYK1_Kpep).*ADP.*PEP)./...
        ((1+ADP./p.PYK1_Kadp).*(1+PEP./p.PYK1_Kpep))).*...
        ((PEP./p.PYK1_Kpep+1).^p.PYK1_hill./(p.PYK1_L.*((ATP./p.PYK1_Katp+1)./(F16BP./p.PYK1_Kf16bp+1)).^p.PYK1_hill+(PEP./p.PYK1_Kpep+1).^p.PYK1_hill)));
end

% plotting options
if setup.parEst.drawnowCF == 1
    figure(1000)
    
    subplot(2,2,1)
    plot(v_PYK_predicted,v_PYK_data,'r*')
    line([0 max(v_PYK_data)*1.2],[0 max(v_PYK_data)*1.2])
    xlabel('predicted flux (mM/s)')
    ylabel('measured flux (mM/s)')
    title('PYK') 
    
    subplot(2,2,3)
    plot(darray, v_PYK_data, 'r*')
    hold on
    plot(darray, v_PYK_predicted, 'k-')
    hold off
    xlabel('dilution rate [h^{-1}]')
    ylabel('PYK flux [mM s^{-1}]')
    legend('experimental data','simulated','location','northwest')

    subplot(2,2,2)
    bar(x_temp)
    title('parameter estimates')
    
    delete(findall(gcf,'type','annotation'))
%     dim = [.3 .5 0 0.1];
    dim = [0.65 .35 0 .1];
    
%     str = ['  p_{1}=  ', num2str(x_temp(1)), ',  p_{2}=  ', num2str(x_temp(2)), ',  p_{3}=  ', num2str(x_temp(3)), ',  p_{4}=  ', num2str(x_temp(4))];
    str1 = ['  p_{1}=  ', num2str(x_temp(1)),'  p_{2}=  ', num2str(x_temp(2))];
    str2 = ['  p_{3}=  ', num2str(x_temp(3)),'  p_{4}=  ', num2str(x_temp(4))];
    str3 = ['  p_{5}=  ', num2str(x_temp(5)),'  p_{6}=  ', num2str(x_temp(6))];
    str4 = ['  p_{7}=  ', num2str(x_temp(7))];
    str = {str1, str2, str3, str4};
    
    annotation('textbox',dim,'string',str,'FitBoxToText','on')
    
    drawnow()
end
% if setup.parEst.drawnow == 1
%     figure(1000)
%     
%     subplot(2,2,1)
%     plot(v_PYK_predicted,v_PYK_data,'r*')
%     line([0 max(v_PYK_data)*1.2],[0 max(v_PYK_data)*1.2])
%     xlabel('predicted flux (mM/s)')
%     ylabel('measured flux (mM/s)')
%     title('PYK') 
%     
%     if setup.parEst.costfun2 == 1
%     subplot(2,2,3)
%     plot(PYR_predicted,PYR,'r*')
%     line([0 max(PYR)*1.2],[0 max(PYR)*1.2])
%     xlabel('predicted concentration [mM]')
%     ylabel('measured concentration [mM]')
%     title('PYR') 
%     
%     subplot(2,2,4)
%     plot(PEP_predicted,PEP,'r*')
%     line([0 max(PEP)*1.2],[0 max(PEP)*1.2])
%     xlabel('predicted concentration [mM]')
%     ylabel('measured concentration [mM]')
%     title('PEP')
%     end
%     disp('visualization to be updated here.')
% end

% different types of cost function
if setup.parEst.costfun == 1
    error = [(v_PYK_data-v_PYK_predicted),lambda.*x_temp([1:plength])];
end

end
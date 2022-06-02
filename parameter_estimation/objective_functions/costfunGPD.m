function [error]=costfunGPD(x_temp,canelas_SS,setup,x)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

% adjust parameters
x(setup.caseStudy.parameters) = x_temp;

% set parameter structure
setParameterStructure_Y3M1;

% select data
ADP=canelas_SS.mtlD.ADP;
ATP=canelas_SS.mtlD.ATP;
F16BP=canelas_SS.mtlD.FBP;
DHAP=canelas_SS.mtlD.DHAP;
NADH=canelas_SS.mtlD.NADH;
NAD=canelas_SS.mtlD.NAD;
GLYC3P=canelas_SS.mtlD.G3P;
v_G3PDH_data=canelas_SS.mtlD.v_G3PDH;
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
    v_G3PDH_predicted=((((p.GPD1_kcat.*f.GPD1)./(p.GPD1_Kdhap.*p.GPD1_Knadh)).*...
        (DHAP.*NADH-(GLYC3P.*NAD)./p.GPD1_Keq))./...
        ((1+F16BP./p.GPD1_Kf16bp+ATP./p.GPD1_Katp+ADP./p.GPD1_Kadp).*...
        (1+DHAP./p.GPD1_Kdhap+GLYC3P./p.GPD1_Kglyc3p).*(1+NADH./p.GPD1_Knadh+NAD./p.GPD1_Knad)));
end

% plotting options
if setup.parEst.drawnowCF == 1
    figure(1000)
    
    subplot(2,2,1)
    plot(v_G3PDH_predicted,v_G3PDH_data,'r*')
    line([0 max(v_G3PDH_data)*1.2],[0 max(v_G3PDH_data)*1.2])
    xlabel('predicted flux (mM/s)')
    ylabel('measured flux (mM/s)')
    title('G3PDH') 
    
    subplot(2,2,3)
    plot(darray, v_G3PDH_data, 'r*')
    hold on
    plot(darray, v_G3PDH_predicted, 'k-')
    hold off
    xlabel('dilution rate [h^{-1}]')
    ylabel('G3PDH flux [mM s^{-1}]')
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
    str4 = ['  p_{7}=  ', num2str(x_temp(7)),'  p_{8}=  ', num2str(x_temp(8))];
    str5 = ['  p_{9}=  ', num2str(x_temp(9))];
   str = {str1, str2, str3, str4, str5};
    
    annotation('textbox',dim,'string',str,'FitBoxToText','on')
    
    drawnow()
end
% if setup.parEst.drawnow == 1
%     figure(1000)
%     
%     subplot(2,2,1)
%     plot(v_G3PDH_predicted,v_G3PDH_data,'r*')
%     line([0 max(v_G3PDH_data)*1.2],[0 max(v_G3PDH_data)*1.2])
%     xlabel('predicted flux (mM/s)')
%     ylabel('measured flux (mM/s)')
%     title('GPD') 
%     
%     if setup.parEst.costfun2 == 1
% %         subplot(2,2,3)
% %         plot(P2G_predicted,P2G,'r*')
% %         line([0 max(P2G)*1.2],[0 max(P2G)*1.2])
% %         xlabel('predicted concentration [mM]')
% %         ylabel('measured concentration [mM]')
% %         title('P2G') 
% % 
% %         subplot(2,2,4)
% %         plot(PEP_predicted,PEP,'r*')
% %         line([0 max(PEP)*1.2],[0 max(PEP)*1.2])
% %         xlabel('predicted concentration [mM]')
% %         ylabel('measured concentration [mM]')
% %         title('PEP')
%     end
% end

% different types of cost function
if setup.parEst.costfun == 1
    error = [(v_G3PDH_data-v_G3PDH_predicted),lambda.*x_temp([1:9])];
end

end
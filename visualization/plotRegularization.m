function [h] = plotRegularization(output,setup)
lambdalist          = output.regularization.ltot;
errorParams_array   = output.regularization.errorParams_array;
errorData_array     = output.regularization.errorData_array;

h = figure;
% h = 1;
% figure(100);
[AX,H1,H2] = plotyy(lambdalist,errorParams_array,lambdalist,errorData_array,@loglog);
%setting colors
left_color = [0, 0, 1];
right_color = [1, 0, 0];
set(H1,'Marker','o','Color',left_color)
set(H2,'Marker','o','Color',right_color)
set(AX(1),'yColor',left_color)
set(AX(2),'yColor',right_color)
set(AX,'defaultAxesColorOrder',[left_color; right_color]);
legend('SSE parameters','SSE data','Location','northeast')
%pointing the lambda value chosen
% % % % ylim([min(plRes)-0.1 max(plRes)+0.1])
% % % % y1 = ylim;
% % % % line([par(i) par(i)],y1,'Color','magenta','LineStyle',':') % this is the right parameter valu
lam = setup.parEst.lambda1;
line([lam lam],[AX(1).YLim(1) AX(1).YLim(2)],'Color','black') 
% display the value of lam. textbox. put next to black line on top left or
% below the plot.
textlam = ['lambda = ', num2str(lam)];
ylim = get(gca,'ylim');
xlim = get(gca,'xlim');
% text(xlim(1)+0.025,ylim(2)-0.1,textlam,'fontsize',12);
text(xlim(1)+0.000001,ylim(2),textlam,'fontsize',12);
title(['regularization for enzyme ', setup.legenda.fluxes{setup.caseStudy.fluxes}])

end
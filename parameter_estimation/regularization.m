function [h,output] = regularization(canelas_SS,setup)
%UNTITLED13 Summary of this function goes here
%   Detailed explanation goes here
refpars     = setup.caseStudy.parameters;
npars       = size(refpars);
refvals     = setup.refVals.parameters(refpars);
lambdalist  = setup.parEst.lambdalist;

% select setup options
load('parameters_blank.mat');
options = optimset('Display','iter');
plength = length(setup.caseStudy.parameters);
x_temp(plength)=0;
lb = -3*ones(1,plength);
ub =  3*ones(1,plength);
ltot = length(setup.parEst.lambdalist);

xres_array          = zeros(ltot,plength);
% resnorm_array       = zeros(ltot,1);
exitflag_array      = zeros(ltot,1);
errorData_array     = zeros(ltot,1);
errorLambda_array   = zeros(ltot,1);
errorParams         = zeros(npars(2),ltot);
errorParams_array   = zeros(ltot,1);

for i = 1:ltot
    setup.parEst.lambda = setup.parEst.lambdalist(i);
    disp(['estimation for lambda =', num2str(setup.parEst.lambda)]);
    if setup.caseStudy.PGI == 1
        [xres,resnorm,residual,exitflag] = lsqnonlin(@costfunPGI,x_temp,lb,ub,options,canelas_SS,setup,x);
    elseif setup.caseStudy.PGM == 1
        [xres,resnorm,residual,exitflag] = lsqnonlin(@costfunPGM,x_temp,lb,ub,options,canelas_SS,setup,x);
    elseif setup.caseStudy.ENO == 1
        [xres,resnorm,residual,exitflag] = lsqnonlin(@costfunENO,x_temp,lb,ub,options,canelas_SS,setup,x);
    elseif setup.caseStudy.PYK == 1
        [xres,resnorm,residual,exitflag] = lsqnonlin(@costfunPYK,x_temp,lb,ub,options,canelas_SS,setup,x);
    elseif setup.caseStudy.PDC == 1
        [xres,resnorm,residual,exitflag] = lsqnonlin(@costfunPDC,x_temp,lb,ub,options,canelas_SS,setup,x);
    elseif setup.caseStudy.PFK == 1
        [xres,resnorm,residual,exitflag] = lsqnonlin(@costfunPFK,x_temp,lb,ub,options,canelas_SS,setup,x);
    elseif setup.caseStudy.FBA == 1
        [xres,resnorm,residual,exitflag] = lsqnonlin(@costfunFBA,x_temp,lb,ub,options,canelas_SS,setup,x);
    elseif setup.caseStudy.TPI == 1
        [xres,resnorm,residual,exitflag] = lsqnonlin(@costfunTPI,x_temp,lb,ub,options,canelas_SS,setup,x);
    elseif setup.caseStudy.GPD == 1
        [xres,resnorm,residual,exitflag] = lsqnonlin(@costfunGPD,x_temp,lb,ub,options,canelas_SS,setup,x);
    elseif setup.caseStudy.GLT_HXK == 1
        [xres,resnorm,residual,exitflag] = lsqnonlin(@costfunGLT_HXK,x_temp,lb,ub,options,canelas_SS,setup,x);
    elseif setup.caseStudy.GAPDH_PGK == 1
        [xres,resnorm,residual,exitflag] = lsqnonlin(@costfunGAPDH_PGK,x_temp,lb,ub,options,canelas_SS,setup,x);
    elseif setup.caseStudy.PDC_ADH == 1
        [xres,resnorm,residual,exitflag] = lsqnonlin(@costfunPDC_ADH,x_temp,lb,ub,options,canelas_SS,setup,x);
    elseif setup.caseStudy.PDC_ADH_branches == 1
        [xres,resnorm,residual,exitflag] = lsqnonlin(@costfunPDC_ADH_BRANCHES,x_temp,lb,ub,options,canelas_SS,setup,x);
    elseif setup.caseStudy.GLT_HXK_old == 1
        [xres,resnorm,residual,exitflag] = lsqnonlin(@costfunGLT_HXK_old,x_temp,lb,ub,options,canelas_SS,setup,x);
    %elseif
    end
    
    %error for parameters
    for j = 1:npars(2)
        errorParams(j,i)     = refvals(j) .* 10 .^ (xres(j));
    end
    
    % put that inside the arrays again
    xres_array(i,:)          = xres;
%     resnorm_array(i,:)       = resnorm;
    exitflag_array(i,:)      = exitflag;
    errorData_array(i,:)     = sum(residual(1:end-4).^2);
    errorLambda_array(i,:)   = sum(residual(end-3:end).^2);
    errorParams_array(i,:)   = sum(errorParams(:,i)); 
%     disp(errorParams_array(i,:));
end

% % plots1
% h =figure;
% left_color = [0, 0, 1];
% right_color = [1, 0, 0];
% set(h,'defaultAxesColorOrder',[left_color; right_color]);
% yyaxis right
% loglog(lambdalist,errorData_array,'ro-')
% % ylabel('SSE parameters')
% hold on
% % loglog(setup.parEst.lambdalist,errorLambda_array,'bo-')
% yyaxis left
% loglog(lambdalist,errorParams_array,'bo-');
% % legend('SSE data','SSE regularisation','Location','southeast')
% legend('SSE data','SSE parameters','Location','northeast')
% xlabel('lambda')
% % ylabel('SSE data')
% title(['regularization for enzyme ', setup.legenda.fluxes{setup.caseStudy.fluxes}])

if setup.parEst.drawnow == 1
    %plots2
    figure(100)
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

    h = AX;
else
    h = 0;
end

output.ltot                 = setup.parEst.lambdalist;
output.xres_array           = xres_array;
output.exitflag_array       = exitflag_array;
output.errorData_array      = errorData_array;
output.errorLambda_array    = errorLambda_array;
output.errorParams_array    = errorParams_array;
end

% %% Alternative plotting start
% % % to chage
% % lambdalist          = outputPFKss.regularization.ltot;
% % errorParams_array   = outputPFKss.regularization.errorLambda_array;
% % errorData_array     = outputPFKss.regularization.errorData_array;
% lambdalist          = outputENOss.regularization.ltot;
% errorParams_array   = outputENOss.regularization.errorLambda_array;
% errorData_array     = outputENOss.regularization.errorData_array;
% % % to run
% figure(100)
% [AX,H1,H2] = plotyy(lambdalist,errorParams_array,lambdalist,errorData_array,@loglog);
% %setting colors
% left_color = [0, 0, 1];
% right_color = [1, 0, 0];
% set(H1,'Marker','o','Color',left_color)
% set(H2,'Marker','o','Color',right_color)
% set(AX(1),'yColor',left_color)
% set(AX(2),'yColor',right_color)
% set(AX,'defaultAxesColorOrder',[left_color; right_color]);
% legend('SSE parameters','SSE data','Location','northeast')
% %pointing the lambda value chosen
% % % % % ylim([min(plRes)-0.1 max(plRes)+0.1])
% % % % % y1 = ylim;
% % % % % line([par(i) par(i)],y1,'Color','magenta','LineStyle',':') % this is the right parameter valu
% lam = setup.parEst.lambda1;
% line([lam lam],[AX(1).YLim(1) AX(1).YLim(2)],'Color','black') 
% % display the value of lam. textbox. put next to black line on top left or
% % below the plot.
% textlam = ['lambda = ', num2str(lam)];
% ylim = get(gca,'ylim');
% xlim = get(gca,'xlim');
% % text(xlim(1)+0.025,ylim(2)-0.1,textlam,'fontsize',12);
% text(xlim(1)+0.000001,ylim(2),textlam,'fontsize',12);
% title(['regularization for enzyme ', setup.legenda.fluxes{setup.caseStudy.fluxes}])
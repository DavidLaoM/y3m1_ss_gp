function [h,pProfile] = runPLA(x_temp,x,canelas_SS,setup)
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here

% bypass cost function
% % % % setup.parEst.costfun1   = 0;
% % % % setup.parEst.costfun2   = 1;
% setup.parEst.costfun1   = 1;
% setup.parEst.costfun2   = 0;

if setup.parEst.nonreg == 1 && setup.parEst.reg == 0
    setup.parEst.lambda = setup.parEst.lambda0;
elseif setup.parEst.reg == 1 && setup.parEst.nonreg == 0
    setup.parEst.lambda = setup.parEst.lambda1;
end


ptot = length(setup.caseStudy.parameters);
names = {'p1', 'p2', 'p3', 'p4', 'p5', 'p6', 'p7', 'p8', 'p9', 'p10', 'p11', 'p12', 'p13', 'p14', 'p15', 'p16', 'p17', 'p18', 'p19', 'p20'};
for i = 1:ptot
% % % % for i = 1:2
% for i = 4
    par         = x_temp; 
    index       = setup.caseStudy.parameters(i);
    parname     = setup.legenda.parameters{index};%name
    setup.i     = i;

    threshold   = chi2inv(0.5,size(par)); % threshold - chi square distribution, threshold = chi2inv(0.5,size(par));
    lb          = -3 * ones(1,ptot); %lower bounds for parameter estimation (lsqnonlin)
    ub          = +3 * ones(1,ptot); %upper bounds for parameter estimation (lsqnonlin)  
    
    if setup.caseStudy.PGI == 1
        func        = @(par)costfunPGI(par,canelas_SS,setup,x);
    elseif setup.caseStudy.PGM == 1
        func        = @(par)costfunPGM(par,canelas_SS,setup,x);
    elseif setup.caseStudy.ENO == 1
        func        = @(par)costfunENO(par,canelas_SS,setup,x);
    elseif setup.caseStudy.PYK == 1
        func        = @(par)costfunPYK(par,canelas_SS,setup,x);
    elseif setup.caseStudy.PDC == 1
        func        = @(par)costfunPDC(par,canelas_SS,setup,x);
    elseif setup.caseStudy.PFK == 1
        func        = @(par)costfunPFK(par,canelas_SS,setup,x);
    elseif setup.caseStudy.FBA == 1
        func        = @(par)costfunFBA(par,canelas_SS,setup,x);
    elseif setup.caseStudy.TPI == 1
        func        = @(par)costfunTPI(par,canelas_SS,setup,x);
    elseif setup.caseStudy.GPD == 1
        func        = @(par)costfunGPD(par,canelas_SS,setup,x);
    elseif setup.caseStudy.GLT_HXK == 1
        func        = @(par)costfunGLT_HXK(par,canelas_SS,setup,x);
    elseif setup.caseStudy.GAPDH_PGK == 1
        func        = @(par)costfunGAPDH_PGK(par,canelas_SS,setup,x);
    elseif setup.caseStudy.PDC_ADH == 1
        func        = @(par)costfunPDC_ADH(par,canelas_SS,setup,x);
    elseif setup.caseStudy.PDC_ADH_branches == 1
        func        = @(par)costfunPDC_ADH_BRANCHES(par,canelas_SS,setup,x);
    elseif setup.caseStudy.HOR2_GlycT == 1
        func        = @(par)costfunHOR2_GlycT(par,canelas_SS,setup,x);
    % %     elseif
    end

%     threshold   = chi2inv(0.5,size(par)); % threshold - chi square distribution, threshold = chi2inv(0.5,size(par));
%     lb          = -3 * ones(1,ptot); %lower bounds for parameter estimation (lsqnonlin)
%     ub          = +3 * ones(1,ptot); %upper bounds for parameter estimation (lsqnonlin)
    
    Optimoptions= optimoptions('lsqnonlin','Display','off'); % options for optimization
    minStep     = 0.01; % minimal step factor
    maxStep     = 0.1; % maximal step factor
    minChange   = 0.001; % minimal change of resnorm
    maxChange   = 0.05; % maximal change of resnorm
    nr          = 100; % no. of samples in profile likelihood

    maxPar  = par(i) * 1.1;

%     if i == 7
%         maxPar_up   = 1;
%         maxPar_down = -0.9;
%     else
        if setup.parEst.nonreg == 1 && setup.parEst.reg == 0
%             if i == 3
%                 maxPar_up   = 2.4;
%                 maxPar_down = -3;
            if (setup.caseStudy.PYK == 1 && i==6)
                maxPar_up   = 1;
                maxPar_down = -1;
            elseif (setup.caseStudy.PYK == 1 && i==7)
                maxPar_up   = 0.8;
                maxPar_down = -0.8;                
            elseif (setup.caseStudy.PDC == 1 && i==3)
                maxPar_up   = 2;
                maxPar_down = -3;
            else
                maxPar_up   = 3;
                maxPar_down = -3;
            end
        elseif setup.parEst.reg == 1 && setup.parEst.nonreg == 0
            if(setup.caseStudy.PYK == 1 && i==7)
                maxPar_up   = 0.8;
                maxPar_down = -0.8;
            elseif(setup.caseStudy.PYK == 1 && i==6)
                maxPar_up   = 1.8;
                maxPar_down = -0.8;
            elseif(setup.caseStudy.PFK == 1 && (i == 1 || i == 3 || i == 9))
                maxPar_up   = 3;
                maxPar_down = -3;
            else
%                 maxPar_up   = 1;
%                 maxPar_down = -1;
                maxPar_up   = 3;
                maxPar_down = -3;
                
            end
        end

%     end

%     [plPar,plRes]=PLA(func,par,i,maxPar,threshold,lb,ub,Optimoptions,minStep,maxStep,minChange,maxChange,nr);
%     [plPar,plRes]=PLA2(func,par,i,maxPar,threshold,lb,ub,Optimoptions,minStep,maxStep,minChange,maxChange,nr,setup);
%     [plPar,plRes]=PLA3(func,par,i,maxPar,threshold,lb,ub,Optimoptions,minStep,maxStep,minChange,maxChange,nr);
    [plPar,plRes]=PLA4(func,par,i,maxPar,threshold,lb,ub,Optimoptions,minStep,maxStep,minChange,maxChange,nr,maxPar_up,maxPar_down,setup);
    pProfile.pPar.(names{i})     = plPar';
    pProfile.pRes.(names{i})     = plRes';
    
    if setup.parEst.drawnow == 1
        h=figure(99);
        if ((setup.caseStudy.PYK == 1) || (setup.caseStudy.FBA == 1) || (setup.caseStudy.GPD == 1) || (setup.caseStudy.GLT_HXK == 1))
            subplot(3,3,i)
        elseif ((setup.caseStudy.PFK == 1) || (setup.caseStudy.GAPDH_PGK == 1) || (setup.caseStudy.PDC_ADH == 1))
            subplot(4,4,i)
        else
            subplot(2,2,i)
        end
        plot(plPar, plRes)
        title([parname])
        xlabel('Parameter value')
        ylabel('Fit residual')
% %         if setup.parEst.nonreg == 1 && setup.parEst.reg == 0
            xlim([-3 3])
% %         elseif setup.parEst.reg == 1 && setup.parEst.nonreg == 0
% %             if(setup.caseStudy.PFK == 1 && (i == 1 || i == 3 || i == 9))
% %                 xlim([-3 3])
% %             else
% %                 xlim([-1 1])
% %             end
% %         end
        hold on
        [val,ind] = min(plRes(~ismember(plRes,0)));
        ylim([min(plRes)-0.1 max(plRes)+0.1])
        y1 = ylim;
        line([par(i) par(i)],y1,'Color','magenta','LineStyle',':') % this is the right parameter value
        hold on
        line([plPar(ind) plPar(ind)],y1,'Color','red','LineStyle','--') % this is the current minimum
        hold on
        drawnow()
    else
        h = 0;
    end
end

if setup.parEst.drawnow == 1
    if setup.parEst.nonreg == 1
        suptitle(['NO reg. Profile Likelihood Analysis (PLA)', setup.legenda.fluxes(setup.caseStudy.fluxes)])
    elseif setup.parEst.reg == 1
        suptitle(['reg. Profile Likelihood Analysis (PLA)', setup.legenda.fluxes(setup.caseStudy.fluxes)])
    end
    hold off
end

end

% %% Alternative plotting start
% % to chage
% % Par = struct2cell(outputENOss.PLAnonReg.pPar);
% % Res = struct2cell(outputENOss.PLAnonReg.pRes);
% % par = outputENOss.MSparestNonReg.b;
% Par = struct2cell(outputPFKss.PLAnonReg.pPar);
% Res = struct2cell(outputPFKss.PLAnonReg.pRes);
% par = outputPFKss.MSparestNonReg.b;
% h=figure(99);
% % for i = 1:4
% for i = 1:14
% % to run
%     plPar = Par{i};
%     plRes = Res{i};
%     index       = setup.caseStudy.parameters(i);
%     parname     = setup.legenda.parameters{index};%name
%     
%     if ((setup.caseStudy.PYK == 1) || (setup.caseStudy.FBA == 1) || (setup.caseStudy.GPD == 1) || (setup.caseStudy.GLT_HXK == 1))
%         subplot(3,3,i)
%     elseif ((setup.caseStudy.PFK == 1) || (setup.caseStudy.GAPDH_PGK == 1) || (setup.caseStudy.PDC_ADH == 1))
%         subplot(4,4,i)
%     else
%         subplot(2,2,i)
%     end
%     plot(plPar, plRes)
%     title([parname])
%     xlabel('Parameter value')
%     ylabel('Fit residual')
%     if setup.parEst.nonreg == 1 && setup.parEst.reg == 0
%         xlim([-3 3])
%     elseif setup.parEst.reg == 1 && setup.parEst.nonreg == 0
%         if(setup.caseStudy.PFK == 1 && (i == 1 || i == 3 || i == 9))
%             xlim([-3 3])
%         else
%             xlim([-1 1])
%         end
%     end
%     hold on
%     [val,ind] = min(plRes(~ismember(plRes,0)));
%     ylim([min(plRes)-0.1 max(plRes)+0.1])
%     y1 = ylim;
%     line([par(i) par(i)],y1,'Color','magenta','LineStyle',':') % this is the right parameter value
%     hold on
%     line([plPar(ind) plPar(ind)],y1,'Color','red','LineStyle','--') % this is the current minimum
%     hold on
%     drawnow()
% end
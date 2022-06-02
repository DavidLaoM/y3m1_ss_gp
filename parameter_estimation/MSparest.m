function [b,fval,exitflag,output,solutions] = MSparest(canelas_SS,setup,x)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

load('parameters_blank.mat');
plength = length(setup.caseStudy.parameters);
x_temp(plength)=0;
lb = -3*ones(1,plength); %lb = -1*ones(1,plength); %lb(3) = -0.5; lb(4) = -1.5;
ub =  3*ones(1,plength); %ub =  1*ones(1,plength); %ub(3) = -0.5; ub(4) = 1.5;
% lb = -1*ones(1,plength); %lb(3) = -0.5; lb(4) = -1.5;
% ub =  1*ones(1,plength); %ub(3) = -0.5; ub(4) = 1.5;

% set the lambda to regularized or not regularized
if setup.parEst.MSnonReg == 1
    setup.parEst.lambda = setup.parEst.lambda0;
elseif setup.parEst.MSReg == 1
	setup.parEst.lambda = setup.parEst.lambda1;
end
% select enzyme for study
if setup.caseStudy.PGI == 1
    model = @(xres)costfunPGI(xres,canelas_SS,setup,x);
elseif setup.caseStudy.PGM == 1
    model = @(xres)costfunPGM(xres,canelas_SS,setup,x);   
elseif setup.caseStudy.ENO == 1
    model = @(xres)costfunENO(xres,canelas_SS,setup,x);
elseif setup.caseStudy.PYK == 1
    model = @(xres)costfunPYK(xres,canelas_SS,setup,x);
elseif setup.caseStudy.PDC == 1
    model = @(xres)costfunPDC(xres,canelas_SS,setup,x);
elseif setup.caseStudy.PFK == 1
    model = @(xres)costfunPFK(xres,canelas_SS,setup,x);
elseif setup.caseStudy.FBA == 1
    model = @(xres)costfunFBA(xres,canelas_SS,setup,x);
elseif setup.caseStudy.TPI == 1
    model = @(xres)costfunTPI(xres,canelas_SS,setup,x);
elseif setup.caseStudy.GPD == 1
    model = @(xres)costfunGPD(xres,canelas_SS,setup,x);
elseif setup.caseStudy.GLT_HXK == 1
    model = @(xres)costfunGLT_HXK(xres,canelas_SS,setup,x);
elseif setup.caseStudy.GAPDH_PGK == 1
    model = @(xres)costfunGAPDH_PGK(xres,canelas_SS,setup,x);
elseif setup.caseStudy.HOR2_GlycT == 1
    model = @(xres)costfunHOR2_GlycT(xres,canelas_SS,setup,x);
elseif setup.caseStudy.PDC_ADH == 1
    model = @(xres)costfunPDC_ADH(xres,canelas_SS,setup,x);
elseif setup.caseStudy.PDC_ADH_branches == 1
    model = @(xres)costfunPDC_ADH_BRANCHES(xres,canelas_SS,setup,x);
%elseif
end

%prepare and run multistart itself 
problem = createOptimProblem('lsqnonlin', 'objective', model, ...
    'xdata', canelas_SS, 'ydata', canelas_SS, 'x0', x_temp, 'lb', lb, ...
    'ub', ub);
ms = MultiStart('Display','iter');
% ms = MultiStart('Display','iter', 'PlotFcns', {@gsplotbestf,@gsplotfunccount});

if setup.parEst.MSnonReg == 1
    if setup.caseStudy.PFK == 1
        [b,fval,exitflag,output,solutions] = run(ms, problem, 500);
    else
        [b,fval,exitflag,output,solutions] = run(ms, problem, 50);
    end
elseif setup.parEst.MSReg == 1
	[b,fval,exitflag,output,solutions] = run(ms, problem, 50);
end

end


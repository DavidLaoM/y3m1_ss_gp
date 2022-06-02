function [xres,resnorm,residual,exitflag] = runLSQsimple(canelas_SS,setup)
%UNTITLED11 Summary of this function goes here
%   Detailed explanation goes here

% select setup options
load('parameters_blank.mat');
options = optimset('Display','iter');
% options = optimset('Display','iter','OutputFcn',@saveParameterIterates);
% if setup.caseStudy.ENO == 1
%     options = optimset('Display','iter','PlotFcns',@plotEstimates);
% end
plength = length(setup.caseStudy.parameters);
x_temp(plength)=0;
% x_temp=[0 0 0 0];
lb = -3*ones(1,plength);
ub =  3*ones(1,plength);
setup.parEst.lambda = setup.parEst.lambda0;
if((setup.parEst.MSReg == 1)&&(setup.parEst.MSnonReg == 0))
    setup.parEst.lambda = setup.parEst.lambda1;
end

% run lsqnonlinwith objective function
if setup.caseStudy.PGM == 1
    [xres,resnorm,residual,exitflag] = lsqnonlin(@costfunPGM,x_temp,lb,ub,options,canelas_SS,setup,x);
elseif setup.caseStudy.ENO == 1
%     [xres,resnorm,residual,exitflag] = lsqnonlin(@costfunENO,x_temp,lb,ub,options,canelas_SS,setup,x);
    % testing the output function here
    global history
    history.x = [];%.x; x];
    history.iteration = [];%.iteration; optimValues.iteration];
    history.funccount = [];%.funccount; optimValues.funccount];
    history.stepsize = [];%.stepsize; optimValues.stepsize];
    history.gradient = [];%.gradient; optimValues.gradient'];
    history.firstorderopt = [];%.firstorderopt; optimValues.firstorderopt];
    history.cgiterations = [];%.cgiteration; optimValues.cgiteration];
    history.positivedefinite = [];%.positivedefinite; optimValues.positivedefinite];
    history.ratio = [];%.ratio; optimValues.ratio];
    history.degenerate = [];%.degenerate; optimValues.degenerate];
    history.trustregionradius = [];%.trustregionradius; optimValues.trustregionradius];
    history.residual = [];%.residual; optimValues.residual'];
    history.resnorm = [];%.resnorm; optimValues.resnorm];
%     history.x = [];
%     history.fval = [];
%     searchdir = [];
%     options = optimset('Display','iter');
    options = optimset('Display','iter','OutputFcn',{@saveIterationsENO});
    [xres,resnorm,residual,exitflag] = lsqnonlin(@costfunENO,x_temp,lb,ub,options,canelas_SS,setup,x);
    exitflag2 = exitflag;
    clear exitflag
    exitflag.exitflag = exitflag2;
    exitflag.history = history;
elseif setup.caseStudy.PYK == 1
    [xres,resnorm,residual,exitflag] = lsqnonlin(@costfunPYK,x_temp,lb,ub,options,canelas_SS,setup,x);
elseif setup.caseStudy.PDC == 1
    [xres,resnorm,residual,exitflag] = lsqnonlin(@costfunPDC,x_temp,lb,ub,options,canelas_SS,setup,x);
elseif setup.caseStudy.PGI == 1
    [xres,resnorm,residual,exitflag] = lsqnonlin(@costfunPGI,x_temp,lb,ub,options,canelas_SS,setup,x);
elseif setup.caseStudy.FBA == 1
    [xres,resnorm,residual,exitflag] = lsqnonlin(@costfunFBA,x_temp,lb,ub,options,canelas_SS,setup,x);
elseif setup.caseStudy.TPI == 1
    [xres,resnorm,residual,exitflag] = lsqnonlin(@costfunTPI,x_temp,lb,ub,options,canelas_SS,setup,x);
elseif setup.caseStudy.PFK == 1
    [xres,resnorm,residual,exitflag] = lsqnonlin(@costfunPFK,x_temp,lb,ub,options,canelas_SS,setup,x);
elseif setup.caseStudy.GPD == 1
    [xres,resnorm,residual,exitflag] = lsqnonlin(@costfunGPD,x_temp,lb,ub,options,canelas_SS,setup,x);
elseif setup.caseStudy.GLT_HXK == 1
    [xres,resnorm,residual,exitflag] = lsqnonlin(@costfunGLT_HXK,x_temp,lb,ub,options,canelas_SS,setup,x);
elseif setup.caseStudy.GAPDH_PGK == 1
    [xres,resnorm,residual,exitflag] = lsqnonlin(@costfunGAPDH_PGK,x_temp,lb,ub,options,canelas_SS,setup,x);
elseif setup.caseStudy.HOR2_GlycT == 1
    [xres,resnorm,residual,exitflag] = lsqnonlin(@costfunHOR2_GlycT,x_temp,lb,ub,options,canelas_SS,setup,x);
elseif setup.caseStudy.PDC_ADH == 1
    [xres,resnorm,residual,exitflag] = lsqnonlin(@costfunPDC_ADH,x_temp,lb,ub,options,canelas_SS,setup,x);
elseif setup.caseStudy.PDC_ADH_branches == 1
    [xres,resnorm,residual,exitflag] = lsqnonlin(@costfunPDC_ADH_BRANCHES,x_temp,lb,ub,options,canelas_SS,setup,x);
    % elseif
end
end


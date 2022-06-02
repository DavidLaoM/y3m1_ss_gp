% SETUP FOR STUDY

%% case study data
setup.caseStudy.PGI         = 1;
setup.caseStudy.PGM         = 0;
setup.caseStudy.PYK         = 0;
setup.caseStudy.ENO         = 0;
setup.caseStudy.PDC         = 0;
setup.caseStudy.PFK         = 0;
setup.caseStudy.FBA         = 0;
setup.caseStudy.GPD         = 0;
setup.caseStudy.GLT_HXK     = 0;
setup.caseStudy.GAPDH_PGK   = 0;
setup.caseStudy.PDC_ADH     = 0;
selectCaseStudy;

%% single PSA
setup.sPSA.numEvals         = 101;

%% subsystem analysis
setup.subsystem.numEvals    = 150; %400; %250;

setup.subsystem.all         = 0;
% setup.subsystem.PGI         = 1;

%% parameterEstimation
setup.parEst.drawnow        = 0; %0;
setup.parEst.drawnowCF      = 0; %0;
setup.parEst.lambda0        = 0;
setup.parEst.lambdalist     = [1E-5, 2E-5, 5E-5,1E-4, 2E-4, 5E-4, ...
                                1E-3, 2E-3, 5E-3, 1E-2, 2E-2, 5E-2, ...
                                1E-1, 2E-1, 5E-1, 1E0, 2E0, 5E0, ...
                                1E1, 2E1, 5E1, 1E2, 2E2, 5E2,...
                                1E3, 2E3, 5E3];
setup.parEst.lambda1        = 0.5;
                            
setup.parEst.costfun       = 1; % cost function with only SS fluxes

setup.parEst.PGI            = 1;

setup.parEst.MSnonReg       = 1;
setup.parEst.MSReg          = 0;
setup.parEst.nonreg         = 1; 
setup.parEst.reg            = 0;

%% parameters index
setup.indexParameters.PGI  = 57:60;

%% 2021/01/16 laters addition
setup.startParamsOFF = 1;
setup.treCycleON = 0;
setup.importPiON = 0;
setup.conditionsSS = 0;
setup.conditionsGP = 0;
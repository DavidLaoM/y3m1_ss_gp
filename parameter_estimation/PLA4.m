function [plPar,plRes]=PLA3(func,par,i,maxPar,threshold,lb,ub,Optimoptions,minStep,maxStep,minChange,maxChange,nr,maxPar_up,maxPar_down,setup)

disp(' ');disp(['Profile Likelihood calculation for parameter ' int2str(i)]);disp(' ');

% set default values
if nargin < 4, maxPar = 10*par(i);end
if nargin < 5, threshold = chi2inv(0.5,size(par));end   
if nargin < 6, lb=[]; end
if nargin < 7, ub=[]; end
if nargin < 8, Optimoptions=[]; end
if nargin < 9, minStep=0.01; end
if nargin < 10, maxStep=0.1; end
if nargin < 11, minChange=0.001; end
if nargin < 12, maxChange=0.05; end
if nargin < 13, nr=100; end

minChange = minChange * threshold;
maxChange = maxChange * threshold;
% %----Choice 1----------%
% minStep = par(i)*minStep;
% maxStep = par(i)*maxStep;
% %----Choice 2---------%
minStep = abs(par(i)-maxPar)*minStep;
maxStep = abs(par(i)-maxPar)*maxStep;

step = minStep;

% specify the direction: 1 increase -1 decrease
flag = sign(maxPar - par(i));

% select parameters to be re-optimized
nPar = [par(1:i-1) par(i+1:end)];
if isempty(lb), lb=[]; else lb = [lb(1:i-1) lb(i+1:end)]; end
if isempty(ub), ub=[]; else ub = [ub(1:i-1) ub(i+1:end)]; end

% function handle for optimization
optFun = @(nPar,iPar)feval(func,[nPar(1:i-1) iPar nPar(i:end)]); 
res    = sum(func(par).^2);

% % % INCREASE PROCESS % % % 
k_inc = 1;
count_inc = 1; % Lijkt dubbel op met k
plPar_inc(1) = par(i);
plRes_inc(1) = res; % was res
step_inc = step;

[x_start, resnorm_start] = lsqnonlin(@(nPar)optFun(nPar, par(i)), nPar(k_inc,:), lb, ub, Optimoptions);

% % % % % Late additions to simplify teh code (David 2019-09-17)
% maxPar_up   = 1;
% maxPar_down = -1;
% maxPar_up       = 3;
% maxPar_down     = -3;
if ((setup.caseStudy.GLT_HXK==1) || (setup.caseStudy.PDC_ADH == 1))%||())
    step = 0.1;
elseif(setup.caseStudy.GAPDH_PGK == 1)
    step = 0.5;
else
    step = 0.01;
end
step_inc        = step;
step_dec        = step;

while ((plPar_inc(end)<=maxPar_up && logical(flag+1)) || (plPar_inc(end)>=maxPar_down && ~logical(flag+1)))

    tempPar = flag * step_inc + plPar_inc(end);
% % % %     disp(tempPar);
    [x,resnorm] = lsqnonlin( @(nPar)optFun(nPar,tempPar),nPar(k_inc,:),lb,ub,Optimoptions);
    
    resChange = resnorm - plRes_inc(end);
    count_inc = count_inc+1;
    k_inc = k_inc+1;
    plPar_inc(k_inc) = tempPar; 
    plRes_inc(k_inc) = resnorm;
    nPar(k_inc,:)= x;
end

% % % DECREASE PROCESS % % % 
k_dec = 1;
count_dec = 1; % Lijkt dubbel op met k
plPar_dec(1) = par(i);
plRes_dec(1) = res; % was res
step_dec = step;

nPar = [par(1:i-1) par(i+1:end)];

while ((plPar_dec(end)>=maxPar_down && logical(flag+1)) || (plPar_dec(end)<=maxPar_up && ~logical(flag+1)))
    tempPar = - flag * step_dec + plPar_dec(end);
% % % %     disp(tempPar);
%     if setup.i == 7
%         disp(tempPar)
%     end
    
    [x,resnorm] = lsqnonlin( @(nPar)optFun(nPar,tempPar),nPar(k_dec,:),lb,ub,Optimoptions);

    resChange = resnorm - plRes_dec(end);
    count_dec = count_dec+1;
    k_dec = k_dec+1;
    plPar_dec(k_dec) = tempPar; 
    plRes_dec(k_dec) = resnorm;
    nPar(k_dec,:)= x;
end


plPar = [flip(plPar_dec) par(i) plPar_inc];
plRes = [flip(plRes_dec) resnorm_start plRes_inc];

end
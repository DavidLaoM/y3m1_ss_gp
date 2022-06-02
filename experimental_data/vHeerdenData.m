% this file loads the 'vHeerden_trehalose_data_micromolgdw.mat' file,
% changes the units and inferrs some variables.
% Instead of micromolgdw, now the units are mM.

% load experimental data
load('vHeerden_trehalose_data_micromolgdw') %load vdHeerden GP dataset

%convert to mmol/L assuming 1 gdw = 2 mL
data.nucleotides=data.nucleotides/2;
data.metabolites=data.metabolites/2;
data.fluxes=data.fluxes/2;

% derive AXP dynamics from exp data
AXP=sum(data.nucleotides(:,2:4)');
time=[0:1:340];
AXP_intp=interp1(data.time_nucleotides,AXP,time,'linear','extrap');
dAXPdt=diff(AXP_intp);
data.time_dAXPdt=time(2:end);
data.AXP=AXP;
data.dAXPdt=dAXPdt;

vHeerden_GP = data;

clear AXP time AXP_intp dAXPdt data ans
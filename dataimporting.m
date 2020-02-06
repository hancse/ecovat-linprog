% This script loads data from excel sheet and stores it into mat files
% This saves run time, loading data from excel is slow

%% Irradiation in J/cm2:
Irr_J_per_cm2 = xlsread('Climate_Data',3,'BA11:BA8770');
save('Irr_J_per_cm2.mat','Irr_J_per_cm2');

%% Heat demand:
tabwater = xlsread('Climate_Data',2,'K11:K8770');
save('tabwater.mat','tabwater');
space_heating = xlsread('Climate_Data',2,'O11:O8770');
save('space_heating.mat','space_heating');

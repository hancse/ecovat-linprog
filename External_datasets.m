%% Load and pre-process extrenal datasets:
% This script loads the external datasets used in ECOVAT model
% Hourly samples for a duration of one year (8760 samples)

%% Solar irradiation [W/m2]:
disp('> Loading Solar Irradiation Dataset...')
G = xlsread('Climate_Data',3,'BA11:BA8770');
G = G.*((10^4)/(3600)); %convert to [w/m2]
disp('> Loading Solar Irradiation Dataset Completed.')

%% Ambient Temperature [C]:
disp('> Loading Ambient Temperature Dataset...')
load('NEN5060_A2');
Ta = Tout(:,2);
disp('> Loading Solar Irradiation Dataset Completed.')

%% Heat Demand [W] [Wh]:
Nhouses = input('Enter the Number of Households Connected to ECOVAT:');
Aoffices= input('Enter the Area (in m2) of the Office Space:');
disp('> Generating Heat demand Dataset...')
load('space_heating.mat');
load('tabwater.mat');
% Household heating demand per household [GJ/household]
HHr= 9.5;
% Hot tapwater demand per household [GJ/household]
HTr= 12.52;
% Office space heating demand per m2 [GJ/m2]
Ofr = 0.1143;
% Normalized values:
Norm_SH = space_heating./sum(space_heating);
Norm_TW = tabwater./sum(tabwater);
% Datasets in [GJ]:
Q_heating = (Nhouses*HHr + Ofr*Aoffices).*Norm_SH;
Q_tab = (Nhouses*HTr).*Norm_TW;
%Total heat demand [Wh]
Q_dem = (Q_heating +Q_tab).*((10^9)/3600);
disp('> Heat demand Dataset Generated')

% The temperature at at which the demand should be supplied [C]:
Tdem =40;

% Total number of PVTs:
Npvt =Nhouses*6;
%% Electricity prices [Euro/MWh]:
disp('> Loading Electricity Prices Dataset...')
load('HM_fEAPX');
ePrice = fEAPX;
disp('> Loading Electricity Prices Dataset Completed.')
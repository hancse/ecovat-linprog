%% Problem notation and parameters:
disp('> Loading ECOVAT System Parameters...')

% Time horizon (Hourly for one year period):
T_horizon = 24*4;
% Time step (15 minutes)
dt = 900; %[seconds]
% Number of segments (1= highest, 5=lowest):
Nseg = 5; 
%Number of devices:
Ndev =6;
% Heat Sources/load:
Dev = ["PVT","AW","WW1","WW2","RES","DEM"];

% Maximum temperature in each segment:
Tmax =[90 90 90 90 90]; %[C]
% Masses of segments 1 to 5 [Kg]:
Ms = [1.04e6 1.04e6 1.04e6 9.11e5 9.11e5];
% Area of a PVT panel [m2]:
Apvt = 1.8;
% Specific Heat capacity of water [J/(kg K)]:
Cp= 4168;
% Flow rate of water through the PVT [Kg/s]* 3600;
Fpvt = 0.018*3600;
% Thermal efficiency at a reduced temperature of zero [m2??CW]:
nth0 = 0.73;
% COP of the air/water heat pump:
COPaw = 2.686;
% COP of water/water heat pump 1:
COPww1 = 2.851;
% COP of water/water heat pump 2:
COPww2 = 3.681;
% Initial temperature in each segment [C]:
Temp0 = [95 90 90 90 90];
% PVT thermal loss coefficient:
ath = 7.25;
% PVT electrical loss coefficient:
ael = 0.44;
% maximum thermal efficiency of thePVT panels:
nth_max = 0.75;
% PVT electrical efficiency at a reduced temperature of zero:
nel0 = 0.1;
% maximum electrical efficiency of the PVT panels:
nel_max = 0.15;
% air/water heat pump minimum temperature [C]:
Taw_min = 0;
%  air/water heat pump maximum temperature [C]:
Taw_max = 59;
% water/water heat pump 1 minimum temperature [C]:
Tww1_min = 0;
% water/water heat pump 1 maximum temperature [C]:
Tww1_max = 49;
% water/water heat pump 2 minimum temperature [C]:
Tww2_min = 48;
% water/water heat pump 2 maximum temperature [C]:
Tww2_max = 79;
%Heat loss factor over 6 moths:
B = 0.08;
% Temperature of the ground water [C]:
Tgw = 15;
% Temperature at which heat is demanded [C]:
Tdem = 40;
%Big M factor for linearization of nonlinear constraints:
M =200;
% Power demanded by the AW heat pump when turned on [W]:
Caw= 9000;
% Power demanded by the WW1 heat pump when turned on [W]:
Cww1 = 1500;
% Power demanded by the WW2 heat pump when turned on [W]:
Cww2 = 1500;
% Power consumed by the resistance heater when turned on [W]:
Cres = 1000e3;
disp('> Loading ECOVAT System Parameters Completed')

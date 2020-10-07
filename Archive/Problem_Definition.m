%%
Clear all
close all 
clc

%% Load Model inputs:

% Solar irradiation [W/m2]:
G = zeros(time);
% Ambient Temperature [C]:
Ta = zeros(time);
% Heat Demand [W] [Wh]:
Qd = zeros(time);
% Electricity prices [Euro/W]:
Ke = zeros(T);

%% Load the System Parameters and initial conditions:
Model_Parameters;


%% Decision variable:
% x(t,s,d): The state of device d at time t with respect to segment s.
x = optimvar('x',{time,Nseg,Dev},'Type','integer','LowerBound',0,'UpperBound',1);

%% Model Constraints:


%% Any segment can have only one device connected to it at a time:
segmentconstr = optimconstr(time,Nseg);
for j=1:Nseg
    for i=1:time
        segmentconstr(i,j)= sum(x(i,j,:)) <=1;
    end
end
       
            
%% A device can only be connected to a single segment at a time:

deviceconstr = optimconstr(time,Ndev);
for j=1:Ndev
    for i=1:time
        deviceconstr(i,k) = sum(x(i,:,K)) <= 1;
    end
end

%% Temperature in each segment should not exceed a maximum predefined temperature:
T = optimexpr(time, Nseg);
% T(t,S) is the temperature of segment S at time t
temp_max_const = optimconstr(time,Nseg);
for j=1:Nseg
    for i=1:time
       temp_max_const(i,j)= T(i,j) <= Tmax(j);
    end
end


%% Temperature gradient constraint:
% Temperature should always be decreasing from top to bottom
temp_grad_constr = optimconstr(time,4);

for i=1:time
    temp_grad_constr(i,1) = T(i,1) > T(i,2);
    temp_grad_constr(i,2) = T(i,2) > T(i,3);
    temp_grad_constr(i,3) = T(i,3) > T(i,4);
    temp_grad_constr(i,4) = T(i,4) > T(i,5);
end
    
%% PVT Model:

%The output temperature of the PVT panel:
Tpvtout =optimexpr(T_horizon);

% Decision variable for when to connect the PVT to  the buffer:
wt = optimvar('wt',{T,Nseg},'Type','integer','LowerBound',0,'UpperBound',1);

% Consider changing into "for" loops later (Here just for clariity)
% pvtconstr1 = Tpvtout >= T(:,1).*wt(:,1) ;
% pvtconstr2 = Tpvtout >= T(:,2).*wt(:,2) ;
% pvtconstr3 = Tpvtout >= T(:,3).*wt(:,3);
% pvtconstr4 = Tpvtout >= T(:,4).*wt(:,4);
% pvtconstr5 = Tpvtout >= T(:,5).*wt(:,5);
% 
% pvtconstr11 = Tpvtout <= T(:,1) + M.*wt(:,1);
% pvtconstr22 = Tpvtout <= T(:,2) + M.*wt(:,2);
% pvtconstr33 = Tpvtout <= T(:,3) + M.*wt(:,3);
% pvtconstr44 = Tpvtout <= T(:,4) + M.*wt(:,4);
% pvtconstr55 = Tpvtout <= T(:,5) + M.*wt(:,5);







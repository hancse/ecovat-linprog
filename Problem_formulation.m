
%% Create the optimization problem :
disp('> Creating Optimization Problem...')
ECOVAT = optimproblem('ObjectiveSense','minimize');
%Integer linear programming
%Suppress iterative display

%% Define the Decision variable:
% x(t,s,d): The state of device d at time t with respect to segment s.
% 1= ON  0=OFF
x = optimvar('x',T_horizon,Nseg,Dev,'Type','integer','LowerBound',0,'UpperBound',1);

%% Problem Constraints:
disp('> Creating Constraints...')

%% Any segment can have only one device connected to it at a time:
segmentconstr = optimconstr(T_horizon,Nseg);
for j=1:Nseg
    for i=1:T_horizon
        segmentconstr(i,j)= sum(x(i,j,:)) <=1;
    end
end
ECOVAT.Constraints.segmentconstr= segmentconstr;
  disp('>> Segment Constraints created...')          
%% A device can only be connected to a single segment at a time:

deviceconstr = optimconstr(T_horizon,Ndev);
for j=1:Ndev
    for i=1:T_horizon
        deviceconstr(i,j) = sum(x(i,:,j)) <= 1;
    end
end
ECOVAT.Constraints.deviceconstr= deviceconstr;
disp('>> Devices Constraints created...') 

%% Temperature in each segment should not exceed a maximum predefined temperature:
T = optimexpr(T_horizon, Nseg);
% T(t,S) is the temperature of segment S at time t
temp_max_const = optimconstr(T_horizon,Nseg);
for j=1:Nseg
    for i=1:T_horizon
       temp_max_const(i,j)= T(i,j) <= Tmax(j);
    end
end
ECOVAT.Constraints.temp_max_const= temp_max_const;
disp('>> Temperature limits Constraints created...') 

%% Temperature gradient constraint:
% Temperature should always be decreasing from top to bottom
temp_grad_constr = optimconstr(T_horizon,4);

for i=1:T_horizon
    temp_grad_constr(i,1) = T(i,1) >= T(i,2);
    temp_grad_constr(i,2) = T(i,2) >= T(i,3);
    temp_grad_constr(i,3) = T(i,3) >= T(i,4);
    temp_grad_constr(i,4) = T(i,4) >= T(i,5);
end
ECOVAT.Constraints.temp_grad_constr= temp_grad_constr;
disp('>> Temperature gradients constraints created...') 

%% PVT Temperature output (Equality constraint) :

%The output temperature of the PVT panel:
Tpvtout =optimexpr(T_horizon);

%Equality constraint for the PVT temperature output:
PVT_temp_constr = optimconstr(T_horizon);

for i=1:T_horizon
    PVT_temp_constr(i)= Tpvtout(i)== ((2.*Fpvt.*Cp.*T(i,1))...
                                    - (ath.*Apvt.*T(i,1))...
                                    +(2.*Apvt.*nth0.*G(i))...
                                    + (2.*ath.*Apvt.*Ta(i)))...
                                    /((ath.*Apvt)+(2*Fpvt*Cp));
end 

ECOVAT.Constraints.PVT_temp_constr = PVT_temp_constr;
disp('>> PVT Temperature constraints created...') 


%% Connecting PVT to the lower segment:
% Decision variable for when to connect the PVT to  the buffer:
wt = optimvar('wt',T_horizon,'Type','integer','LowerBound',0,'UpperBound',1);
PVT_connect_constr1 = optimconstr(T_horizon);
PVT_connect_constr2 = optimconstr(T_horizon);
PVT_segment_constr = optimconstr(T_horizon);

for i=1:T_horizon
   PVT_connect_constr1(i)= Tpvtout(i) >= T(i,1)*wt(i);
   PVT_connect_constr2(i)= Tpvtout(i) <= T(i,1)+ M*wt(i);
   PVT_segment_constr(i) = x(i,1,'PVT') <= wt(i);
end

ECOVAT.Constraints.PVT_connect_constr1= PVT_connect_constr1;
ECOVAT.Constraints.PVT_connect_constr2= PVT_connect_constr2;
ECOVAT.Constraints.PVT_segment_constr= PVT_segment_constr;
disp('>> PVT connection constraints created...') 


%% The (reduced) temperature:
PVT_red_constr = optimconstr(T_horizon);
T_red = optimexpr(T_horizon);

for i=1:T_horizon
    if G(i)==0
       PVT_red_constr(i)= T_red(i)==0;
    else
       PVT_red_constr(i)=T_red(i)== (((T(i,5)+ Tpvtout(i))/2) - Ta(i))/G(i);
    end
end

ECOVAT.Constraints.PVT_red_constr = PVT_red_constr;
disp('>> PVT Temperature reduction constraints created...') 

%% Thermal efficiency of the PVT:
nth = optimexpr(T_horizon);
nth_constr1= optimconstr(T_horizon);
nth_constr2= optimconstr(T_horizon);
nth_constr3= optimconstr(T_horizon);


%yth is a binary indicator:
yth = optimvar('yth',T_horizon,'Type','integer','LowerBound',0,'UpperBound',1);

for i=1:T_horizon
    nth_constr1(i)= nth(i)<= nth0 -ath*T_red(i)+ M*yth(i);
    nth_constr2(i)= nth(i)<= nth_max*(1-yth(i));
    nth_constr3(i)= nth(i)>=0;
end

ECOVAT.Constraints.nth_constr1= nth_constr1;
ECOVAT.Constraints.nth_constr2= nth_constr2;
ECOVAT.Constraints.nth_constr3= nth_constr3;

disp('>> PVT thermal efficiency constraints created...') 

%% Electrical efficiency of the PVT:
nel = optimexpr(T_horizon);
nel_constr1= optimconstr(T_horizon);
nel_constr2= optimconstr(T_horizon);
nel_constr3= optimconstr(T_horizon);

%yth is a binary indicator:
yel = optimvar('yel',T_horizon,'Type','integer','LowerBound',0,'UpperBound',1);

for i=1:T_horizon
    nel_constr1(i)= nel(i)<= nth0 -ael*T_red(i)+ M*yel(i);
    nel_constr2(i)= nel(i)<= nel_max*(1-yel(i));
    nel_constr3(i)= nel(i)>=0;
end

ECOVAT.Constraints.nel_constr1= nel_constr1;
ECOVAT.Constraints.nel_constr2= nel_constr2;
ECOVAT.Constraints.nel_constr3= nel_constr3;

disp('>> PVT electrical efficiency constraints created...') 

%% Air/Water heat pump model:
AW_min_constr = optimconstr(T_horizon,5);

AW_max_constr1 = optimconstr(T_horizon,5);
AW_max_constr2 = optimconstr(T_horizon,5);
AW_max_constr3 = optimconstr(T_horizon,5);
AW_max_constr4 = optimconstr(T_horizon,5);

y=optimexpr(T_horizon,5);

for i=1:T_horizon
    for j=1:5
    AW_min_constr(i,j) = (1-x(i,j,'AW'))*M + T(i,j) >= Taw_min;
    end 
end

ECOVAT.Constraints.AW_min_constr = AW_min_constr;

for i=1:T_horizon 
    for j =1:5
        AW_max_constr1(i,j) = y(i,j) <= Taw_max* x(i,j,'AW');
        AW_max_constr2(i,j) = y(i,j) <= T(i,j);
        AW_max_constr3(i,j) = y(i,j) >= T(i,j) - Taw_max*(1- x(i,j,'AW'));
        AW_max_constr4(i,j) = y(i,j) >=0;
    end
end

ECOVAT.Constraints.AW_max_constr1 = AW_max_constr1;
ECOVAT.Constraints.AW_max_constr2 = AW_max_constr2;
ECOVAT.Constraints.AW_max_constr3 = AW_max_constr3;
ECOVAT.Constraints.AW_max_constr4 = AW_max_constr4;

disp('>> A/W heat pump constraints created...') 

%% Water Water heat pump 1 Model:

% Decision variable for when connecting segment s to the source of HP 1
z1_src = optimvar('z1_src',T_horizon,Nseg,'Type','integer','LowerBound',0,'UpperBound',1);

% Decision variable for when connecting segment s to the sink of HP 1
z1_sink = optimvar('z1_sink',T_horizon,Nseg,'Type','integer','LowerBound',0,'UpperBound',1);

WW1_range_constr1 = optimconstr(T_horizon,5);
WW1_range_constr2 = optimconstr(T_horizon,5);
WW1_range_constr3=  optimconstr(T_horizon,5);
WW1_range_constr4 = optimconstr(T_horizon,5);
WW1_range_constr5 = optimconstr(T_horizon,5);


for i=1:T_horizon
    for j=1:5
    WW1_range_constr1(i,j) = (1-z1_src(i,j))*M + T(i,j) >= Tww1_min;
  
    WW1_range_constr2(i,j) = y(i,j) <= Tww1_max* z1_sink(i,j);
    WW1_range_constr3(i,j) = y(i,j) <= T(i,j);
    WW1_range_constr4(i,j) = y(i,j) >= T(i,j) - Tww1_max*(1- z1_sink(i,j));
    WW1_range_constr5(i,j) = y(i,j) >=0;
    end
end

ECOVAT.Constraints.WW1_range_constr1= WW1_range_constr1;
ECOVAT.Constraints.WW1_range_constr2= WW1_range_constr2;
ECOVAT.Constraints.WW1_range_constr3= WW1_range_constr3;
ECOVAT.Constraints.WW1_range_constr4= WW1_range_constr4;
ECOVAT.Constraints.WW1_range_constr5= WW1_range_constr5;

disp('>> W/W 1 heat pump constraints created...') 


%% Water/Water Heat pump 2 Model:

% Decision variable for when connecting segment s to the source of HP 1
z2_src = optimvar('z2_src',T_horizon,Nseg,'Type','integer','LowerBound',0,'UpperBound',1);

% Decision variable for when connecting segment s to the source of HP 1
z2_sink = optimvar('z2_sink',T_horizon,Nseg,'Type','integer','LowerBound',0,'UpperBound',1);

WW2_range_constr1 = optimconstr(T_horizon,5);
WW2_range_constr2 = optimconstr(T_horizon,5);
WW2_range_constr3 = optimconstr(T_horizon,5);
WW2_range_constr4 = optimconstr(T_horizon,5);
WW2_range_constr5 = optimconstr(T_horizon,5);


for i=1:T_horizon
    for j=1:5
    WW2_range_constr1(i,j) = (1-z2_src(i,j))*M + T(i,j) >= Tww2_min;
  
    WW2_range_constr2(i,j) = y(i,j) <= Tww2_max* z2_sink(i,j);
    WW2_range_constr3(i,j) = y(i,j) <= T(i,j);
    WW2_range_constr4(i,j) = y(i,j) >= T(i,j) - Tww2_max*(1- z2_sink(i,j));
    WW2_range_constr5(i,j) = y(i,j) >=0;
    end
end

ECOVAT.Constraints.WW2_range_constr1= WW2_range_constr1;
ECOVAT.Constraints.WW2_range_constr2= WW2_range_constr2;
ECOVAT.Constraints.WW2_range_constr3= WW2_range_constr3;
ECOVAT.Constraints.WW2_range_constr4= WW2_range_constr4;
ECOVAT.Constraints.WW2_range_constr5= WW2_range_constr5;

disp('>> W/W 2 heat pump constraints created...') 


%% Heat Demand Model:
Demand_constr1 = optimconstr(T_horizon,5);
Demand_constr2 = optimconstr(T_horizon,5);
Demand_constr3 = optimconstr(T_horizon,5);
Demand_constr4 = optimconstr(T_horizon,5);


% Ydem is a variable introduced to linearize the nonlinear constraint
ydem= optimexpr(T_horizon,5);

for i=1:T_horizon
    for j=1:5
        Demand_constr1(i,j)= ydem(i,j) <= T(i,j)*x(i,j,'DEM');
        Demand_constr2(i,j)= ydem(i,j) <= Tdem;
        Demand_constr3(i,j)= ydem(i,j) >= Tdem - T(i,j)*(1-x(i,j,'DEM'));
        Demand_constr4(i,j)= ydem(i,j) >= 0;
    end
end

ECOVAT.Constraints.Demand_constr1= Demand_constr1;
ECOVAT.Constraints.Demand_constr2= Demand_constr2;
ECOVAT.Constraints.Demand_constr3= Demand_constr3;
ECOVAT.Constraints.Demand_constr4= Demand_constr4;

disp('>> heat demand constraints created...') 


%% Model of the Heat losses for each segment:

Qloss = optimexpr(T_horizon,Nseg);

for i=1:T_horizon
    for j=1:Nseg
        Qloss(i,j) = (1-(1-B)^(1/4380))*(T(i,j)-Tgw)*(Ms(j)*Cp/3600);
    end
end

%% Temperature evolution inside the buffer:
Temp_constr = optimconstr(T_horizon-1,Nseg);

   
for i=1:T_horizon-1
    
   % For the bottom segment:
   Temp_constr(i,5)= T(i+1,5)== T(i,5)+ (dt/(Ms(5)*Cp))*(...
               nth(i)*G(i)*Apvt*Npvt*x(i,5,'PVT')-...
               Cww1*(COPww1-1)*z1_src(i,5)-...
               Cww2*(COPww2-1)*z2_src(i,5)-...
               Q_dem(i)*x(i,1,'DEM')-...
               (1-((1-B)^(1/4380))*(T(i,5)-Tgw)*Ms(5)*Cp/3600));
               
    
    % For the rest of the segments:
    for j=1:4
     Temp_constr(i,j)= T(i+1,j)== T(i,j) +...
                    (dt/(Ms(j)*Cp))*(...
                    Caw*COPaw*x(i,j,'AW')+...
                    Cww1*COPww1*z1_sink(i,j)-...
                    Cww1*(COPww1-1)*z1_src(i,j)+...
                    Cww2*COPww2*z2_sink(i,j)-...
                    Cww2*(COPww2-1)*z2_src(i,j)+...
                    Cres*x(i,j,'RES')-...
                    Q_dem(i)*x(i,j,'DEM')-...
                    (1-((1-B)^(1/4380))*(T(i,j)-Tgw)*Ms(j)*Cp/3600));
    end
end     

ECOVAT.Constraints.Temp_constr = Temp_constr;

disp('>> Temperature evolution constraints created...') 


%% Create the objective function:
disp('> Creating The Objective Function...')
% The cost of operating the heat pumps at time t:
HP_Cost = optimexpr(T_horizon);
% The electricity generated by the PVTs at time t:
PVT_gen=optimexpr(T_horizon);
% The total cost of operation at time t:
tot_Cost= optimexpr(T_horizon);

for i=1:T_horizon
    PVT_gen(i)=nel(i)*G(i)*Apvt*Npvt;
    HP = 0;
    for j=1: Nseg
    HP = HP+((Caw*x(i,j,'AW'))+(Cww1*z1_sink(i,j))+(Cww2*z2_sink(i,j)));
    end
    HP_Cost(i) = ePrice(i)*HP;
    tot_Cost(i) = HP_Cost(i) - PVT_gen(i);
end

TOTAL_COST =sum(tot_Cost);

dispatch.Objective= TOTAL_COST;
disp('>> Objective function created...') 


%% Solve the Optimization Problem:
%Integer linear programming
%Suppress iterative display
disp('> Solving the optimization problem...') 
options = optimoptions('intlinprog','Display','iter');
[ECOVATsol,fval,exitflag,output] = solve(ECOVAT,'options',options);










%% Create the optimization problem :
disp('> Initializing the Optimization Problem...')
ECOVAT = optimproblem('ObjectiveSense','min');
%Integer linear programming
%Suppress iterative display

%% Define the Decision variable:
% x(t,s,d): The state of device d at time t with respect to segment s.
% 1= Device ON  0= Device OFF
x = optimvar('x',T_horizon,Nseg,Dev,'Type','integer','LowerBound',0,'UpperBound',1);

%% Problem Constraints:
disp('> Creating Constraints...')

%% Any segment can have only one device connected to it at a time:
segmentconstr = optimconstr(T_horizon,Nseg);
for j=1:Nseg
    for i=1:T_horizon
        segmentconstr(i,j)= x(i,j,'PVT')+x(i,j,'AW')+ x(i,j,'WW1')+ x(i,j,'RES')+x(i,j,'DEM')  <=1;
    end
end
ECOVAT.Constraints.segmentconstr= segmentconstr;
  disp('>> Segment Constraints created...')   
  
%% A device can only be connected to a single segment at a time:
% 
deviceconstr = optimconstr(T_horizon,Ndev);
for j=1:Ndev
    for i=1:T_horizon
        deviceconstr(i,j) = x(i,1,j)+ x(i,2,j)+x(i,3,j)...
                          + x(i,4,j)+ x(i,5,j)<= 1;
    end
end

ECOVAT.Constraints.deviceconstr= deviceconstr;
disp('>> Devices Constraints created...') 

%% Temperature in each segment should not exceed a maximum predefined temperature:
T = optimvar('T',T_horizon, Nseg,'Type','continuous');
% T(t,S) is the temperature of segment S at time t
% temp_max_constr = optimconstr(T_horizon,Nseg);
% for j=1:Nseg
%     for i=1:T_horizon
%        temp_max_constr(i,j)= T(i,j) <= Tmax(j);
%     end
% end
% 
% ECOVAT.Constraints.temp_max_const= temp_max_constr;
% disp('>> Maximum Temperatures Constraints created...') 

%% Temperature gradient constraint:
% %Temperature should always be decreasing from top to bottom
% temp_grad_constr = optimconstr(T_horizon,4);
% 
% for i=1:T_horizon
%     temp_grad_constr(i,1) = T(i,1) >= T(i,2);
%     temp_grad_constr(i,2) = T(i,2) >= T(i,3);
%     temp_grad_constr(i,3) = T(i,3) >= T(i,4);
%     temp_grad_constr(i,4) = T(i,4) >= T(i,5);
% end
% ECOVAT.Constraints.temp_grad_constr= temp_grad_constr;
% disp('>> Temperature gradients constraints created...') 

%% PVT Temperature output (Equality constraint) :

%The output temperature of the PVT panel:
Tpvtout =optimvar('Tpvtout',T_horizon,'Type','continuous');


%Equality constraint for the PVT temperature output:
PVT_temp_constr = optimconstr(T_horizon);

for i=1:T_horizon
    PVT_temp_constr(i)= Tpvtout(i)== ((2.*Fpvt.*Cp.*T(i,5))...
                                    - (ath.*Apvt.*T(i,5))...
                                    +(2.*Apvt.*nth0.*G(i))...
                                    + (2.*ath.*Apvt.*Ta(i)))...
                                    /((ath.*Apvt)+(2*Fpvt*Cp));
end 

ECOVAT.Constraints.PVT_temp_constr = PVT_temp_constr;
disp('>> PVT Temperature constraints created...') 


%% Connecting PVT to the lower segment:
% %Decision variable for when to connect the PVT to the buffer:
% wt = optimvar('wt',T_horizon,'Type','integer','LowerBound',0,'UpperBound',1);
% PVT_connect_constr1 = optimconstr(T_horizon);
% PVT_connect_constr2 = optimconstr(T_horizon);
% PVT_segment_constr3 = optimconstr(T_horizon);
% 
% 
% 
% for i=1:T_horizon
%    PVT_connect_constr1(i)= Tpvtout(i) >= T(i,5)+ (1-wt(i))*M;
%    PVT_connect_constr2(i)= Tpvtout(i) <= T(i,5)+ M*wt(i);
%    PVT_segment_constr3(i) = x(i,5,'PVT') <= wt(i);
% end
% 
% ECOVAT.Constraints.PVT_connect_constr1= PVT_connect_constr1;
% ECOVAT.Constraints.PVT_connect_constr2= PVT_connect_constr2;
% ECOVAT.Constraints.PVT_segment_constr3= PVT_segment_constr3;
% 
% disp('>> PVT connection constraints created...') 


%% The (reduced) temperature:
% PVT_red_constr = optimconstr(T_horizon);
% T_red = optimexpr(T_horizon);
% 
% for i=1:T_horizon
%     if G(i)==0
%        PVT_red_constr(i)= T_red(i)==0;
%     else
%        PVT_red_constr(i)=T_red(i)== (((T(i,5)+ Tpvtout(i))/2) - Ta(i))/G(i);
%     end
% end
% 
% ECOVAT.Constraints.PVT_red_constr = PVT_red_constr;
% disp('>> PVT Temperature reduction constraints created...') 

%% Thermal efficiency of the PVT:
% nth = optimexpr(T_horizon);
% nth_constr1= optimconstr(T_horizon);
% nth_constr2= optimconstr(T_horizon);
% nth_constr3= optimconstr(T_horizon);
% 
% 
% %yth is a binary indicator:
% yth = optimvar('yth',T_horizon,'Type','integer','LowerBound',0,'UpperBound',1);
% 
% for i=1:T_horizon
%     nth_constr1(i)= nth(i)<= nth0 -ath*T_red(i)+ M*yth(i);
%     nth_constr2(i)= nth(i)<= nth_max*(1-yth(i));
%     nth_constr3(i)= nth(i)>=0;
% end
% 
% ECOVAT.Constraints.nth_constr1= nth_constr1;
% ECOVAT.Constraints.nth_constr2= nth_constr2;
% ECOVAT.Constraints.nth_constr3= nth_constr3;
% 
% disp('>> PVT thermal efficiency constraints created...') 

%% Electrical efficiency of the PVT:
% nel = optimexpr(T_horizon);
% nel_constr1= optimconstr(T_horizon);
% nel_constr2= optimconstr(T_horizon);
% nel_constr3= optimconstr(T_horizon);
% 
% %yth is a binary indicator:
% yel = optimvar('yel',T_horizon,'Type','integer','LowerBound',0,'UpperBound',1);
% 
% for i=1:T_horizon
%     nel_constr1(i)= nel(i)<= nel0 -ael*T_red(i)+ M*yel(i);
%     nel_constr2(i)= nel(i)<= nel_max*(1-yel(i));
%     nel_constr3(i)= nel(i)>=0;
% end
% 
% ECOVAT.Constraints.nel_constr1= nel_constr1;
% ECOVAT.Constraints.nel_constr2= nel_constr2;
% ECOVAT.Constraints.nel_constr3= nel_constr3;
% 
% disp('>> PVT electrical efficiency constraints created...') 

%% Air/Water heat pump model:
AW_min_constr = optimconstr(T_horizon,5);

AW_max_constr1 = optimconstr(T_horizon,5);
AW_max_constr2 = optimconstr(T_horizon,5);
AW_max_constr3 = optimconstr(T_horizon,5);
AW_max_constr4 = optimconstr(T_horizon,5);
AW_max_constr5 = optimconstr(T_horizon,5);

% y_aw is a variable introduced to linearize the product of x*T:
y_aw= optimvar('y_aw',T_horizon,5,'Type','continuous');


for i=1:T_horizon
    for j=1:5
    AW_min_constr(i,j) = (1-x(i,j,'AW'))*M + T(i,j) >= Taw_min;
    end 
end

ECOVAT.Constraints.AW_min_constr = AW_min_constr;

for i=1:T_horizon 
    for j =1:5
       
        AW_max_constr1(i,j) = y_aw(i,j) <= M* x(i,j,'AW');
        AW_max_constr2(i,j) = y_aw(i,j) <= T(i,j);
        AW_max_constr3(i,j) = y_aw(i,j) >= T(i,j) - M*(1- x(i,j,'AW'));
        AW_max_constr4(i,j) = y_aw(i,j) >=0;
        AW_max_constr5(i,j) = y_aw(i,j) <= Taw_max;
    end
end

ECOVAT.Constraints.AW_min_constr = AW_min_constr;
ECOVAT.Constraints.AW_max_constr1 = AW_max_constr1;
ECOVAT.Constraints.AW_max_constr2 = AW_max_constr2;
ECOVAT.Constraints.AW_max_constr3 = AW_max_constr3;
ECOVAT.Constraints.AW_max_constr4 = AW_max_constr4;
ECOVAT.Constraints.AW_max_constr5 = AW_max_constr5;
disp('>> A/W heat pump constraints created...') 

%% Water Water heat pump 1 range Model:

% Decision variable for when connecting segment s to the source of HP 1
HP1_src = optimvar('HP1_src',T_horizon,Nseg,'Type','integer','LowerBound',0,'UpperBound',1);

%  Decision variable for when connecting segment s to the sink of HP 1
HP1_sink = optimvar('HP1_sink',T_horizon,Nseg,'Type','integer','LowerBound',0,'UpperBound',1);

% y_ww1 is a variable introduced to linearize the product of x*T:
y_ww1= optimvar('y_ww1',T_horizon,5,'Type','continuous');
% 
WW1_range_constr1 = optimconstr(T_horizon,5);
WW1_range_constr2 = optimconstr(T_horizon,5);
WW1_range_constr3 = optimconstr(T_horizon,5);
WW1_range_constr4 = optimconstr(T_horizon,5);
WW1_range_constr5 = optimconstr(T_horizon,5);
WW1_range_constr6 = optimconstr(T_horizon,5);

for i=1:T_horizon
    for j=1:5
    WW1_range_constr1(i,j) = (1-HP1_src(i,j))*M + T(i,j) >= Tww1_min;
  
    WW1_range_constr2(i,j) = y_ww1(i,j) <= M* HP1_sink(i,j);
    WW1_range_constr3(i,j) = y_ww1(i,j) <= T(i,j);
    WW1_range_constr4(i,j) = y_ww1(i,j) >= T(i,j) - M*(1- HP1_sink(i,j));
    WW1_range_constr5(i,j) = y_ww1(i,j) >=0;
    
    WW1_range_constr6(i,j) = y_ww1(i,j) <= Tww1_max;
    
    end
end

ECOVAT.Constraints.WW1_range_constr1= WW1_range_constr1;
ECOVAT.Constraints.WW1_range_constr2= WW1_range_constr2;
ECOVAT.Constraints.WW1_range_constr3= WW1_range_constr3;
ECOVAT.Constraints.WW1_range_constr4= WW1_range_constr4;
ECOVAT.Constraints.WW1_range_constr5= WW1_range_constr5;
ECOVAT.Constraints.WW1_range_constr6= WW1_range_constr6;

disp('>> W/W 1 heat pump range constraints created...')
 
 %% WW1 Heat pump segments constraints:
% To make sure the heat sink has a higher temperature than the heat source

% z_ww1 is a variable introduced to linearize the product of x*T:
z_ww1=optimvar('z_ww1',T_horizon,5,'Type','continuous');

WW1_sink_constr1 = optimconstr(T_horizon,5);
WW1_sink_constr2 = optimconstr(T_horizon,5);
WW1_sink_constr3 = optimconstr(T_horizon,5);
WW1_sink_constr4 = optimconstr(T_horizon,5);

% z_ww1 is a variable introduced to linearize the product of x*T:
r_ww1= optimvar('r_ww1',T_horizon,5,'Type','continuous');

WW1_source_constr1 = optimconstr(T_horizon,5);
WW1_source_constr2 = optimconstr(T_horizon,5);
WW1_source_constr3 = optimconstr(T_horizon,5);
WW1_source_constr4 = optimconstr(T_horizon,5);

% constraint To insure the heat sink has a higher temperature than the
% heat source:
WW1_tempdiff_constr = optimconstr(T_horizon);

for i=1:T_horizon 
    for j =1:5
        WW1_sink_constr1(i,j) = z_ww1(i,j) <= M* HP1_sink(i,j);
        WW1_sink_constr2(i,j) = z_ww1(i,j) <= T(i,j);
        WW1_sink_constr3(i,j) = z_ww1(i,j) >= T(i,j) - M*(1- HP1_sink(i,j));
        WW1_sink_constr4(i,j) = z_ww1(i,j) >=0;
        
        WW1_source_constr1(i,j) = r_ww1(i,j) <= M* HP1_src(i,j);
        WW1_source_constr2(i,j) = r_ww1(i,j) <= T(i,j);
        WW1_source_constr3(i,j) = r_ww1(i,j) >= T(i,j) - M*(1- HP1_src(i,j));
        WW1_source_constr4(i,j) = r_ww1(i,j) >=0;   
        
    end
    WW1_tempdiff_constr(i)= z_ww1(i,1)+ z_ww1(i,2) + z_ww1(i,3)+ z_ww1(i,4)+ z_ww1(i,5)...
                      >=  r_ww1(i,1)+ r_ww1(i,2) + r_ww1(i,3)+ r_ww1(i,4)+ r_ww1(i,5);
end

ECOVAT.Constraints.WW1_sink_constr1= WW1_sink_constr1;
ECOVAT.Constraints.WW1_sink_constr2= WW1_sink_constr2;
ECOVAT.Constraints.WW1_sink_constr3= WW1_sink_constr3;
ECOVAT.Constraints.WW1_sink_constr4= WW1_sink_constr4;

ECOVAT.Constraints.WW1_source_constr1= WW1_source_constr1;
ECOVAT.Constraints.WW1_source_constr2= WW1_source_constr2;
ECOVAT.Constraints.WW1_source_constr3= WW1_source_constr3;
ECOVAT.Constraints.WW1_source_constr4= WW1_source_constr4;

ECOVAT.Constraints.WW1_tempdiff_constr= WW1_tempdiff_constr;

disp('>> W/W 1 heat pump source/sink constraints created...')

%% WW1 heat pump connection constraint:

%To ensure that either both a heat source and heat sink are selected, or neither is.

WW1_conn_constr = optimconstr(T_horizon);

for i=1:T_horizon
    WW1_conn_constr(i)= HP1_sink(i,1)+ HP1_sink(i,2)+ HP1_sink(i,3)+ HP1_sink(i,4)+ HP1_sink(i,5)...
                     == HP1_src(i,1)+ HP1_src(i,2) + HP1_src(i,3)+ HP1_src(i,4)+ HP1_src(i,5);
end

ECOVAT.Constraints.WW1_conn_constr= WW1_conn_constr;
    
disp('>> W/W 1 heat  pump connection constraints created...')


%% Water/Water Heat pump 2 Model:
% 
% % Decision variable for when connecting segment s to the source of HP 2
% HP2_src = optimvar('HP2_src',T_horizon,Nseg,'Type','integer','LowerBound',0,'UpperBound',1);
% 
% % Decision variable for when connecting segment s to the sink of HP 2
% HP2_sink = optimvar('HP2_sink',T_horizon,Nseg,'Type','integer','LowerBound',0,'UpperBound',1);
% 
% % y_ww2 is a variable introduced to linearize the product of x*T:
% y_ww2=optimexpr(T_horizon,5);
% 
% WW2_range_constr1 = optimconstr(T_horizon,5);
% WW2_range_constr2 = optimconstr(T_horizon,5);
% WW2_range_constr3=  optimconstr(T_horizon,5);
% WW2_range_constr4 = optimconstr(T_horizon,5);
% WW2_range_constr5 = optimconstr(T_horizon,5);
% WW2_range_constr6 = optimconstr(T_horizon,5);
% 
% for i=1:T_horizon
%     for j=1:5
%     WW2_range_constr1(i,j) = (1-HP2_src(i,j))*M + T(i,j) >= Tww2_min;
%   
%     WW2_range_constr2(i,j) = y_ww2(i,j) <= M* HP2_sink(i,j);
%     WW2_range_constr3(i,j) = y_ww2(i,j) <= T(i,j);
%     WW2_range_constr4(i,j) = y_ww2(i,j) >= T(i,j) - M*(1- HP2_sink(i,j));
%     WW2_range_constr5(i,j) = y_ww2(i,j) >=0;
%     
%     WW2_range_constr6(i,j) = y_ww2(i,j) <= Tww2_max;
%     
%     end
% end
% 
% ECOVAT.Constraints.WW2_range_constr1= WW2_range_constr1;
% ECOVAT.Constraints.WW2_range_constr2= WW2_range_constr2;
% ECOVAT.Constraints.WW2_range_constr3= WW2_range_constr3;
% ECOVAT.Constraints.WW2_range_constr4= WW2_range_constr4;
% ECOVAT.Constraints.WW2_range_constr5= WW2_range_constr5;
% ECOVAT.Constraints.WW2_range_constr6= WW2_range_constr6;
% 
% disp('>> W/W 2 heat pump range constraints created...')
% 
% %% WW2 Heat pump segments constraints:
% % To make sure the heat sink has a higher temperature than the heat source
% 
% % z_ww2 is a variable introduced to linearize the product of x*T:
% z_ww2=optimexpr(T_horizon,5);
% 
% WW2_sink_constr1 = optimconstr(T_horizon,5);
% WW2_sink_constr2 = optimconstr(T_horizon,5);
% WW2_sink_constr3 = optimconstr(T_horizon,5);
% WW2_sink_constr4 = optimconstr(T_horizon,5);
% 
% % r_ww2 is a variable introduced to linearize the product of x*T:
% r_ww2=optimexpr(T_horizon,5);
% 
% WW2_source_constr1 = optimconstr(T_horizon,5);
% WW2_source_constr2 = optimconstr(T_horizon,5);
% WW2_source_constr3 = optimconstr(T_horizon,5);
% WW2_source_constr4 = optimconstr(T_horizon,5);
% 
% % constraint To insure the heat sink has a higher temperature than the
% % heat source:
% WW2_tempdiff_constr = optimconstr(T_horizon);
% 
% for i=1:T_horizon 
%     for j =1:5
%         WW2_sink_constr1(i,j) = z_ww2(i,j) <= M* HP2_sink(i,j);
%         WW2_sink_constr2(i,j) = z_ww2(i,j) <= T(i,j);
%         WW2_sink_constr3(i,j) = z_ww2(i,j) >= T(i,j) - M*(1- HP2_sink(i,j));
%         WW2_sink_constr4(i,j) = z_ww2(i,j) >=0;
%         
%         WW2_source_constr1(i,j) = r_ww2(i,j) <= M* HP2_src(i,j);
%         WW2_source_constr2(i,j) = r_ww2(i,j) <= T(i,j);
%         WW2_source_constr3(i,j) = r_ww2(i,j) >= T(i,j) - M*(1- HP2_src(i,j));
%         WW2_source_constr4(i,j) = r_ww2(i,j) >=0;   
%         
%     end
%     WW2_tempdiff_constr(i)= z_ww2(i,1)+ z_ww2(i,2) + z_ww2(i,3)+ z_ww2(i,4)+ z_ww2(i,5)...
%                       >=  r_ww2(i,1)+ r_ww2(i,2) + r_ww2(i,3)+ r_ww2(i,4)+ r_ww2(i,5);
% end
% 
% ECOVAT.Constraints.WW2_sink_constr1= WW2_sink_constr1;
% ECOVAT.Constraints.WW2_sink_constr2= WW2_sink_constr2;
% ECOVAT.Constraints.WW2_sink_constr3= WW2_sink_constr3;
% ECOVAT.Constraints.WW2_sink_constr4= WW2_sink_constr4;
% 
% ECOVAT.Constraints.WW2_source_constr1= WW2_source_constr1;
% ECOVAT.Constraints.WW2_source_constr2= WW2_source_constr2;
% ECOVAT.Constraints.WW2_source_constr3= WW2_source_constr3;
% ECOVAT.Constraints.WW2_source_constr4= WW2_source_constr4;
% 
% ECOVAT.Constraints.WW2_tempdiff_constr= WW2_tempdiff_constr;
% 
% disp('>> W/W 2 heat pump source/sink constraints created...')
% 
% %% WW2 heat pump connection constraint:
% %To ensure that either both a heat source and heat sink are selected, or neither is.
% 
% WW2_conn_constr = optimconstr(T_horizon);
% 
% for i=1:T_horizon
%     WW2_conn_constr(i)= HP2_sink(i,1)+ HP2_sink(i,2)+ HP2_sink(i,3)+ HP2_sink(i,4)+ HP2_sink(i,5)...
%                      == HP2_src(i,1)+ HP2_src(i,2) + HP2_src(i,3)+ HP2_src(i,4)+ HP2_src(i,5);
% end
% 
% ECOVAT.Constraints.WW2_conn_constr= WW2_conn_constr;
%     x(i,j,'AW')+
% disp('>> W/W 2 heat pump connection constraints created...')


%% Heat Demand Model:
Demand_constr1 = optimconstr(T_horizon,5);
Demand_constr2 = optimconstr(T_horizon,5);
Demand_constr3 = optimconstr(T_horizon,5);
Demand_constr4 = optimconstr(T_horizon,5);

Demand_constr5 = optimconstr(T_horizon,5);

Demand_constr6 = optimconstr(T_horizon);



% Ydem is a variable introduced to linearize the nonlinear constraint
ydem= optimvar('ydem',T_horizon,5,'Type','continuous');

for i=1:T_horizon
    for j=1:5
        Demand_constr1(i,j)= ydem(i,j) <= M*x(i,j,'DEM');
        Demand_constr2(i,j)= ydem(i,j) <= Tdem;
        Demand_constr3(i,j)= ydem(i,j) >= Tdem - M*(1-x(i,j,'DEM'));
        Demand_constr4(i,j)= ydem(i,j) >= 0;
        
        Demand_constr5(i,j)= T(i,j) >= ydem(i,j);
    end
    
    if Q_dem(i)==0
    Demand_constr6(i) =x(i,1,'DEM')+ x(i,2,'DEM')+ x(i,3,'DEM')+...
                       x(t,4,'DEM')+ x(t,5,'DEM')==0;
    else
    Demand_constr6(i) =x(i,1,'DEM')+ x(i,2,'DEM')+ x(i,3,'DEM')+...
                       x(i,4,'DEM')+ x(i,5,'DEM')==1;  
    end
end

ECOVAT.Constraints.Demand_constr1= Demand_constr1;
ECOVAT.Constraints.Demand_constr2= Demand_constr2;
ECOVAT.Constraints.Demand_constr3= Demand_constr3;
ECOVAT.Constraints.Demand_constr4= Demand_constr4;
ECOVAT.Constraints.Demand_constr5= Demand_constr5;
ECOVAT.Constraints.Demand_constr6= Demand_constr6;


disp('>> Heat demand constraints created...') 


%% Model of the Heat losses for each segment:

Qloss = optimexpr(T_horizon,Nseg);

for i=1:T_horizon
    for j=1:Nseg
        Qloss(i,j) = (1-(1-B)^(1/4380))*(T(i,j)-Tgw)*(Ms(j)*Cp/3600);
    end
end

%% Temperature evolution inside the buffer:

% 
% for i=1:T_horizon-1
%     
% %    T(i+1,5)= T(i,5)+ (dt/(Ms(5)*Cp))*(...
% %                nth0*G(i)*Apvt*Npvt*x(i,5,'PVT')-...
% %                Cww1*(COPww1-1)*HP1_src(i,5)-...
% %                Cww2*(COPww2-1)*HP2_src(i,5)-...
% %                Q_dem(i)*x(i,5,'DEM')-...
% %                (1-((1-B)^(1/4380))*(T(i,5)-Tgw)*Ms(5)*Cp/3600));
%                
%             T(i+1,5)= T(i,5)+ (dt/(Ms(5)*Cp))*(...
%                -Q_dem(i)*x(i,5,'DEM')-...
%                (1-((1-B)^(1/4380))*(T(i,5)-Tgw)*Ms(5)*Cp/3600));
%     % For the rest of the segments:
%     for j=1:4
%                
% %       T(i+1,j)= T(i,j) +...
% %                     (dt/(Ms(j)*Cp))*(...
% %                     Caw*COPaw*x(i,j,'AW')+...
% %                     Cww1*COPww1*HP1_sink(i,j)-...
% %                     Cww1*(COPww1-1)*HP1_src(i,j)+...
% %                     Cww2*COPww2*HP2_sink(i,j)-...
% %                     Cww2*(COPww2-1)*HP2_src(i,j)+...
% %                     Cres*x(i,j,'RES')-...
% %                     Q_dem(i)*x(i,j,'DEM')-...
% %                     (1-((1-B)^(1/4380))*(T(i,j)-Tgw)*Ms(j)*Cp/3600));
%                 
%                     T(i+1,j)= T(i,j) +...
%                     (dt/(Ms(j)*Cp))*(...
%                     Caw*COPaw*x(i,j,'AW')+...
%                     Cres*x(i,j,'RES')-...
%                     Q_dem(i)*x(i,j,'DEM')-...
%                     (1-((1-B)^(1/4380))*(T(i,j)-Tgw)*Ms(j)*Cp/3600));
%     end
% end
%      

Temp_constr = optimconstr(T_horizon,Nseg);

    
    Temp_constr(1,1)= T(1,1)==90;
    Temp_constr(1,2)= T(1,2)==75;
    Temp_constr(1,3)= T(1,3)==50;
    Temp_constr(1,4)= T(1,4)==30;
    Temp_constr(1,5)= T(1,5)==5;
   
   
for i=1:T_horizon-1
    
   Temp_constr(i+1,5)= T(i+1,5)== T(i,5)+ (dt/(Ms(5)*Cp))*(-...
                Cww1*(COPww1-1)*HP1_src(i,5)+...
               nth0*G(i)*Apvt*Npvt*x(i,5,'PVT')-...
               Q_dem(i)*x(i,5,'DEM'));
           
  
               
    % For the rest of the segments:
    for j=1:4
               
    Temp_constr(i+1,j)= T(i+1,j)== T(i,j) +...
                    (dt/(Ms(j)*Cp))*(+...
                    Cww1*COPww1*HP1_sink(i,j)-...
                    Cww1*(COPww1-1)*HP1_src(i,j)+...
                    Caw*COPaw*x(i,j,'AW')+...
                    Cres*x(i,j,'RES')-...
                    Q_dem(i)*x(i,j,'DEM'));
    end
end
     

ECOVAT.Constraints.Temp_constr = Temp_constr;

disp('>> Temperature evolution constraints created...') 


%% Create the objective function:
disp('> Creating The Objective Function...')
% The cost of operating ECOVAT at time t:
% OP_Cost = optimexpr(T_horizon);
% % The electricity generated by the PVTs at time t:
% PVT_gen=optimexpr(T_horizon);
% % The total cost of operation at time t:
% tot_Cost= optimexpr(T_horizon);
% 
% term2= optimexpr(T_horizon);
% term3= optimexpr(T_horizon);

for i=1:T_horizon
    PVT_gen(i)=nel0*G(i)*Apvt*Npvt;
    HP = 0;
    for j=1: Nseg
    HP = HP+((Caw*x(i,j,'AW'))+ (Cres*x(i,j,'RES')));
    end
    OP_Cost(i) = HP;
    tot_Cost(i) = (OP_Cost(i) - PVT_gen(i))*ePrice(i)*1e-6;
end

% Coefficients for tuning the cost function:
c2= 1e-5;
c3= 1e-5;

for i=1:T_horizon
term2(i) = c2.* [5 4 3 2 1]*[T(i,1); T(i,2); T(i,3); T(i,4) ; T(i,5)];
term3(i) = c3 * nth0*G(i)*Apvt*Npvt;
end

 ECOVAT.Objective= sum(sum(tot_Cost) - sum(term2) - sum(term3));
%ECOVAT.Objective= 1;
disp('>> Objective function created...') 


%% Solve the Optimization Problem:
%Integer linear programming
%Suppress iterative display
disp('> Solving the optimization problem...') 
%options = optimoptions('intlinprog','Display','iter');
% ECOVATproblem = prob2struct(ECOVAT);
% 
% ECOsol = intlinprog(ECOVATproblem)
%[ECOVATsol,fval,exitflag,output] = solve(ECOVAT,'options',options);
[ECOVATsol ,fval,exitflag,output] = solve(ECOVAT);




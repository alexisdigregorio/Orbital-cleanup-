
clc
clear
close all
% appendix 
r_earth = 6378;
mu_earth = 398600;


%% earth 

[x,y,z] = sphere;
x = x*6378;
y = y*6378;
z = z*6378;
figure(1)
surf(x,y,z)

%% Pegasus Data object 1

%Pegasus Rocket Body TLE Data
Pegasus.h = 5.4201e+04;
Pegasus.inc = 12.9994*(pi/180);
Pegasus.RAAN = 325.2850*(pi/180);
Pegasus.ecc =  0.094107;
Pegasus.w = 261.0031*(pi/180);
Pegasus.ME = 97.9;
Pegasus.EOD = 21327.86782693;
Pegasus.tsp = ((Pegasus.ME*(Pegasus.h^3/mu_earth^2))^(2/3))/(1-(Pegasus.ecc^2)); % time since prigee passage

[PegRvect, PegVvect] = PerigeeRandV(Pegasus.h, Pegasus.ecc, Pegasus.RAAN, Pegasus.inc, Pegasus.w);

%Checking to see if R and V vector match with TLE COEs
COES_1 = COEs(PegRvect, PegVvect);
Pegasus.Period = COES_1(7);
%:) COEs from the ECI r and v vector match with TLE COE data
state = [PegRvect PegVvect];
tspan = [0 Pegasus.tsp];
options = odeset('RelTol', 1e-8, 'AbsTol', 1e-8);
[time_peg,PegasusRV] = ode45(@EOM, tspan, state,options);

hold on 
plot3(PegasusRV(:,1), PegasusRV(:,2), PegasusRV(:,3));
hold on
plot3(PegasusRV(end,1), PegasusRV(end,2), PegasusRV(end,3),'b.', 'MarkerSize',10)

% poition vector at starting date: 
Peg_startR = [PegasusRV(end,1), PegasusRV(end,2), PegasusRV(end,3)]; 
Peg_startV = [PegasusRV(end,4), PegasusRV(end,5), PegasusRV(end,6)];


%% Falcon Data OBJECT 2

%Falcon 1 Rocket Body TLE Data
Falcon.h = 53984.81099;
Falcon.inc = 9.0452*(pi/180);
Falcon.RAAN = 264.0368*(pi/180);
Falcon.ecc =  0.046397;
Falcon.w = 331.4079*(pi/180);
Falcon.ME = 28.3575;
% time since perigee passage for Falcon
Falcon.tsp = ((Falcon.ME*(Falcon.h^3/mu_earth^2))^(2/3))/(1-(Falcon.ecc^2));Falcon.EOD = 21327.81008810;

%Solving for time passed since the Epoch of date
Falcon.EOD = 21327.81008810;
Fal_EOD_diff = (Pegasus.EOD-Falcon.EOD)*86400;

[FalRvect, FalVvect] = PerigeeRandV(Falcon.h, Falcon.ecc, Falcon.RAAN, Falcon.inc, Falcon.w);

%Checking to see if R and V vector match with TLE COEs
COES = COEs(FalRvect, FalVvect);
Falcon.Period = COES(7);
%:) COEs from the ECI r and v vector match with TLE COE data
state = [FalRvect FalVvect];

%Transfer time for Lambert's
t_transfer_addition = 87.6*60; %Seconds

tspan = [0 Falcon.tsp+Fal_EOD_diff+t_transfer_addition];

[time_fal,FalconRV] = ode45(@EOM, tspan, state,options);


hold on 
plot3(FalconRV(:,1), FalconRV(:,2), FalconRV(:,3));
hold on 
plot3(FalconRV(end,1), FalconRV(end,2), FalconRV(end,3), 'g.','MarkerSize',10); 


% poition vector at s: 
FalR_meet = [FalconRV(end,1), FalconRV(end,2), FalconRV(end,3)]; 
FalV_meet = [FalconRV(end,4), FalconRV(end,5), FalconRV(end,6)]; 


%% Ayame Data OBJECT 3

%Ayame 2 Debris TLE Data
Ayame.inc = 7.8540*(pi/180);
Ayame.RAAN = 111.3732*(pi/180);
Ayame.ecc =  0.4561612;
Ayame.w = 148.0010*(pi/180);
Ayame.ME = 250.5539;
Ayame.h = sqrt(mu_earth*(1+Ayame.ecc)*15674);
Ayame.tsp = ((Ayame.ME*(Ayame.h^3/mu_earth^2))^(2/3))/(1-(Ayame.ecc^2)); % time since prigee passage
Ayame.EOD = 21327.81529923;

% time from peg to Fal
Peg2Aye_EOD =  (Pegasus.EOD-Ayame.EOD)*86400;
Ayame_time = Ayame.tsp +  Peg2Aye_EOD ;% time to progpogate to to obtain position vecotr of Falcon at starting date 

[AyaRvect, AyaVvect] = PerigeeRandV(Ayame.h, Ayame.ecc, Ayame.RAAN, Ayame.inc, Ayame.w);

% Checking to see if R and V vector match with TLE COEs
COES_3 = COEs(AyaRvect, AyaVvect);
Ayame.Period = COES_3(7);
%:) COEs from the ECI r and v vector match with TLE COE data
state = [AyaRvect AyaVvect];

%Transfer time for Lambert's
t_transfer_addition = 87.6*60; %Seconds

tspan = [0 Ayame_time + t_transfer_addition];
    

[time_Aye,AyameRV] = ode45(@EOM, tspan, state, options);

hold on 
plot3(AyameRV(:,1), AyameRV(:,2), AyameRV(:,3));
hold on 
plot3(AyameRV(end,1), AyameRV(end,2), AyameRV(end,3),'k.', 'MarkerSize',10)



%% Titan Data OBJECT 4

%Titan 3C Transtage Debris TLE Data
Titan.inc = 2.6060*(pi/180);
Titan.RAAN = 303.4416*(pi/180);
Titan.ecc =  0.0046856;
Titan.w = 78.3306*(pi/180);
Titan.ME = 45.3306;
Titan.h = sqrt(mu_earth*(1+Titan.ecc)*41471);
Titan.tsp = ((Titan.ME*(Titan.h^3/mu_earth^2))^(2/3))/(1-(Titan.ecc^2)); % time since prigee passage
Titan.EOD = 21314.04823075;

% time from peg to Fal
Peg2Tit_EOD = (Pegasus.EOD-Titan.EOD)*24*60*60;
Titan_time = Titan.tsp + Peg2Tit_EOD; % time to progpogate to to obtain position vecotr of Falcon at starting date 

[TitRvect, TitVvect] = PerigeeRandV(Titan.h, Titan.ecc, Titan.RAAN, Titan.inc, Titan.w);

%Checking to see if R and V vector match with TLE COEs
COES_4 = COEs(TitRvect, TitVvect);
Titan.Period = COES_4(7);
%:) COEs from the ECI r and v vector match with TLE COE data
state = [TitRvect TitVvect];

%Transfer time for Lambert's
t_transfer_addition = 87.6*60; %Seconds

tspan = [0 Titan_time + t_transfer_addition];

% options = odeset('RelTol', 1e-8, 'AbsTol', 1e-8);

[time_Tit,TitanRV] = ode45(@EOM, tspan, state, options);

hold on 
plot3(TitanRV(:,1), TitanRV(:,2), TitanRV(:,3));
hold on
plot3(TitanRV(end,1), TitanRV(end,2), TitanRV(end,3),'c.', 'MarkerSize',10)
xlabel('x')
ylabel('y')
zlabel('z')
legend('Earth Radius', "Pegasus orbit", 'Peg pos', "Falcon orbit", 'Fal pos', "Ayame 2 Debris orbit", 'Ayame 2 Debri',"Titan 3C Debris orbit", 'Titan 3C Debris')

%% transfer from object 1 to 2: LAMBERTS
% poition an velocity vector: 
  
dtime = 87.6*60; % [seconds]
[v1_t1, v2_t1] = Lamb(Peg_startR,FalR_meet,dtime) ;

% dela v to get on transfer
dV_on = Peg_startV - v1_t1;

% delta v to get off transfer
dV_off = FalV_meet - v2_t1;

% total delta v 
deltaV1 = abs(norm(dV_on)) + abs(norm(dV_off)); 

disp('the total delta V to get from Pegasus to Flacon in 87.6 min is: [km/s]')
disp(deltaV1)

%% plot transfer 1: Lamberts for Pegasus to Falcon

tran1_R = Peg_startR; 
tran1_V = v1_t1;
COEs_tran1 = COEs(tran1_R,tran1_V); 
tran1_time = 87.6*60; 

%:) COEs from the ECI r and v vector match with TLE COE data
state = [tran1_R tran1_V];
tspan = [0 tran1_time];
options = odeset('RelTol', 1e-8, 'AbsTol', 1e-8);

[time,Tran1RV] = ode45(@EOM, tspan, state, options);

% plot transfer: 
figure (2)
[x,y,z] = sphere;
x = x*6378;
y = y*6378;
z = z*6378;
% Peg 
surf(x,y,z)
hold on 
plot3(PegasusRV(end,1), PegasusRV(end,2), PegasusRV(end,3), 'g.','markersize',10)
hold on 
plot3(PegasusRV(:,1), PegasusRV(:,2), PegasusRV(:,3));
hold on 
plot3(FalconRV(:,1), FalconRV(:,2), FalconRV(:,3));
hold on 
plot3(FalconRV(end,1), FalconRV(end,2), FalconRV(end,3), 'b.', 'markersize',10); 
hold on 
plot3(Tran1RV(:,1), Tran1RV(:,2), Tran1RV(:,3),'k');
grid on
title('transfer 1')
legend('Earth','Peg Pos','Pegasus','Falcon','Fal Pos','transfer')


% transfer from two to three
% delta t up to this point: 
% start mission on 21327.86782693-> Tuesday, Novemeber 23, 2021 at 20:49:40
% arrive at Flacon: 11-23-21 22:07:16
% leave Falcon: 11-24-21 06:14:16

%% transfer 2: Falcon to Ayame


% circularize orbit: need to get to rp
% where we are now, where rp is:
figure(3)
plot3(FalconRV(:,1), FalconRV(:,2), FalconRV(:,3));
hold on 
plot3(FalconRV(end,1), FalconRV(end,2), FalconRV(end,3), 'g.','MarkerSize',10); % after 5 orbits
hold on 
plot3(FalconRV(1,1), FalconRV(1,2), FalconRV(1,3), 'b.','MarkerSize',10); % after 5 orbits
legend('orbit','where we are','perigee')
grid on 


% FIND TIME TAKEN TO REACH PEIRGEE TO ADD TO TOTAL TIME
[FalRvect, FalVvect] = PerigeeRandV(Falcon.h, Falcon.ecc, Falcon.RAAN, Falcon.inc, Falcon.w);
COES = COEs(FalRvect, FalVvect);
Falcon.Period = COES(7);
state = [FalRvect FalVvect];

% Transfer time for Lambert's
t_transfer_addition = 87.6*60; %Seconds
tspan = [0 Falcon.tsp+Fal_EOD_diff+t_transfer_addition+(5*Falcon.Period)];
[time_transfer2,FalconRV] = ode45(@EOM, tspan, state,options);

% time taken to go from peri on starting date to location after 5 orbit
total_time = time_transfer2(end);
num_per_passed = 7*Falcon.Period;

t_left_sec = total_time - num_per_passed; 

time_left = abs(t_left_sec); %time left to reach perigee in orbit [min] 

% R and V vect at perigee:CHECK THAT THEY MATCH 
% find time left to reach perigee
start_pos = [FalconRV(end,1), FalconRV(end,2), FalconRV(end,3)];
start_vel = [FalconRV(end,4), FalconRV(end,5), FalconRV(end,6)];

state = [start_pos,start_vel];

tspan = [total_time (time_left+total_time+.001)];
[time_peri,Falcon_perRV] = ode45(@EOM, tspan, state,options);

Falcon_rp = [Falcon_perRV(end,1), Falcon_perRV(end,2), Falcon_perRV(end,3)];
Fal_rp = norm(Falcon_rp);
Falcon_vp = [Falcon_perRV(end,4), Falcon_perRV(end,5), Falcon_perRV(end,6)];
Fal_vp = norm(Falcon_vp);

% 1) circualrize at perigee 
ecc_fal = 0.046397; 
ecc_circ = 0; 

h_fal = Falcon.h;
h_circ = sqrt(Fal_rp*mu_earth);

v_circ = sqrt(mu_earth/Fal_rp); 

dV5 = sqrt((Fal_vp^2)+(v_circ^2)-(2*Fal_vp)*(v_circ)*cosd(0));

disp('delta V needed to circularize Falcon orbit:')
disp(dV5)
disp('km/s')

% 2) inc and RAAN change 
Falcon.inc = 9.0452*(pi/180);
Falcon.RAAN = 264.0368*(pi/180);

Ayame.inc = 7.8540*(pi/180);
Ayame.RAAN = 111.3732*(pi/180);

delta_RAAN = Falcon.RAAN - Ayame.RAAN;
alpha = acos(cos(Falcon.inc)*cos(Ayame.inc)+sin(Falcon.inc)*sin(Ayame.inc)*cos(delta_RAAN));  

dV_incRAAN = 2*v_circ*sin(alpha/2); 

disp('delta v to change inc and RAAN to same as Ayame is: [km/s]')
disp(dV_incRAAN)

% expand circularized orbit to have r same as rp of Ayame

rp = AyaRvect;
vp = AyaVvect;

r_init = norm(Falcon_rp); % r_circ 
v_init = v_circ; % already have 

% transfer orbit: 
rpt = r_init; % rad of circ orbit 
rat =  norm(AyaRvect); % rad of perigee for ayame
ecc_HT = (rat-rpt)/(rat+rpt); % (:
h_t = sqrt(rpt*mu_earth*(1+ecc_HT)); 
vat = h_t/rat ;
vpt = h_t/rpt ;

v_fin = sqrt(mu_earth/rat); 

dV_HT = abs(v_init-vpt)+abs(v_fin-vat);


% CEOs of new big circ orbit
BigCirc.inc = 7.8540*(pi/180);
BigCirc.RAAN = 111.3732*(pi/180);
BigCirc.w = 148.0010*(pi/180);
BigCirc.ecc = 0; 
BigCirc.h = sqrt(norm(AyaRvect)*mu_earth);
BigCirc.Period = (2*pi*norm(AyaRvect))/v_fin;

[BigCirc_Rvect, BigCirc_Vvect] = PerigeeRandV(BigCirc.h, BigCirc.ecc, BigCirc.RAAN, BigCirc.inc, BigCirc.w);

state = [BigCirc_Rvect, BigCirc_Vvect];
tspan = [0 0.5*BigCirc.Period];
    
[time_BigCirc,BigCircRV] = ode45(@EOM, tspan, state, options);

figure(4)
plot3(BigCircRV(:,1), BigCircRV(:,2), BigCircRV(:,3));
hold on
plot3(AyameRV(:,1), AyameRV(:,2), AyameRV(:,3));
grid on 
title('HT check with Ayame')
legend('large circualrized orbit','Ayame Orbit')

 




 











%% Space Debris Removal Project
%  Aero 351
%  Due Dec 3, 2021
%  Group 14
%  Carlin Sherman-Shannon, Selene Svencheko, Carlos Lopez, Alexis
%  Digregorio

clear 
close 
clc

rearth = 6378;
muearth = 398600;
options = odeset('RelTol',1e-8,'AbsTol',1e-8); % tolerances for ODE45

%% Epochs for objects

% object 1 - Pegasus R/B

EpPeg = 21327.86782693;

% object 2 - Falcon 1 R/B

EpFal = 21327.81008810;

% object 3 - Ayame 2 Deb

EpAya = 21327.81529923;

% object 4 - Titan 3C Deb

EpTit = 21314.04823075;

Eps = [EpPeg EpFal EpAya EpTit];

Epoch = max(Eps);

dtPegh = Epoch - EpPeg; % in fraction of day
dtFalh = Epoch - EpFal;
dtAyah = Epoch - EpAya;
dtTith = Epoch - EpTit;

dtPeg = dtPegh*24*60*60; % converted to seconds
dtFal = dtFalh*24*60*60;
dtAya = dtAyah*24*60*60;
dtTit = dtTith*24*60*60;

%% Plotting radius of earth 
% r = rearth;
% x = 0;
% y = 0;
% th = 0:pi/50:2*pi;
% xunit = r * cos(th) + x;
% yunit = r * sin(th) + y;
% plot(xunit, yunit);
[x,y,z] = sphere;
x = x*6378;
y = y*6378;
z = z*6378;

surf(x,y,z)

%% Pegasus Data

%Pegasus Rocket Body TLE Data
Pegasus.h = 5.4201e+04;
Pegasus.inc = 12.9994*(pi/180);
Pegasus.RAAN = 325.2850*(pi/180);
Pegasus.ecc =  0.094107;
Pegasus.w = 261.0031*(pi/180);
Pegasus.ME = 97.9;
Pegasus.tsp = ((Pegasus.ME*(Pegasus.h^3/muearth^2))^(2/3))/(1-(Pegasus.ecc^2)); % time since perigee passage

[PegRvect, PegVvect] = PerigeeRandV(Pegasus.h, Pegasus.ecc, Pegasus.RAAN, Pegasus.inc, Pegasus.w);

%Checking to see if R and V vector match with TLE COEs
COES = COEs(PegRvect, PegVvect);
Pegasus.Period = COES(7);
%:) COEs from the ECI r and v vector match with TLE COE data
state = [PegRvect PegVvect];
tspan = [0 Pegasus.tsp];

[time,PegasusRV] = ode45(@EOM, tspan, state, options);

hold on 
plot3(PegasusRV(:,1), PegasusRV(:,2), PegasusRV(:,3));
plot3(PegasusRV(end,1),PegasusRV(end,2),PegasusRV(end,3),'*')

PegRi = [PegasusRV(end,1),PegasusRV(end,2),PegasusRV(end,3)]; % Pegasus R and V vectors at epoch/start date
PegVi = [PegasusRV(end,4),PegasusRV(end,5),PegasusRV(end,6)];

%% Falcon Data

%Falcon 1 Rocket Body TLE Data
Falcon.h = 53984.81099;
Falcon.inc = 9.0452*(pi/180);
Falcon.RAAN = 264.0368*(pi/180);
Falcon.ecc =  0.046397;
Falcon.w = 331.4079*(pi/180);
Falcon.ME = 28.3575;
Falcon.tsp = ((Falcon.ME*(Falcon.h^3/muearth^2))^(2/3))/(1-(Falcon.ecc^2)); % time since perigee passage

[FalRvect, FalVvect] = PerigeeRandV(Falcon.h, Falcon.ecc, Falcon.RAAN, Falcon.inc, Falcon.w);

%Checking to see if R and V vector match with TLE COEs
COES = COEs(FalRvect, FalVvect);
Falcon.Period = COES(7);
% :) COEs from the ECI r and v vector match with TLE COE data
state = [FalRvect FalVvect];
tspan = [0 Falcon.tsp+dtFal];

[time,FalconRV] = ode45(@EOM, tspan, state, options);

hold on 
plot3(FalconRV(:,1), FalconRV(:,2), FalconRV(:,3));
plot3(FalconRV(end,1),FalconRV(end,2),FalconRV(end,3),'*')

FalRi = [FalconRV(end,1),FalconRV(end,2),FalconRV(end,3)]; % Falcon R and V vectors at epoch/start date
FalVi = [FalconRV(end,4),FalconRV(end,5),FalconRV(end,6)];

%% Ayame Data

%Ayame 2 Debris TLE Data
Ayame.inc = 7.8540*(pi/180);
Ayame.RAAN = 111.3732*(pi/180);
Ayame.ecc =  0.4561612;
Ayame.w = 148.0010*(pi/180);
Ayame.ME = 250.5539;
Ayame.h = sqrt(muearth*(1+Ayame.ecc)*15674);
Ayame.tsp = ((Ayame.ME*(Ayame.h^3/muearth^2))^(2/3))/(1-(Ayame.ecc^2)); % time since perigee passage

[AyaRvect, AyaVvect] = PerigeeRandV(Ayame.h, Ayame.ecc, Ayame.RAAN, Ayame.inc, Ayame.w);

%Checking to see if R and V vector match with TLE COEs
ACOES = COEs(AyaRvect, AyaVvect);
Ayame.Period = ACOES(7);
%:) COEs from the ECI r and v vector match with TLE COE data
state = [AyaRvect AyaVvect];
tspan = [0 Ayame.tsp+dtAya];

[time,AyameRV] = ode45(@EOM, tspan, state, options);

hold on 
plot3(AyameRV(:,1), AyameRV(:,2), AyameRV(:,3));
plot3(AyameRV(end,1),AyameRV(end,2),AyameRV(end,3),'*')

AyaRi = [AyameRV(end,1),AyameRV(end,2),AyameRV(end,3)]; % Ayame R and V vectors at epoch/start date
AyaVi = [AyameRV(end,4),AyameRV(end,5),AyameRV(end,6)];

%% Titan Data

%Titan 3C Transtage Debris TLE Data
Titan.inc = 2.6060*(pi/180);
Titan.RAAN = 303.4416*(pi/180);
Titan.ecc =  0.0046856;
Titan.w = 78.3306*(pi/180);
Titan.ME = 45.3306;
Titan.h = sqrt(muearth*(1+Titan.ecc)*41471);
Titan.tsp = ((Titan.ME*(Titan.h^3/muearth^2))^(2/3))/(1-(Titan.ecc^2)); % time since perigee passage

[TitRvect, TitVvect] = PerigeeRandV(Titan.h, Titan.ecc, Titan.RAAN, Titan.inc, Titan.w);

%Checking to see if R and V vector match with TLE COEs
COES = COEs(TitRvect, TitVvect);
Titan.Period = COES(7);
%:) COEs from the ECI r and v vector match with TLE COE data
state = [TitRvect TitVvect];
tspan = [0 Titan.tsp+dtTit];

[time,TitanRV] = ode45(@EOM, tspan, state, options);

hold on 
plot3(TitanRV(:,1), TitanRV(:,2), TitanRV(:,3));
plot3(TitanRV(end,1),TitanRV(end,2),TitanRV(end,3),'*')
xlabel('x')
ylabel('y')
zlabel('z')
legend('Earth Radius', "Pegasus R/B orbit", "Pegasus R/B at Epoch", "Falcon 1 R/B orbit", "Falcon 1 R/B at Epoch", "Ayame 2 Debris orbit", "Ayame 2 Debris at Epoch", "Titan 3C Debris orbit", "Titan 3c Debris at Epoch")

TitRi = [TitanRV(end,1),TitanRV(end,2),TitanRV(end,3)]; % Titan R and V vectors at epoch/start date
TitVi = [TitanRV(end,4),TitanRV(end,5),TitanRV(end,6)];

%% Pegasus R/B to Falcon 1 debris

% try lamberts first

dt1 = 87.6*60; % transfer time - 87.6 minutes converted to seconds

% find Falcon's position at the end of transfer time
state = [FalRvect FalVvect];
tspan = [0 Falcon.tsp+dtFal+dt1];

[time,FalconRV1] = ode45(@EOM, tspan, state, options);

FalR1 = [FalconRV1(end,1),FalconRV1(end,2),FalconRV1(end,3)]; % Falcon R and V vectors at end time of transfer 1
FalV1 = [FalconRV1(end,4),FalconRV1(end,5),FalconRV1(end,6)];

R1 = PegRi; % position of Pegasus at t = 0
R2 = FalR1; % position of Falcon at t = 87.6 [min]
V1 = PegVi; % velocity of Pegasus at t = 0 
V2 = FalV1; % velocity of Falcon at t = 87.6

[V1t,V2t] = LambertswSC(R1,R2,dt1); % lamberts function
% outputs transfer orbit velocities at t = 0 and t = 87.6

v1 = norm(V1);
v2 = norm(V2);
v1t = norm(V1t);
v2t = norm(V2t);

deltaV1 = abs(v1-v1t);
deltaV2 = abs(v2-v2t);

deltaV = deltaV1 + deltaV2; 
%deltav1 = norm(deltaV); % total delta V for first transfer, from object 1 (Pegasus) to object 2 (Falcon)

% transfer orbit
t1COEs = COEs(R2,V2t);

[t1Rvect, t1Vvect] = PerigeeRandV(t1COEs(1),t1COEs(2),t1COEs(4),t1COEs(3),t1COEs(6));

state = [t1Rvect t1Vvect];
tspan = [0 dt1];

[time,LambRV] = ode45(@EOM, tspan, state,options);

figure(2)
plot3(PegasusRV(:,1), PegasusRV(:,2), PegasusRV(:,3));
hold on
plot3(PegasusRV(end,1), PegasusRV(end,2), PegasusRV(end,3),'g.', 'MarkerSize',10)

hold on
plot3(FalconRV(:,1), FalconRV(:,2), FalconRV(:,3));
hold on
plot3(FalconRV(end,1), FalconRV(end,2), FalconRV(end,3),'b.', 'MarkerSize',10)

hold on
plot3(LambRV(:,1), LambRV(:,2), LambRV(:,3));
hold on
plot3(LambRV(end,1), LambRV(end,2), LambRV(end,3),'b.', 'MarkerSize',10)
grid on
grid minor
title("Transfer orbit from Pegasus R/B to Falcon R/B")
xlabel('x')
ylabel('y')
zlabel('z')
legend("Pegasus R/B Orbit", "Pegasus R/B ", "Falcon 1 R/B Orbit", "Falcon 1 R/B", "Lambert's Maneuver")

% find Ayame and Titan's position at end of first transfer
state = [AyaRvect AyaVvect];
tspan = [0 Ayame.tsp+dtAya+dt1];

[time,AyameRV1] = ode45(@EOM, tspan, state, options);

AyaR1 = [AyameRV1(end,1),AyameRV1(end,2),AyameRV1(end,3)]; % Ayame R and V vectors at end of first transfer
AyaV1 = [AyameRV1(end,4),AyameRV1(end,5),AyameRV1(end,6)];

state = [TitRvect TitVvect];
tspan = [0 Titan.tsp+dtTit+dt1];

[time,TitanRV1] = ode45(@EOM, tspan, state, options);

TitR1 = [TitanRV1(end,1),TitanRV1(end,2),TitanRV1(end,3)]; % Titan R and V vectors at end of first transfer
TitV1 = [TitanRV1(end,4),TitanRV1(end,5),TitanRV1(end,6)];

%% Falcon 1 R/B to Ayame 2 debris

% Falcon's location is the same as at the end of the first transfer
% (exactly 5 periods later)

% Falcon goes to perigee (+ 21.89 minutes)

state = [FalRvect FalVvect];
tspan = [0 Falcon.tsp+dtFal+dt1+(5*Falcon.Period)+(21.892*60)+0.001];

[time,FalconRV2] = ode45(@EOM, tspan, state, options);
FalR2p = [FalconRV2(end,1),FalconRV2(end,2),FalconRV2(end,3)]; % Falcon R and V vectors at perigee after 5 periods after transfer 1
FalV2p = [FalconRV2(end,4),FalconRV2(end,5),FalconRV2(end,6)]; 
Falr2p = norm(FalR2p);
Falv2p = norm(FalV2p);

% circularize Falcon orbit at perigee
ecc1 = Falcon.ecc;
Falv2pperp = Falcon.h/Falr2p;
Falv2pr = 0;
flightpath1 = atand(Falv2pr/Falv2pperp);
ecc2 = 0;
h2 = sqrt(Falr2p*muearth*(1+ecc2*cos(0)));
vperp2 = h2/Falr2p;
vr2 = (muearth/h2)*ecc2*sin(0);
v2 = sqrt(vr2^2 + vperp2^2);
flightpath2 = atand(vr2/vperp2);
deltafp = flightpath2-flightpath1;
deltav = sqrt(Falv2p^2+v2^2-2*Falv2p*v2*cosd(deltafp));

Falr2 = Falr2p;
Falv2 = v2;
Falv2perp = Falv2;
Falv2r = 0;

% Falcon 1 Rocket Body TLE Data for circularized orbit
CircFalcon.h = h2;
CircFalcon.inc = 9.0452*(pi/180);
CircFalcon.RAAN = 264.0368*(pi/180);
CircFalcon.ecc =  0;
CircFalcon.w = 331.4079*(pi/180);

[CircFalRvect, CircFalVvect] = PerigeeRandV(CircFalcon.h, CircFalcon.ecc, CircFalcon.RAAN, CircFalcon.inc, CircFalcon.w);

COES = COEs(CircFalRvect, CircFalVvect);
CircFalcon.Period = COES(7);

deltav2 = norm(CircFalVvect-FalVvect);

% inc and RAAN change
inci = Falcon.inc;
incf = Ayame.inc;

deltaRAAN = Ayame.RAAN - Falcon.RAAN; % [rad]

alpha = acos(cos(inci)*cos(incf) + sin(inci)*sin(incf)*cos(deltaRAAN));

deltav3 = 2*Falv2*sin(alpha/2);

% Hohmann transfer to an orbit w same radius as Ayame
rp = AyaRvect;
vp = AyaVvect;

r_init = norm(CircFalRvect); % r_circ 
v_init = norm(CircFalVvect); % already have 

% transfer orbit: 
rpt = r_init; % rad of circ orbit 
rat =  norm(AyaRvect); % rad of perigee for ayame
ecc_HT = (rat-rpt)/(rat+rpt); % (:
h_t = sqrt(rpt*muearth*(1+ecc_HT)); 
vat = h_t/rat ;
vpt = h_t/rpt ;

v_fin = sqrt(muearth/rat); 

dV_HT = abs(v_init-vpt)+abs(v_fin-vat);

aHT = (rat + rpt)/2;
THT = ((2*pi)/sqrt(muearth))*aHT^(3/2); % [sec]

HTtime = (THT/2); % [sec]

hf = sqrt(rat*muearth); % h of the post-transfer circular orbit

% Ayame location after the HT
state = [AyaRvect AyaVvect];
tspan = [0 Ayame.tsp+dtAya+dt1+(5*Falcon.Period)+(21.892*60)+0.001+HTtime];
%          origin   +epoch+t1 +5 per w Falcon   + Fal @ perigee   +HT

[time,AyameRVpHT] = ode45(@EOM, tspan, state, options);

AyaRpHT = [AyameRVpHT(end,1),AyameRVpHT(end,2),AyameRVpHT(end,3)]; % Ayame R and V vectors at end of the Hohmann transfer
AyaVpHT = [AyameRVpHT(end,4),AyameRVpHT(end,5),AyameRVpHT(end,6)];

AyapHT = COEs(AyaRpHT, AyaVpHT);

% % apse line rotation 
%     % technically, the circular orbit does not have an apse line
%     % we have to make one up 
% etadeg = 41; % [deg] measured from Ayame apse line to CircFalcon "apse line"
%     % CircFalcon apse line from "perigee" (where s/c is after HT) to "apogee" (halfway around the circle)
% eta = etadeg*(pi/180); % [rad]
% a = ((Ayame.h^2)*0 - (hf^2)*Ayame.ecc*cos(eta));
% b = (-(hf^2)*Ayame.ecc*sin(eta));
% c = (hf^2) - (Ayame.h^2);
% alpha = atan(b/a);
% TA1 = alpha + acos((c/a)*cos(alpha)); % [rad]
% TA2 = TA1 - eta; % [rad]
% 
% TA1 = eta;
% TA2 = 0; % read: the orbits intersect at Ayame's perigee (where the HT was set)

% decircularize (get our s/c onto Ayame's orbit)
ecc1 = 0;
h1 = hf; % h of circular orbit
vperp1 = v_fin;
vr1 = 0; % circular orbit - no radial velocity
v1 = vperp1;
flightpath1 = atan(vr1/vperp1);
ecc2 = Ayame.ecc;
h2 = Ayame.h;
vperp2 = (muearth/h2)*(1+ecc2); % at perigee - TA = 0 -> cosTA = 1
vr2 = 0; % at perigee - no radial velocity
v2 = sqrt(vr2^2 + vperp2^2);
flightpath2 = atan(vr2/vperp2);
deltafp = flightpath2-flightpath1;
deltavdcirc = sqrt(v1^2+v2^2-2*v1*v2*cos(deltafp));


% phasing maneuver - inner transfer
T = Ayame.Period;
ecc = Ayame.ecc;
TAb = AyapHT(5)*(pi/180); % TA of Ayame at the end of the Hohmann Transfer
TAa = 0; % our s/c is at perigee
Eb = 2*atan(sqrt(((1 - ecc)/(1 + ecc))*tan(TAb/2)));
Meb = Eb - ecc*sin(Eb);
tab = (T/(2*pi))*Meb;
tba = T - tab;

Ttrans = tba; % transfer period
at = (Ttrans*sqrt(muearth)/(2*pi))^(2/3); % semi major axis of transfer orbit
rpt = norm(AyaRvect); % radius of perigee for Ayame
rat = 2*at - rpt; % :) larger than rp
ecct = (rat - rpt)/(rat + rpt);
ht = sqrt(rpt*muearth*(1 + ecct)); 

vpt = ht/rpt;
vat = ht/rat;
vpi = v2; % :) matches Ayame perigee V magnitude
vpf = vpi; % both at perigee

deltaVphase = abs(vpt - vpi) + abs(vpf - vat);




% % lambert's
% dt3 = 565*60; % transfer time - 87.6 minutes converted to seconds
% 
% % find Ayame's position at the end of transfer time
% state = [AyaRvect AyaVvect];
% tspan = [0 Ayame.tsp+dtAya+dt1+dt3];
% 
% [time,AyameRV4] = ode45(@EOM, tspan, state, options);
% 
% AyaR4 = [AyameRV4(end,1),AyameRV4(end,2),AyameRV4(end,3)]; % Ayame R and V vectors at end time of transfer 3
% AyaV4 = [AyameRV4(end,4),AyameRV4(end,5),AyameRV4(end,6)];
% 
% R1 = CircFalRvect; 
% R2 = AyaR4; 
% V1 = CircFalVvect;  
% V2 = AyaV4; 
% 
% [V1t,V2t] = LambertswSC(R1,R2,dt3); % lamberts function
% 
% v1 = norm(V1);
% v2 = norm(V2);
% v1t = norm(V1t);
% v2t = norm(V2t);
% 
% deltaV1 = abs(v1-v1t);
% deltaV2 = abs(v2-v2t);
% 
% deltaV4 = deltaV1 + deltaV2; 
% %deltav4 = norm(deltaV); % total delta V for fourth transfer
% 
% % transfer orbit
% t4COEs = COEs(R2,V2t);


%% Last Transfer - Ayame to Titan

% dt5 = ; % transfer time -
% 
% % find Falcon's position at the end of transfer time
% state = [TitRvect TitVvect];
% tspan = [0 Titan.tsp+dtTit+dt1+dt2+dt3+dt4+dt5];
% 
% [time,TitanRV5] = ode45(@EOM, tspan, state, options);
% 
% TitR5 = [TitanRV5(end,1),TitanRV5(end,2),TitanRV5(end,3)]; % Titan R and V vectors at end time of transfer 5
% TitV5 = [TitanRV5(end,4),TitanRV5(end,5),TitanRV5(end,6)];
% 
% Ri = AyaR4; 
% Rf = TitR5;
% Vi = AyaV4; 
% Vf = TitV5;  
% 
% [Vit,Vft] = LambertswSC(Ri,Rf,dt5); % lamberts function
% % outputs transfer orbit velocities at t = 0 and t = 
% 
% vi = norm(Vi);
% vf = norm(Vf);
% vit = norm(Vit);
% vft = norm(Vft);
% 
% deltaV1 = abs(vi-vit);
% deltaV2 = abs(vf-vft);
% 
% deltaV5 = deltaV1 + deltaV2; 
% 
% % transfer orbit
% t5COEs = COEs(Rf,Vft);



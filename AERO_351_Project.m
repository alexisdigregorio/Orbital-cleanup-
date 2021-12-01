clc
clear
close all

rearth = 6378;
muearth = 398600;

% %Plotting radius of earth 
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
% time since perigee passage for Pegauses 
Pegasus.tsp = ((Pegasus.ME*(Pegasus.h^3/muearth^2))^(2/3))/(1-(Pegasus.ecc^2));

Pegasus.EOD = 21327.86782693;

[PegRvect, PegVvect] = PerigeeRandV(Pegasus.h, Pegasus.ecc, Pegasus.RAAN, Pegasus.inc, Pegasus.w);

%Checking to see if R and V vector match with TLE COEs
COES = COEs(PegRvect, PegVvect);
Pegasus.Period = COES(7);
%:) COEs from the ECI r and v vector match with TLE COE data
state = [PegRvect PegVvect];
tspan = [0 Pegasus.tsp];
options = odeset('RelTol', 1e-8, 'AbsTol', 1e-8);
[time,PegasusRV] = ode45(@EOM, tspan, state,options);

hold on 
plot3(PegasusRV(:,1), PegasusRV(:,2), PegasusRV(:,3));
hold on
plot3(PegasusRV(end,1), PegasusRV(end,2), PegasusRV(end,3),'g.', 'MarkerSize',10)

%START OF TRANSFERS
LamStart = [PegasusRV(end,1), PegasusRV(end,2), PegasusRV(end,3)];
PegTransStartVvect = [PegasusRV(end,4), PegasusRV(end,5), PegasusRV(end,6)];

%% Falcon Data

%Falcon 1 Rocket Body TLE Data
Falcon.h = 53984.81099;
Falcon.inc = 9.0452*(pi/180);
Falcon.RAAN = 264.0368*(pi/180);
Falcon.ecc =  0.046397;
Falcon.w = 331.4079*(pi/180);
Falcon.ME = 28.3575;
% time since perigee passage for Falcon
Falcon.tsp = ((Falcon.ME*(Falcon.h^3/muearth^2))^(2/3))/(1-(Falcon.ecc^2));

%Solving for time passed since the Epoch of date
Falcon.EOD = 21327.81008810;
FalconDifference = (Pegasus.EOD-Falcon.EOD)*24*60*60;

[FalRvect, FalVvect] = PerigeeRandV(Falcon.h, Falcon.ecc, Falcon.RAAN, Falcon.inc, Falcon.w);

%Checking to see if R and V vector match with TLE COEs
COES = COEs(FalRvect, FalVvect);
Falcon.Period = COES(7);
%:) COEs from the ECI r and v vector match with TLE COE data
state = [FalRvect FalVvect];

%Transfer time for Lambert's
tLamb = 87.6*60; %Seconds

tspan = [0 Falcon.tsp+FalconDifference+tLamb];

[time,FalconRV] = ode45(@EOM, tspan, state,options);

hold on 
plot3(FalconRV(:,1), FalconRV(:,2), FalconRV(:,3));
hold on
plot3(FalconRV(end,1), FalconRV(end,2), FalconRV(end,3),'b.', 'MarkerSize',10)

%Getting the final position and final velocity
LamEnd = [FalconRV(end,1), FalconRV(end,2), FalconRV(end,3)];
FalTransEndVvect = [FalconRV(end,4), FalconRV(end,5), FalconRV(end,6)];
%% Ayame Data

%Ayame 2 Debris TLE Data
Ayame.inc = 7.8540*(pi/180);
Ayame.RAAN = 111.3732*(pi/180);
Ayame.ecc =  0.4561612;
Ayame.w = 148.0010*(pi/180);
Ayame.ME = 250.5539;
Ayame.h = sqrt(muearth*(1+Ayame.ecc)*15674);
% time since perigee passage for Ayame 
Ayame.tsp = ((Ayame.ME*(Ayame.h^3/muearth^2))^(2/3))/(1-(Ayame.ecc^2));

Ayame.EOD =  21327.81529923;
AyameDifference = (Pegasus.EOD-Ayame.EOD)*24*60*60;

[AyaRvect, AyaVvect] = PerigeeRandV(Ayame.h, Ayame.ecc, Ayame.RAAN, Ayame.inc, Ayame.w);

%Checking to see if R and V vector match with TLE COEs
COES = COEs(AyaRvect, AyaVvect);
Ayame.Period = COES(7);
%:) COEs from the ECI r and v vector match with TLE COE data
state = [AyaRvect AyaVvect];
tspan = [0 Ayame.tsp+AyameDifference];
[time,AyameRV] = ode45(@EOM, tspan, state,options);

hold on 
plot3(AyameRV(:,1), AyameRV(:,2), AyameRV(:,3));
hold on
plot3(AyameRV(end,1), AyameRV(end,2), AyameRV(end,3),'c.', 'MarkerSize',10)
%% Titan Data

%Titan 3C Transtage Debris TLE Data
Titan.inc = 2.6060*(pi/180);
Titan.RAAN = 303.4416*(pi/180);
Titan.ecc =  0.0046856;
Titan.w = 78.3306*(pi/180);
Titan.ME = 45.3306;
Titan.h = sqrt(muearth*(1+Titan.ecc)*41471);
% time since perigee passage for Ayame 
Titan.tsp = ((Titan.ME*(Titan.h^3/muearth^2))^(2/3))/(1-(Titan.ecc^2));

Titan.EOD = 21314.04823075;
TitanDifference = (Pegasus.EOD-Titan.EOD)*24*60*60;

[TitRvect, TitVvect] = PerigeeRandV(Titan.h, Titan.ecc, Titan.RAAN, Titan.inc, Titan.w);

%Checking to see if R and V vector match with TLE COEs
COES = COEs(TitRvect, TitVvect);
Titan.Period = COES(7);
%:) COEs from the ECI r and v vector match with TLE COE data
state = [TitRvect TitVvect];
tspan = [0 Titan.tsp+TitanDifference];

[time,TitanRV] = ode45(@EOM, tspan, state,options);

hold on 
plot3(TitanRV(:,1), TitanRV(:,2), TitanRV(:,3));
hold on
plot3(TitanRV(end,1), TitanRV(end,2), TitanRV(end,3),'b.', 'MarkerSize',10)
xlabel('x')
ylabel('y')
zlabel('z')
legend('Earth Radius', "Pegasus R/B Orbit", "Pegasus R/B ", "Falcon 1 R/B Orbit", "Falcon 1 R/B", "Ayame 2 Debris Orbit", "Ayame 2 Debris","Titan 3C Debris Orbit","Titan 3C Debris")
%% Lambert's Transfer
 
%Solving for the transfer between Pegasus and Falcon
[v1vect,v2vect] = Lamberts(LamStart, LamEnd, tLamb);

deltaVstart = PegTransStartVvect-v1vect;
deltaVend = FalTransEndVvect-v2vect;

deltaV = abs(norm(deltaVstart)) + abs(norm(deltaVend));
disp("Delta V for a Lambert's Transfer = " + double(deltaV) + "km/s")

state = [LamStart v1vect];
tspan = [0 tLamb];
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

%% Falcon to Ayame
%Initial Falcon time at the start of transfer orbit 2
Faltr2t = Falcon.tsp+FalconDifference+tLamb+(5*Falcon.Period);

tspan = [0 Faltr2t];
state = [FalRvect FalVvect]
[Ftime,FalconRV] = ode45(@EOM, tspan, state,options);

plot3(FalconRV(:,1), FalconRV(:,2), FalconRV(:,3),'b');
hold on
plot3(FalconRV(end,1), FalconRV(end,2), FalconRV(end,3),'b.', 'MarkerSize',10)

Ayatr2t = Ayame.tsp+AyameDifference+tLamb+(5*Falcon.Period);
tspan = [0 Ayatr2t];
[~,AyameRV] = ode45(@EOM, tspan, state,options);

%Circularize orbit of Falcon
timeleft2perigee = (7*Falcon.Period)-Faltr2t;

tspan = [Faltr2t 7*Falcon.Period+.001];
state = [FalconRV(end,1) FalconRV(end,2) FalconRV(end,3) FalconRV(end,4) FalconRV(end,5) FalconRV(end,6)];
[Ftime,FalconRV] = ode45(@EOM, tspan, state,options);

hold on 
plot3(FalconRV(:,1), FalconRV(:,2), FalconRV(:,3), 'g');
hold on
plot3(FalconRV(end,1), FalconRV(end,2), FalconRV(end,3),'c.', 'MarkerSize',10)
hold on 
plot3(FalRvect(1),FalRvect(2),FalRvect(3),'r.', 'MarkerSize', 10)

rpfalcon=norm(FalconRV(end,1:3));

hnew = sqrt(rpfalcon*muearth);

%Falcon 1 Rocket Body TLE Data
CircFalcon.h = hnew;
CircFalcon.inc = 9.0452*(pi/180);
CircFalcon.RAAN = 264.0368*(pi/180);
CircFalcon.ecc =  0;
CircFalcon.w = 331.4079*(pi/180);

[CircFalRvect, CircFalVvect] = PerigeeRandV(CircFalcon.h, CircFalcon.ecc, CircFalcon.RAAN, CircFalcon.inc, CircFalcon.w);

deltaVcircl = norm(CircFalVvect-FalVvect);

%Inclination change
deltainc = Ayame.inc- Falcon.inc;

v1 = norm([FalconRV(end,4) FalconRV(end,5) FalconRV(end,6)]);

deltVinc = 2*v1*sin(deltainc/2);

%RAAN Change
deltaRAAN = Ayame.RAAN-Falcon.RAAN;


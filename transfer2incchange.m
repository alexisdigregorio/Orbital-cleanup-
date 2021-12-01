%% Falcon 1 R/B to Ayame 2 debris

% Falcon's location is the same as at the end of the first transfer
% (exactly 5 periods later)

% inclination change

dihed = Falcon.inc - Ayame.inc; % dihedral angle - change in inclination required [rad]
dt2 = 385*60; % 385 minutes converted to seconds

% Falcon velocity components at end of first transfer/start of second
Falv1 = norm(FalV1);
Falr1 = norm(FalR1);
FalTA1 = acos((((Falcon.h)^2)/(muearth*Falr1) - 1)/Falcon.ecc);
Falvperp1 = (muearth/Falcon.h)*(1 + Falcon.ecc*cos(FalTA1));
Falvrad1 = (muearth/Falcon.h)*Falcon.ecc*sin(FalTA1);

% find Ayame's position at end of second transfer
state = [AyaRvect AyaVvect];
tspan = [0 Ayame.tsp+dtAya+dt1+dt2];

[time,AyameRV2] = ode45(@EOM, tspan, state, options);

AyaR2 = [AyameRV2(end,1),AyameRV2(end,2),AyameRV2(end,3)]; % Ayame R and V vectors at end of first transfer
AyaV2 = [AyameRV2(end,4),AyameRV2(end,5),AyameRV2(end,6)]; 

% Ayame velocity components at end of second transfer
Ayav2 = norm(AyaV2);
Ayar2 = norm(AyaR2);
AyaTA2 = acos((((Ayame.h)^2)/(muearth*Ayar2) - 1)/Ayame.ecc);
Ayavperp2 = (muearth/Ayame.h)*(1 + Ayame.ecc*cos(AyaTA2));
Ayavrad2 = (muearth/Ayame.h)*Ayame.ecc*sin(AyaTA2);

deltav2 = sqrt((Ayavrad2 - Falvrad1)^2 + (Falvperp1)^2 + (Ayavperp2)^2 - 2*Falvperp1*Ayavperp2*cos(dihed));

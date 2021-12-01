function [v1, v2] = Lamberts(rvect, r2vect, t)
muearth = 398600;
rearth = 6378;

r1 = norm(rvect);
r2 = norm(r2vect);

dTA = acos(dot(rvect,r2vect)/(r1*r2));

tm = (sin(dTA)/(sqrt(1-(cos(dTA)^2))));
%:) tm is positive

%Solving for A
A = (sqrt(r1*r2)*sin(dTA))/sqrt(1-cos(dTA));
%:) A does not equal 0

C = 1/2;
S = 1/6;
z = 0;

y = r1 + r2 + A*(((z*S)-1)/sqrt(C));

UV = sqrt(y/C);

deltaT = (((UV^3)*S)/sqrt(muearth)) + ((A*sqrt(y))/sqrt(muearth));

f = 1 - (((UV^2)/r1)*C);
g = deltaT - (((UV^3)/sqrt(muearth))*S);
fdot = (sqrt(muearth)/(r1*r2))*UV*((z*S)-1);
gdot = 1 - (((UV^2)/r2)*C);

%Bisection method for z

a = 4*(pi^2);
b = -4*(pi^2);
TOL = 1*10^(-8);

z = zbisection(z, a, b, r1, r2, dTA, deltaT, t, TOL);

if z > 0
S = (sqrt(z)-sin(sqrt(z)))/((sqrt(z)^3));
C = (1-cos(sqrt(z)))/z;
elseif z < 0
S = (sinh(sqrt(-z))-sqrt(-z))/(sqrt(-z)^3);
C = (cosh(sqrt(-z))-1)/(-z);
elseif z == 0
S = 1/6;
C = 1/2;
end

A = (sqrt(r1*r2)*sin(dTA))/sqrt(1-cos(dTA));
y = r1 + r2 + A*(((z*S)-1)/sqrt(C));
UV = sqrt(y/C);
deltaT = (((UV^3)*S)/sqrt(muearth)) + ((A*sqrt(y))/sqrt(muearth));

f = 1 - (((UV^2)/r1)*C);
g = deltaT - (((UV^3)/sqrt(muearth))*S);
fdot = (sqrt(muearth)/(r1*r2))*UV*((z*S)-1);
gdot = 1 - (((UV^2)/r2)*C);

v1 = (1/g).*(r2vect - (f.*rvect));
v2 = (fdot*rvect) + (gdot*v1);

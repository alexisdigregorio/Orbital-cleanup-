function [rvect,vvect] = PerigeeRandV(h, ecc, RAAN, inc, w)
muearth = 398600;
%Solving for the R and V Vector at Periapse
TA = 0;

%Determine r vector in perifocal reference frame
rvect = ((h^2/muearth)*(1/(1+ecc*cos(TA))))*[cos(TA) sin(TA) 0];

%Determine v vector in perifocal reference frame
vvect = (muearth/h)*[-sin(TA) ecc+cos(TA) 0];

%Finding the rotation matrix from ECI to PERI
Q = Euler_Angle_313(RAAN, inc, w);

%Transposing matrix to be the rotation matrix from PERI to ECI
Q = transpose(Q);

%Solving for r and v using the rotation matrix
rECI = rvect.*Q;
vECI = vvect.*Q;

rvect = [rECI(1,1) rECI(2,1) rECI(3,1)];
vvect = [vECI(1,2) vECI(2,2) vECI(3,2)];
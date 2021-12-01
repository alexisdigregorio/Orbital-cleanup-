function coes = COEs(rvect,vvect)
muearth = 398600;
%Direction
r = norm(rvect);

%Speed
v = norm(vvect);

%radial velocity
vr = (dot(rvect, vvect)/r);

%Specific angular momentum , h
hvect = cross(rvect,vvect);
h = norm(hvect);

%Inclination, inc
inc = acosd(hvect(3)/h);

%Node Line, N
%Defining unit vector k
kvect = [0 0 1];

%Cross Product and Magnitude
Nvect = cross(kvect,hvect);
N = norm(Nvect);

%Right Ascension of the Ascending Node, Ω
if Nvect(2) >=0 
    RAAN = acosd(Nvect(1)/N);
else
    RAAN = acosd(Nvect(1)/N);
    RAAN = 360-RAAN;
end

%Eccentricity vector
eccvect = (1/muearth)*(cross(vvect,hvect)-(muearth.*(rvect./r)));

%Eccentricity, ecc
ecc = norm(eccvect);

%Argument of Perigee, ω
if eccvect(3) >= 0
    w = acosd((dot(Nvect,eccvect))/(N*ecc));
else
    w = acosd((dot(Nvect,eccvect))/(N*ecc));
    w = 360 - w;
end

%True Anomoly, Θ
if vr >= 0
    theta = acosd(dot(eccvect,rvect)/(ecc*r));
elseif vr < 0
    theta = 360 - acosd(dot(eccvect,rvect)/(ecc*r));
end

rp=(h^2/muearth)*(1/(1+ecc*cosd(0)));
ra=(h^2/muearth)*(1/(1+ecc*cosd(180)));
a=(.5)*(ra+rp);
period = (2*pi/(sqrt(muearth)))*(sqrt(a))^3; %in secs

coes = [h, inc, ecc, RAAN, w, theta, period];

end
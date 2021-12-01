function [z] = zbisection(z, a, b, r1, r2, dTA, deltaT, deltaTgiven, TOL)

muearth = 398600;

while abs(deltaT - deltaTgiven) > TOL

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
     
    if deltaT < deltaTgiven
        b = z;
    else 
        a = z;
    end
    z = (a+b)/2;

end
end
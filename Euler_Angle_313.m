function [Q] = Euler_Angle_313(phi, theta, psi)
Q = zeros(3,3);

Q(1,1) = (cos(psi)*cos(phi)) - (sin(phi)*cos(theta)*sin(psi));
Q(1,2) = (cos(psi)*sin(phi)) + (cos(phi)*cos(theta)*sin(psi));
Q(1,3) = sin(psi)*sin(theta);
Q(2,1) = (-sin(psi)*cos(phi)) - (sin(phi)*cos(psi)*cos(theta));
Q(2,2) = (-sin(psi)*sin(phi)) + (cos(phi)*cos(psi)*cos(theta));
Q(2,3) = cos(psi)*sin(theta);
Q(3,1) = sin(theta)*sin(phi);
Q(3,2) = -sin(theta)*cos(phi);
Q(3,3) = cos(theta);

end
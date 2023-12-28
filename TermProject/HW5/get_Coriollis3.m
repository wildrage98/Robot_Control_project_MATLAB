function [H] = get_Coriollis3(th2,th3,dth1,dth2,dth3)
global  m2 m3  r2 r3 L1 L2

h11 = - L1*dth2^2*m3*r3*sin(th2 + th3) - L1*dth3^2*m3*r3*sin(th2 + th3) - L1*L2*dth2^2*m3*sin(th2) - L1*dth2^2*m2*r2*sin(th2) - L2*dth3^2*m3*r3*sin(th3) - 2*L1*dth1*dth2*m3*r3*sin(th2 + th3) - 2*L1*dth1*dth3*m3*r3*sin(th2 + th3) - 2*L1*dth2*dth3*m3*r3*sin(th2 + th3) - 2*L1*L2*dth1*dth2*m3*sin(th2) - 2*L1*dth1*dth2*m2*r2*sin(th2) - 2*L2*dth1*dth3*m3*r3*sin(th3) - 2*L2*dth2*dth3*m3*r3*sin(th3);
h21 = L1*dth1^2*m3*r3*sin(th2 + th3) + L1*L2*dth1^2*m3*sin(th2) + L1*dth1^2*m2*r2*sin(th2) - L2*dth3^2*m3*r3*sin(th3) - 2*L2*dth1*dth3*m3*r3*sin(th3) - 2*L2*dth2*dth3*m3*r3*sin(th3);
h31 = m3*r3*(L1*dth1^2*sin(th2 + th3) + L2*dth1^2*sin(th3) + L2*dth2^2*sin(th3) + 2*L2*dth1*dth2*sin(th3));

H = [ h11; h21; h31];

end
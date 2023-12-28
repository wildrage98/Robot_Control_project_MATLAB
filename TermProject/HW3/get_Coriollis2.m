function [H] = get_Coriollis2(th2,dth1,dth2)

global m2 L1 L2

h11 = -m2*L1*L2*sin(th2)*(2*dth1*dth2 + dth2^2);
h21 = m2*L1*L2*sin(th2)*dth1^2;

H = [ h11; h21];
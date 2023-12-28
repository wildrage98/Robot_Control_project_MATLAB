function [D] = get_Inertia2(th2)
global m1 m2 L1 L2

d11 = (m1 + m2)*L1^2 + m2*L2^2 + 2*m2*L1*L2*cos(th2);
d12 = m2*L2^2 + m2*L1*L2*cos(th2);
d21 = m2*L2^2 + m2*L1*L2*cos(th2);
d22 = m2*L2^2;

D = [d11 d12;
     d21 d22];

function dydt = TwoLink(t, y);
 
global I1
global I2
global Im1
global Im2
global m1
global m2
global Fs1
global Fs2
global Fv1
global Fv2
global L1
global L2
global r1
global r2
global g
global tau1
global tau2
 
%% State space respresentation of two link manipulator
dydt = [y(2); 
    ((I2 + cos(y(3))*L1*m2*r2)*(L1*sin(y(3))*m2*r2*y(2)^2 - tau2 + Fv2*y(4) + Fs2*sign(y(4)) - cos(y(1)+y(3))*g*m2*r2))/(I1*I2 + I1*Im2 + I2*Im1 + I2*Im2 + Im1*Im2 + I2*L1^2*m2 + Im2*L1^2*m2 - cos(y(3))^2*L1^2*m2^2*r2^2 + 2*cos(y(3))*Im2*L1*m2*r2) + ((I2 + Im2)*(tau1 - Fv1*y(2) - Fs1*sign(y(2)) + cos(y(1))*L1*g*m2 + cos(y(1))*g*m1*r1 + cos(y(1)+y(3))*g*m2*r2 + L1*sin(y(3))*y(4)*m2*r2*(y(2) + y(4)) + L1*sin(y(3))*y(2)*y(4)*m2*r2))/(I1*I2 + I1*Im2 + I2*Im1 + I2*Im2 + Im1*Im2 + I2*L1^2*m2 + Im2*L1^2*m2 - cos(y(3))^2*L1^2*m2^2*r2^2 + 2*cos(y(3))*Im2*L1*m2*r2);
    y(4);
    - ((m2*L1^2 + 2*cos(y(3))*m2*r2*L1 + I1 + I2 + Im1)*(L1*sin(y(3))*m2*r2*y(2)^2 - tau2 + Fv2*y(4) + Fs2*sign(y(4)) - cos(y(1)+y(3))*g*m2*r2))/(I1*I2 + I1*Im2 + I2*Im1 + I2*Im2 + Im1*Im2 + I2*L1^2*m2 + Im2*L1^2*m2 - cos(y(3))^2*L1^2*m2^2*r2^2 + 2*cos(y(3))*Im2*L1*m2*r2) - ((I2 + cos(y(3))*L1*m2*r2)*(tau1 - Fv1*y(2) - Fs1*sign(y(2)) + cos(y(1))*L1*g*m2 + cos(y(1))*g*m1*r1 + cos(y(1)+y(3))*g*m2*r2 + L1*sin(y(3))*y(4)*m2*r2*(y(2) + y(4)) + L1*sin(y(3))*y(2)*y(4)*m2*r2))/(I1*I2 + I1*Im2 + I2*Im1 + I2*Im2 + Im1*Im2 + I2*L1^2*m2 + Im2*L1^2*m2 - cos(y(3))^2*L1^2*m2^2*r2^2 + 2*cos(y(3))*Im2*L1*m2*r2)];

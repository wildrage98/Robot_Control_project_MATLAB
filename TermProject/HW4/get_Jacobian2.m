function [J] = get_Jacobian2(th1, th2)
global L1 L2

j11 = -L1*sin(th1) - L2*sin(th1 + th2);
j12 = -L2*sin(th1 + th2);
j21 = L1*cos(th1) + L2*cos(th1 + th2);
j22 = L2*cos(th1 + th2);
J = [j11 j12;
     j21 j22];
end
function [X] = Kinematics3(th1,th2,th3)
global L1 L2 L3

x = L1*cos(th1) + L2*cos(th1 + th2) + L2*cos(th1 + th2 + th3);
y = L1*sin(th1) + L2*sin(th1 + th2) + L2*sin(th1 + th2 + th3);

X = [x;y;0];
end
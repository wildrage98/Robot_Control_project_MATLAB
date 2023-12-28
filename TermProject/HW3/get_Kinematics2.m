function [X] = get_Kinematics2(th1,th2)
global L1 L2

x = L1*cos(th1) + L2*cos(th1 + th2);
y = L1*sin(th1) + L2*sin(th1 + th2);

X = [x;y];
end
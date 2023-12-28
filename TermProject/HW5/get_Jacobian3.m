function [J] = get_Jacobian3(th1, th2, th3)
global L1 L2 L3

j11 =  -L2 * sin(th1+th2)-L1*sin(th1)-L3*sin(th1+th2+th3);
j12 =  -L2*sin(th1+th2)-L3*sin(th1+th2+th3);
j13 =  -L3*sin(th1+th2+th3);
j21 =  L2 * cos(th1+th2)+L1*cos(th1)+L3*cos(th1+th2+th3);
j22 =  L2*cos(th1+th2)+L3*cos(th1+th2+th3);
j23 =  L3*cos(th1+th2+th3);
j31 =  1;
j32 =  1;
j33 =  1;


J = [j11 j12 j13;
     j21 j22 j23;
     j31 j32 j33];
end
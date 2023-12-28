function [G] = get_Gravity3(th1, th2, th3)
    global L1 L2 g m1 m2 m3 r1 r2 r3

 g11 = g*(m1*r1*cos(th1) + m3*r3*cos(th1 + th2 + th3) + L2*m3*cos(th1 + th2) + m2*r2*cos(th1 + th2) + L1*m2*cos(th1) + L1*m3*cos(th1));
 g21 = L2*g*m3*cos(th1 + th2) + g*m2*r2*cos(th1 + th2) + g*m3*r3*cos(th1 + th2 + th3);
 g31 = g*m3*r3*cos(th1 + th2 + th3);
 
 G = [g11;g21;g31];
end
function G = get_Gravity2(th1,th2)
global m1 m2 g L1 L2
g11 = m2*g*L2*cos(th1 + th2) + (m1 + m2)*L1*g*cos(th1);
g21 = m2*g*L2*cos(th1 + th2);
G = [ g11; g21];
function [D] = get_Inertia3(th2, th3)
    global Iz1 Iz2 Iz3 L1 L2 L3 m1 m2 m3 r1 r2 r3  
    
    d11= Iz1 + Iz2 + Iz3 - L1^2*m1 + L1^2*m2 + L1^2*m3 - L2^2*m2 + L2^2*m3 - L3^2*m3 + 2*L1*m1*r1 + 2*L2*m2*r2 + 2*L3*m3*r3 + 2*L1*m3*r3*cos(th2 + th3) + 2*L1*L2*m3*cos(th2) + 2*L1*m2*r2*cos(th2) + 2*L2*m3*r3*cos(th3);
    d12= Iz2 + Iz3 - L2^2*m2 + L2^2*m3 - L3^2*m3 + 2*L2*m2*r2 + 2*L3*m3*r3 + L1*m3*r3*cos(th2 + th3) + L1*L2*m3*cos(th2) + L1*m2*r2*cos(th2) + 2*L2*m3*r3*cos(th3);
    d13= -m3*L3^2 + 2*m3*r3*L3 + Iz3 + L1*m3*r3*cos(th2 + th3) + L2*m3*r3*cos(th3);
    d21= Iz2 + Iz3 - L2^2*m2 + L2^2*m3 - L3^2*m3 + 2*L2*m2*r2 + 2*L3*m3*r3 + L1*m3*r3*cos(th2 + th3) + L1*L2*m3*cos(th2) + L1*m2*r2*cos(th2) + 2*L2*m3*r3*cos(th3);
    d22= Iz2 + Iz3 - L2^2*m2 + L2^2*m3 - L3^2*m3 + 2*L2*m2*r2 + 2*L3*m3*r3 + 2*L2*m3*r3*cos(th3);
    d23= -m3*L3^2 + 2*m3*r3*L3 + Iz3 + L2*m3*r3*cos(th3);
    d31= -m3*L3^2 + 2*m3*r3*L3 + Iz3 + L1*m3*r3*cos(th2 + th3) + L2*m3*r3*cos(th3);
    d32= -m3*L3^2 + 2*m3*r3*L3 + Iz3 + L2*m3*r3*cos(th3);
    d33= - m3*L3^2 + 2*m3*r3*L3 + Iz3;

    D = [d11 d12 d13;
         d21 d22 d23;
         d31 d32 d33];
end
function J = calculateJ(Ixx, Iyy, Izz, L, r, m)
    J = [
        1/2 * (-Ixx + Iyy + Izz),   0,                       0,                         -m * (L - r);
        0,                          1/2 * (Ixx - Iyy + Izz), 0,                         0;
        0,                          0,                       1/2 * (Ixx + Iyy - Izz),   0;
        -m * (L - r),               0,                       0,                         m];
end
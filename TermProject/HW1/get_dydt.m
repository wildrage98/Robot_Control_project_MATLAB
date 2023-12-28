clc
close all
clear all

%% 3-DOF dydt를 얻기 위해 각 행렬을 계산, 최종적으로 dydt를 얻고 'three_link.m' 파일에 저장

% 심볼릭 변수 선언
syms L1 L2 L3 m1 m2 m3 Ic1 Ic2 Ic3 
syms Im1 Im2 Im3 r1 r2 r3 Iz1 Iz2 Iz3
syms th1 th2 th3 dth1 dth2 dth3 ddth1 ddth2 ddth3

I1xx = 0;
I1yy = Iz1;
I1zz = Iz1;

I2xx = 0;
I2yy = Iz2;
I2zz = Iz2;

I3xx = 0;
I3yy = Iz3;
I3zz = Iz3;

% DH Parameters
d1 = 0;     d2 = 0;     d3 = 0;
a1 = L1;    a2 = L2;    a3 = L3;
al1 = 0;    al2 = 0;    al3 = 0;

DH = [th1, d1, a1, al1;
      th2, d2, a2, al2;
      th3, d3, a3, al3];

% Homogeneous Transformation matrix
global T01 T02 T03 T12 T23 T13

T01 = HT(th1, d1, a1, al1);
T12 = HT(th2, d2, a2, al2);
T23 = HT(th3, d3, a3, al3);

T02 = T01*T12;
T13 = T12*T23;
T03 = T01*T12*T23;

%% U matrix

% 미분 행렬
Qr = [0 -1  0  0 ;
      1  0  0  0 ;
      0  0  0  0 ;
      0  0  0  0];

global Q1 Q2 Q3

% 1, 2, 3번 모두 회전 관절
Q1 = Qr;
Q2 = Qr;
Q3 = Qr;

global U11 U12 U13 U21 U22 U23 U31 U32 U33

% j가 i보다 크면 영행렬
U11 = Q1*T01;       % i=1 j=1
U12 = zeros(4,4);   % i=1 j=2
U13 = zeros(4,4);   % i=1 j=3

U21 = Q1*T02;       % i=2 j=1
U22 = T01*Q2*T12;   % i=2 j=2
U23 = zeros(4,4);   % i=2 j=3

U31 = Q3*T03;       % i=3 j=1
U32 = T01*Q2*T13;   % i=3 j=2
U33 = T02*Q3*T23;   % i=3 j=3

%% Pseudo-inverse

global J1 J2 J3

% calculateJ 함수 호출
J1 = calculateJ(I1xx, I1yy, I1zz, L1, r1, m1);
J2 = calculateJ(I2xx, I2yy, I2zz, L2, r2, m2);
J3 = calculateJ(I3xx, I3yy, I3zz, L3, r3, m3);
  
%% Inertia matrix

n = 3;

for i=1:n
    for k=1:n
        M(i,k) = Inertia(i,k,n); % Inertia 함수 호출
    end
end

%% dUdq

n = 3;

for i=1:n
    for j=1:n
        for k=1:n
            cmd = sprintf('U%d%d%d = dUdq(i,j,k);',i,j,k);
            eval(cmd);
        end
    end
end

%% Coriolis & Centrifugal (h_ikm Matrix)
                                    %                 dth(k)*dth(m)
h111 = trace(U111*J1*U11.')...      % i=1 k=1 m=1 j=1   
      +trace(U211*J2*U21.')...      % i=1 k=1 m=1 j=2   
      +trace(U311*J3*U31.');        % i=1 k=1 m=1 j=3   dth1*dth1
  
h112 = trace(U212*J2*U21.')...      % i=1 k=1 m=2 j=2   
      +trace(U312*J3*U31.');        % i=1 k=1 m=2 j=3   dth1*dth2

h113 = trace(U313*J3*U31.');        % i=1 k=1 m=3 j=3   dth1*dth3
  

h121 = trace(U221*J2*U21.')...      % i=1 k=2 m=1 j=2   
      +trace(U321*J3*U31.');        % i=1 k=2 m=1 j=3   dth2*dth1

h122 = trace(U222*J2*U21.')...      % i=1 k=2 m=2 j=2   
      +trace(U322*J3*U31.');        % i=1 k=2 m=2 j=3   dth2*dth2
  
h123 = trace(U323*J3*U31.');        % i=1 k=2 m=3 j=3   dth2*dth3

h131 = trace(U331*J3*U31.');        % i=1 k=3 m=1 j=3   dth3*dth1

h132 = trace(U332*J3*U31.');        % i=1 k=3 m=2 j=3   dth3*dth2

h133 = trace(U333*J3*U31.');        % i=1 k=3 m=3 j=3   dth3*dth3


h1 = (dth1*dth1)*(h111) + (dth1*dth2)*(h112) + (dth1*dth3)*(h113)...
    +(dth2*dth1)*(h121) + (dth2*dth2)*(h122) + (dth2*dth3)*(h123)...
    +(dth3*dth1)*(h131) + (dth3*dth2)*(h132) + (dth3*dth3)*(h133);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                                    %                 dth(k)*dth(m)
h211 = trace(U211*J2*U22.')...      % i=2 k=1 m=1 j=2   
      +trace(U311*J3*U32.');        % i=2 k=1 m=1 j=3   dth1*dth1

h212 = trace(U212*J2*U22.')...      % i=2 k=1 m=2 j=2
      +trace(U312*J3*U32.');        % i=2 k=1 m=2 j=3   dth1*dth2

h213 = trace(U313*J3*U32.');        % i=2 k=1 m=3 j=3   dth1*dth3
  
h221 = trace(U221*J2*U22.')...      % i=2 k=2 m=1 j=2   dth2*dth1
      +trace(U321*J3*U32.');        % i=2 k=2 m=1 j=3

h222 = trace(U222*J2*U22.')...      % i=2 k=2 m=2 j=2   dth2*dth2
      +trace(U322*J3*U32.');        % i=2 k=2 m=2 j=3

h223 = trace(U323*J3*U32.');        % i=2 k=2 m=3 j=3   dth2*dth3

h231 = trace(U331*J3*U32.');        % i=2 k=3 m=1 j=3   dth3*dth1

h232 = trace(U332*J3*U32.');        % i=2 k=3 m=2 j=3   dth3*dth2

h233 = trace(U333*J3*U32.');        % i=2 k=3 m=3 j=3   dth3*dth3

h2 = (dth1*dth1)*(h211) + (dth1*dth2)*(h212) + (dth1*dth3)*(h213)...
    +(dth2*dth1)*(h221) + (dth2*dth2)*(h222) + (dth2*dth3)*(h223)...
    +(dth3*dth1)*(h231) + (dth3*dth2)*(h232) + (dth3*dth3)*(h233);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                                    %                 dth(k)*dth(m)
h311 = trace(U311*J3*U33.');        % i=3 k=1 m=1 j=3   dth1*dth1

h312 = trace(U312*J3*U33.');        % i=3 k=1 m=2 j=3   dth1*dth2

h313 = trace(U313*J3*U33.');        % i=3 k=1 m=3 j=3   dth1*dth3

h321 = trace(U321*J3*U33.');        % i=3 k=2 m=1 j=3   dth2*dth1

h322 = trace(U322*J3*U33.');        % i=3 k=2 m=2 j=3   dth2*dth2

h323 = trace(U323*J3*U33.');        % i=3 k=2 m=3 j=3   dth2*dth3

h331 = trace(U331*J3*U33.');        % i=3 k=3 m=1 j=3   dth3*dth1

h332 = trace(U332*J3*U33.');        % i=3 k=3 m=2 j=3   dth3*dth2

h333 = trace(U333*J3*U33.');        % i=3 k=3 m=3 j=3   dth3*dth3

h3 = (dth1*dth1)*(h311) + (dth1*dth2)*(h312) + (dth1*dth3)*(h313)...
    +(dth2*dth1)*(h321) + (dth2*dth2)*(h322) + (dth2*dth3)*(h323)...
    +(dth3*dth1)*(h331) + (dth3*dth2)*(h332) + (dth3*dth3)*(h333);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 위의 각 요소를 통해 얻어진 최종 H matrix

h = simplify([h1; h2; h3]);

%% Gravity

syms g

% 무게중심까지의 거리
r11 = [-(L1-r1); 0; 0; 1];  
r22 = [-(L2-r2); 0; 0; 1];
r33 = [-(L3-r3); 0; 0; 1];

gv = [0 -g 0 0];

G1 = -( m1*gv*U11*r11...   % i=1 j=1
       +m2*gv*U21*r22...   % i=1 j=2
       +m3*gv*U31*r33);    % i=1 j=3
   
   
G2 = -( m2*gv*U22*r22...   % i=2 j=2
       +m3*gv*U32*r33);    % i=2 j=3
   
G3 = -(m3*gv*U33*r33);     % i=3 j=3

G = simplify([G1;G2;G3]);

%% DDTH

syms tau1 tau2 tau3

% [tau1; tau2; tau3] = M * [ddth1; ddth2; ddth3] + h + G
% [ddth1; ddth2; ddth3] = inv(M) * ([tau1; tau2; tau3] - h - G)

DDTH = inv(M)*([tau1; tau2; tau3] -h -G); % theta double dot에 의해 만들어진 행렬

%% dydt

dydt = simplify([dth1; DDTH(1); dth2; DDTH(2); dth3; DDTH(3)]);
% matlabFunction(dydt,'file','three_links.m','Optimize',false); % 위 dydt를 함수로 만들기

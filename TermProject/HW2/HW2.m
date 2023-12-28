clc
clear all;
close all

%% HW2: Using the following Dynamics model, estimate the 2-DOF dynamics parameter.

% 19조_2019741054_이종우_2021741071_최나은

% 구성파일 HW2, TwoLink_Regressor

%% 2-DOF regressor 

syms s dq ddq

global I1 I2 L1 L2 Im1 Im2 m1 m2 g r1 r2 Fs1 Fs2 Fv1 Fv2 tau1 tau2

I1 = 0.05;   I2 = 0.05;     % [kgm^2], 링크 관성
Im1 = 0.05;  Im2 = 0.05;    % [kgm^2], 링크 관성
L1 = 0.5;    L2 = 0.5;      % [m], 링크 길이
m1 = 0.2;    m2 = 0.2;      % [kg], 링크 질량
Fs1 = 0.1;   Fs2 = 0.1;     % [Nm]
Fv1 = 0.1;   Fv2 = 0.1;     % [Nm]
r1 = 0.1;    r2 = 0.1;      % [m], 질량 중심
g = 9.806;                  % [m/s^2], 중력 가속도


dt = 0.01; ft = 20 ;

q1 = 0; dq1 = 0;
q2 = 0; dq2 = 0; % 조인트 각과 각속도 초기화 
n = 1;

theta = [0;0;0;0;0;0;0;0;0;0]; % regressor 초기화

Yn = zeros(10); % 10X10 단위행렬
un = 0;
W1_int = [0,0,0,0,0,0,0,0,0,0;
          0,0,0,0,0,0,0,0,0,0];
W2 = [ 0 0 0 0 0 0 0 0 0 0;
       0 0 0 0 0 0 0 0 0 0];    % 칼만 필터에 필요한 파라미터 초기화
u = [0;
    0];
P = eye(10,10);
R = eye(2,2);

%% 시뮬레이션 루프

for cnt=0:dt:ft
    clc;
    disp(cnt);
    
    
    tau1 = sin(cnt) + cos(10*cnt);
    tau2 = sin(cnt) + cos(10*cnt); % 2 link 회전력의 관계식

    [st,x] = ode45('TwoLink_Regressor',[0,dt],[q1; dq1; q2; dq2]);
    index = size(x);
    
    q1  = x(index(1),1);
    dq1 = x(index(1),2);
    q2  = x(index(1),3);
    dq2 = x(index(1),4); % regressor를 구하기 위해 필요한 각도와 각속도

    % 칼만 필터 기반 알고리즘 
    W1_int = W1_int + [ 0 0                 0            g*cos(q1)   g*cos(q1 + q2)   0  -sign(dq1)      0        -dq1      0;
                        0 0    -sin(q2)*(dq2 + dq2)*dq1      0       g*cos(q1 + q2)   0      0       -sign(dq2)     0     -dq2 ] * dt;

    W2 = [ dq1 dq1 + dq2      cos(q2)*dq2+2*cos(q2)*dq1    0  0    0  0  0  0  0;
            0  dq1 + dq2          cos(q2)*(dq1)        0  0   dq2 0  0  0  0 ];

    Y = W2 - W1_int;

    u = u + [tau1; tau2]*dt; % 칼만 필터 파라미터 업데이트
    
    % 칼만 필터 
    P = P - P*Y' *inv(R + Y*P*Y')*Y*P;
    K = P*Y';
    theta = theta + K*(u-Y*theta);
    % Data save
    save_time(n,:) = cnt;      % 시간 저장
    save_theta(n,:) = theta;   % 결과 저장
    n = n + 1;
end

%% 결과 그래프 도출

clc;
close all;
theta % theta 최종결과

figure(1)
hold on
grid on
plot(save_time, save_theta)
legend({'I_1 + m_2I_1^2 + I_m_1', 'I_2', 'm_2r_2I_1', 'm_1r_1 + m_2I_1', 'm_2r_2', 'I_m_2', 'F_s_1', 'F_s_2', 'F_v_1', 'F_v_2'}, 'location', 'best')
AX2.XLabel.String = 'time';
AX2.YLabel.String = 'parameter';
AX2.FontSize = 12;
xlim([0 ft])
ylim([0 Fs1*2])
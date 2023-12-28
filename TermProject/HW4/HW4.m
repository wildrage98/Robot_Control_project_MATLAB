clc, clear all
close all

%% HW3-2: (2-DOF) Cartesian Space PID CTM Controller

% 19조_2019741054_이종우_2021741071_최나은_HW3_1.m

% 목표 궤적: 반경 0.1[m] / 주기 1[s] 원 그리기 
% 중력 보상 오차 적용

%% 시뮬레이션 설정

% 단위 변환 상수
DR = deg2rad(1);    % Degree -> Radian 
RD = rad2deg(1);    % Radian -> Degree 

% 시뮬레이션 플래그
flag_Sim = 1;
flag_Draw = 1;
flag_Draw_Robot = 1;
flag_Draw_Graph = 1;

% 시뮬레이션 파라미터
delta_t     = 0.005;        % [초] : 샘플링 시간
start_t     = 0.000;        % [초] : 시작 시간
finish_t    = 6.000;        % [초] : 종료 시간

g           = 9.8148;       % [m/s^2] : 중력 가속도

%% 로봇 파라미터

% 전역 변수
global Iz1 Iz2 L1 L2 g m1 m2 r1 r2 tau1 tau2

m1 = 0.20;   m2 = 0.20;                        % [kg]    : 링크 질량
L1 = 0.50;   L2 = 0.50;                        % [m]     : 링크 길이
r1 = 0.10;   r2 = 0.10;                        % [m]     : 질량 중심 길이
Iz1 = 0.05;  Iz2 = 0.05;                       % [kgm^2] : 링크 관성 모멘트                                 
init_q1     = -pi/2;  init_q2 = pi/2;          % [rad]   : 초기 관절 각도
init_dq1    = 0.00;   init_dq2= 0.00;          % [rad/s] : 초기 각속도
q           = [init_q1;
               init_q2];                       % [rad]   : 초기 관절 각도
dq          = [init_dq1;
               init_dq2];                      % [rad/s] : 초기 각속도
          
init_X = get_Kinematics2(q(1),q(2));           % [m]     : 초기 엔드 에펙터 위치
X      = init_X;                               % [m]     : 현재 위치
dX     = [0;0];                                % [m/s]   : 현재 속도
X_d    = init_X;                               % [m]     : 현재 엔드 에펙터 목표 위치
dX_d   = [0;0];                                % [m/s]   : 현재 엔드 에펙터 목표 속도
ddX_d  = [0;0];                                % [m/s^2] : 현재 엔드 에펙터 목표 가속도
Error       = X_d - X;
Error_sum   = 0;
tau1 = 0.0;                                    % [Nm]    : 제어 토크 1
tau2 = 0.0;                                    % [Nm]    : 제어 토크 2
tq = [tau1
      tau2];                                   % [Nm]    : 제어 토크
  
% 컨트롤러 게인
Wn          = 20;           % [rad/s]    : 고유 진동수 (Natural frequency)
Kp          = Wn^2;         % [Nm/rad]   : 비례 게인
Kv          = 2 * Wn;       % [Nm*s/rad] : 미분 게인 (감쇠비(zeta) : 1, 비례 게인에 임계 감쇠)
Ki          = 1000;         % [Nm/rad*s] : 적분 게인

%% 시뮬레이션 실행

if(flag_Sim == 1)
    n = 1;
    sin_t = 0;
    pre_J = 0;
    for (time = start_t:delta_t:finish_t)
        % 목표 궤적 설정
        if(time < 1)
            % 초기 1초 동안은 목표 위치를 초기 상태로 설정
            X_d     = init_X;
            dX_d    = [0;0];
            ddX_d   = [0;0];
        elseif(time < 2.0)
            % 1초 이상 2초 미만인 경우, 목표 위치를 초기 위치로 설정하고 Y 방향으로 0.1 만큼 이동
            X_d(1)    =  init_X(1);
            if(X_d(2) < init_X(2) + 0.1)
                X_d(2) = X_d(2) + (0.1/0.5)*delta_t;
            else
                X_d(2) = init_X(2) + 0.1;
            end
            dX_d    = (X_d  - [simul_X_d_x(n-1)  ; simul_X_d_y(n-1)])./delta_t; 
            ddX_d   = (dX_d - [simul_dX_d_x(n-1) ; simul_dX_d_y(n-1)])./delta_t;
        else
            % 2초 이상인 경우, 원 모양의 궤적을 따라 이동
            X_d     = [ 0.1*sin((2*pi*sin_t)) + init_X(1);
                        0.1*cos((2*pi*sin_t)) + init_X(2)];

            % sin_t 업데이트 및 목표 속도, 목표 가속도 설정
            sin_t   = sin_t + delta_t;
            dX_d    = (X_d  - [simul_X_d_x(n-1)  ; simul_X_d_y(n-1)])./delta_t; 
            ddX_d   = (dX_d - [simul_dX_d_x(n-1) ; simul_dX_d_y(n-1)])./delta_t;
        end

        % 동역학 얻기
        J = get_Jacobian2(q(1), q(2));
        dJ = (J-pre_J)/delta_t;
        pre_J = J;
        X = get_Kinematics2(q(1), q(2));
        dX = J*dq;
        D = get_Inertia2(q(2));
        H = get_Coriollis2(q(2), dq(1), dq(2));
        C = get_Gravity2(q(1),q(2));
        
        % PID 컨트롤러
        Error = X_d - X;
        Error_sum = Error_sum + Error;
        u = ddX_d + Kv*(dX_d - dX) + Kp*Error + Ki*Error_sum*delta_t;
        ddq_ref = inv(J)*(u-dJ*dq);
        tq_ctrl = D*ddq_ref + H + C*0.8; % 20% 중력 보상 오차

        % 로봇 모델
        % 역동역학
        tq      = tq_ctrl;
        tau1 = tq(1);    tau2 = tq(2);
        [t,y]   = ode45('two_link', [0, delta_t], [q(1); dq(1); q(2); dq(2)]); 
        index   = length(y);
        q       = [ y(index,1); y(index,3)];
        dq      = [ y(index,2); y(index,4)];

        % 데이터 저장
        simul_time(n)    = time;        % [sec]
        simul_q1(n)      = q(1);        % [rad]
        simul_q2(n)      = q(2);        % [rad]
        simul_dq1(n)     = dq(1);       % [rad/s]
        simul_dq2(n)     = dq(2);       % [rad/s]
        simul_X_x(n)     = X(1);        % [m]
        simul_X_y(n)     = X(2);        % [m]
        simul_dX_x(n)    = dX(1);       % [m/s]
        simul_dX_y(n)    = dX(2);       % [m/s]
        simul_X_d_x(n)   = X_d(1);      % [m]
        simul_X_d_y(n)   = X_d(2);      % [m]
        simul_dX_d_x(n)  = dX_d(1);     % [m/s]
        simul_dX_d_y(n)  = dX_d(2);     % [m/s]    
        n = n+1;
    end
end

%% 시뮬레이션 그래프 

if(flag_Draw == 1)
    font_size_label = 20;
    font_size_title = 25;
    lineWidth_cur   = 3;
    lineWidth_tar   = 5;

    if(flag_Draw_Robot == 1)
        % Draw Robot
            x1 = L1*cos(init_q1);              % [m]   : Joint 1 X-axis Position
            y1 = L1*sin(init_q1);              % [m]   : Joint 1 Y-axis Position
            x2 = L1*cos(init_q1+init_q2);      % [m]   : Joint 2 X-axis Position
            y2 = L1*sin(init_q1+init_q2);      % [m]   : Joint 2 Y-axis Position
            
        

        FG1 = figure('Position', [500 0 700 700], 'Color',[1 1 1]);
        AX = axes('Parent',FG1);
        hold on

        Px1 = [0 x1];       Py1 = [0 y1];     
        Px2 = [x1 x1+x2];   Py2 = [y1 y1+y2];

        p1 = plot(Px1, Py1, '-or', 'Linewidth',3);
        p2 = plot(Px2, Py2, '-ob', 'Linewidth',3);
        axis([-1.0 1.0 -1.6 0.4]);
        grid on

        xlabel('X-axis (m)','FontSize',font_size_label)
        ylabel('Y-axis (m)','FontSize',font_size_label)
        title('2-DOF Robot', 'FontSize',font_size_title)

        n = 1;
        for(time = start_t:delta_t:finish_t)
            cmd = sprintf('Time : %2.2f' , time);
            clc
            disp(cmd)
            q1 = simul_q1(n);     q2 = simul_q2(n);
            x1 = L1*cos(q1);         % [m]   : Joint 1 X-axis Position
            y1 = L1*sin(q1);         % [m]   : Joint 1 Y-axis Position
            x2 = L1*cos(q1+q2);      % [m]   : Joint 2 X-axis Position
            y2 = L1*sin(q1+q2);      % [m]   : Joint 2 Y-axis Position
            
            Px1 = [0 x1];       Py1 = [0 y1];     
            Px2 = [x1 x1+x2];   Py2 = [y1 y1+y2];
            set(p1,'XData', Px1, 'YData',Py1)
            set(p2,'XData', Px2, 'YData',Py2)
            drawnow
            n = n+1;
        end
    end
    if(flag_Draw_Graph == 1)
        % 위치 그래프
        FG2 = figure('Position', [900 100 600 500], 'Color',[1 1 1]);
        plot(simul_time,simul_X_d_x, ':r','LineWidth',lineWidth_tar); hold on;
        plot(simul_time,simul_X_d_y, ':b','LineWidth',lineWidth_tar); hold on;
        
        plot(simul_time,simul_X_x, 'r','LineWidth',lineWidth_cur); hold on;
        plot(simul_time,simul_X_y, 'b','LineWidth',lineWidth_cur); hold on;
        legend({'tar_x', 'tar_y', 'cur_x', 'cur_y'}, 'location', 'best','orientation','horizontal','fontsize',15)
        axis([start_t finish_t -1.25 1]);
        xticks([start_t:1:finish_t])
        yticks([-1:0.25:1])
        grid on

        xlabel('time (s) ','FontSize',font_size_label)
        ylabel('Position (m) ','FontSize',font_size_label)
        title('Cartesian Space PD CTM Comtroller ','FontSize',font_size_title)


        % 속도 그래프
        FG3 = figure('Position', [900 200 600 500], 'Color',[1 1 1]);
        plot(simul_time,simul_dX_d_x, ':r','LineWidth',lineWidth_tar); hold on;
        plot(simul_time,simul_dX_d_y, ':b','LineWidth',lineWidth_tar); hold on;
        
        plot(simul_time,simul_dX_x, 'r','LineWidth',lineWidth_cur); hold on;
        plot(simul_time,simul_dX_y, 'b','LineWidth',lineWidth_cur); hold on;
        legend({'tar_{dx}', 'tar_{dy}', 'cur_{dx}', 'cur_{dy}'}, 'location', 'best','orientation','horizontal','fontsize',15)

        axis([start_t finish_t -1.25 1.25]);
        xticks([start_t:1:finish_t])
        yticks([-1.25:0.25:1.25])
        grid on

        xlabel('time (s) ','FontSize',font_size_label)
        ylabel('Velocity (m/s) ','FontSize',font_size_label)
        title('Cartesian Space PD CTM Comtroller ','FontSize',font_size_title)
    end
end


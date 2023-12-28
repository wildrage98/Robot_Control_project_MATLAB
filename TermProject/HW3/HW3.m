clc
clear all
close all

%% HW3-1: (1-DOF) Joint Space PID CTM Controller

% 19조_2019741054_이종우_2021741071_최나은_HW3_1.m

% 목표 위치: 0 [deg] → 90 [deg]
% 목표 속도: 30 [deg/s]
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
finish_t    = 5.000;        % [초] : 종료 시간

g           = 9.8148;       % [m/s^2] : 중력 가속도

%% 로봇 파라미터

% 전역 변수
global I L g tq m;
m           = 1.0000;       % [kg]    : 링크 질량
L           = 1.0000;       % [m]     : 링크 길이
I           = (m*L^2)/3;    % [kgm^2] : 링크 관성모멘트
tq          = 0.0000;       % [Nm]    : 제어 토크

init_q      = 0.00;         % [rad]   : 초기 관절 각도
init_dq     = 0.00;         % [rad/s] : 초기 각속도
q           = init_q;       % [rad]   : 현재 관절 각도
dq          = init_dq;      % [rad/s] : 현재 각속도

% 목표 위치 파라미터
q_d         = init_q;       % [rad]     : 목표 관절 각도
dq_d        = 0;            % [rad/s]   : 목표 각속도  (정지 상태)
ddq_d       = 0;            % [rad/s^2] : 목표 각가속도 (정지 상태)
Error       = q_d - q;
Error_sum   = 0;

% 컨트롤러 게인
Wn          = 20;           % [rad/s]    : 고유 진동수 (Natural frequency)
Kp          = Wn^2;         % [Nm/rad]   : 비례 게인
Kv          = 2 * Wn;       % [Nm*s/rad] : 미분 게인 (감쇠비(zeta) : 1, 비례 게인에 임계 감쇠)
Ki          = 5000;         % [Nm/rad*s] : 적분 게인

%% 시뮬레이션 실행

if(flag_Sim == 1)
    n = 1;
    for(time = start_t:delta_t:finish_t)
        % 목표 궤적 설정
        if(time < 1)
            % 초기 1초 동안은 목표 각도를 초기 각도로, 각속도 및 각가속도를 0으로 설정
            q_d     = init_q;
            dq_d    = 0;
            ddq_d   = 0;
        else
            % 1초 이후부터는 목표 각도를 90도로 설정
            if(q_d < 90 * DR)
                % 목표 각도가 90도 미만이면 각도를 45도/1.5초로 증가시킴
                q_d = q_d + (45*DR/1.5)*delta_t;
            else
                % 목표 각도가 90도 이상이면 90도로 설정
                q_d =  90 * DR;
            end
            % 목표 각속도 및 각가속도 설정
            dq_d    = (q_d - simul_q_d(n-1))/delta_t;
            ddq_d   = (dq_d - simul_dq_d(n-1))/delta_t;
        end

        % 동역학 얻기
        G = get_Gravity(q);

        % 컨트롤러
        Error = q_d - q;
        Error_sum = Error_sum + Error;
        u = ddq_d + Kv*(dq_d - dq) + Kp*Error + Ki*Error_sum*delta_t;
        tq_ctrl = I*u + G * 0.8; % 컨트롤러가 생성한 제어 토크, 보상 토크는 중력의 80%
        
        % 로봇 모델
        % 역동역학
        tq      = tq_ctrl; % 20%의 중력 보상 오차 적용
        [t,y]   = ode45('one_link',[0, delta_t],[q; dq]); % 라그랑주 함수 수치 적분
        index   = length(y);
        q       = y(index,1);
        dq      = y(index,2);

        % 데이터 저장
        simul_time(n)   = time;     % [s]
        simul_q(n)      = q;        % [rad]
        simul_dq(n)     = dq;       % [rad/s]
        simul_q_d(n)    = q_d;      % [rad]
        simul_dq_d(n)   = dq_d;     % [rad/s]
        n               = n+1;      % 카운터 증가
    end
end

%% 시뮬레이션 결과 그래프

if(flag_Draw == 1)
    font_size_label = 20;
    font_size_title = 25;
    lineWidth_cur   = 3;
    lineWidth_tar   = 5;

    if(flag_Draw_Robot == 1)

        init_x      = L*sin(init_q);
        init_y      = -L*cos(init_q);

        FG1 = figure('Position', [500 0 700 700], 'Color',[1 1 1]);
        AX = axes('Parent',FG1);
        hold on

        p = plot([0 0],[init_x, init_y], '-ob','Linewidth',lineWidth_cur);

        axis([-1.5 1.5 -1.5 1.5]);
        grid on

        xlabel('X-axis (m)','FontSize',font_size_label)
        ylabel('Y-axis (m)','FontSize',font_size_label)
        title('1-DOF Robot', 'FontSize',font_size_title)

        n = 1;
        for(time = start_t:delta_t:finish_t)
            cmd = sprintf('Time : %2.2f' , time);
            clc
            disp(cmd)
            q = simul_q(n);
            x = L*sin(q);   y = -L*cos(q);
            Px = [0,x];     Py = [0,y];
            set(p,'XData',Px,'YData',Py)
            drawnow
            n = n+1;
        end
    end
    if(flag_Draw_Graph == 1)
        % 각도 그래프
        FG2 = figure('Position', [900 100 600 300], 'Color',[1 1 1]);
        plot(simul_time,simul_q_d * 180/pi, ':k','LineWidth',lineWidth_tar); hold on;
        plot(simul_time,simul_q * 180/pi, 'r','LineWidth',lineWidth_cur); hold on;
        legend('Desired','Current')
        axis([start_t finish_t 0 120]);
        xticks([start_t:1:finish_t])
        yticks([0:45:90])
        grid on

        xlabel('time (s) ','FontSize',font_size_label)
        ylabel('Angle (deg) ','FontSize',font_size_label)
        title('Joint Space PID CTM Comtroller ','FontSize',font_size_title)

        % 각속도 그래프
        FG3 = figure('Position', [900 400 600 300], 'Color',[1 1 1]);
        plot(simul_time,simul_dq_d * 180/pi, ':k','LineWidth',lineWidth_tar); hold on;
        plot(simul_time,simul_dq * 180/pi, 'r','LineWidth',lineWidth_cur);
        hold on
        legend('Desired','Current')

        axis([start_t finish_t 0 120]);
        xticks([start_t:1:finish_t])
        yticks([-90:45:90])
        grid on

        xlabel('time (s) ','FontSize',font_size_label)
        ylabel('Anglular Velocity (deg/s) ','FontSize',font_size_label)
        title('Joint Space PID CTM Comtroller ','FontSize',font_size_title)
    end
end

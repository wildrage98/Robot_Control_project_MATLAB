clc, clear all
close all

%%  (3-DOF) Cartesian Space PID CTM 

% 19조_2019741054_이종우_2021741071_최나은

%% Set simulation parameters

% 시뮬레이션 플래그
flag_Sim        = 1;
flag_Draw       = 1;
flag_Draw_Robot = 1;
flag_Draw_Graph = 1;

% 전역 변수
global Iz1 Iz2 Iz3 L1 L2 L3 g m1 m2 m3 r1 r2 r3 tau1 tau2 tau3

% 파라미터 설정
delta_t     = 0.005;        % [sec] : sampling time
start_t     = 0.000;        % [sec] : start time
finish_t    = 6.000;        % [sec] : finish time

g           = 9.8148;       % [m/s^2] : Gravitational Acc

% 로봇 파라미터
DR = deg2rad(1);
RD = rad2deg(1);
m1  = 0.20;     m2 = 0.20;      m3 = 0.20;       % [kg]  : Link mass
L1  = 0.50;     L2 = 0.50;      L3 = 0.50;       % [m]   : Link Length
r1  = 0.10;     r2 = 0.10;      r3 = 0.10;       % [m]   : Centor of Mass
Iz1 = 0.05;     Iz2= 0.05;      Iz3= 0.05;       % [kgm^2] : Link Inertia
tau1= 0.00;     tau2=0.00;      tau3=0.00;       % [Nm]  : Control Torque

init_q1     = -pi/2;  
init_q2     = pi/3;      
init_q3     = 2*pi/3;      % [rad] : Init Joint Angle
init_dq1    = 0.00;  
init_dq2    = 0.00;       
init_dq3    = 0.00;        % [rad/s] : Init Angular Velocity


q       = [init_q1
           init_q2
           init_q3];                                                % [rad]     : Joint Angle
dq       = [init_dq1
           init_dq2
           init_dq3];                                               % [rad/s]   : Angular Velocity


init_X = get_Kinematics3(q(1),q(2),q(3)); % [m]   : Init End-Effector Position
X      = init_X;% [rad]     : Target Joint Angle
dX     = [0;0;0];
X_d    = init_X;
dX_d   = [0;0;0];% [rad/s]   : Target Angular Velocity
ddX_d  = [0;0;0];% [rad/s^2] : Target Angular Acc
tau1 = 0.0;  tau2 = 0.0;    tau3 = 0.0;  
Error  = X_d - X;
ErrorSum = 0;

tq = [tau1 tau2 tau3];% [Nm]      :Control Torque 1, 2

% Controller Gain
Wn = 20;           % [rad/s] : Natural frequency
Kp = Wn^2;         % [Nm/rad] : Propotional Gain
Kv = 2 * Wn;       % [Nm*s/rad] : Derivative Gain
Ki = 2000.0;       % [Nm/rad*s] : Intergrator Gain

%% Simulation
if(flag_Sim == 1)
    n = 1;
    sin_t = 0;
    pre_J = 0;
    for(time = start_t:delta_t:finish_t)
         % Set Target Trajectory
        if(time < 1) % 처음 1초 동안 처음 상태 유지
            X_d     = init_X;
            dX_d    = [0;0;0];
            ddX_d   = [0;0;0];
        elseif(time < 2.0)
            X_d(1) = init_X(1);
            if(X_d(2) < init_X(2) + 0.2)
                X_d(2) = X_d(2) + (0.2/0.5)*delta_t;
            else
                X_d(2) = init_X(2) + 0.2;
            end
            dX_d    = (X_d - [simul_X_d_x(n-1) ; simul_X_d_y(n-1); 0])./delta_t;
            ddX_d    = (dX_d - [simul_dX_d_x(n-1) ; simul_dX_d_y(n-1); 0])./delta_t;
        else
            X_d     = [ 0.1 * sin((2*pi*sin_t)/2)+init_X(1);
                0.1 * cos((2*pi*sin_t)/2)+init_X(2); 0;];
            sin_t     = sin_t + (delta_t*2);
            dX_d    = (X_d - [simul_X_d_x(n-1) ; simul_X_d_y(n-1); 0])./delta_t;
            ddX_d    = (dX_d - [simul_dX_d_x(n-1) ; simul_dX_d_y(n-1); 0])./delta_t;
        end
        %동적 모델 설계에 필요한 함수 
         
        J = get_Jacobian3(q(1),q(2),q(3));
        dj = (J - pre_J)/delta_t;
        pre_J = J;
        X = get_Kinematics3(q(1),q(2),q(3));
        dX = J*dq;

        D = get_Inertia3(q(2),q(3));                      % 관성 계산 
        H = get_Coriollis3(q(2),q(3),dq(1),dq(2),dq(3));  % 전향력 계산 
        C = get_Gravity3(q(1),q(2),q(3));                 % 중력 계산 

        % 오차 계산 
        Error = X_d - X;
        ErrorSum = ErrorSum + Error;
        % 제어 
        u = ddX_d + Kv*(dX_d - dX) + Kp*(X_d - X) + Ki*ErrorSum*delta_t;
        ddq_ref = inv(J)*(u - dj*dq);

        tq_ctrl = D*ddq_ref + H + C*0.8;
        tq      = tq_ctrl;
        tau1 = tq(1); tau2 = tq(2); tau3 = tq(3);

        [t,y]   = ode45('three_link',[0, delta_t],[q(1); dq(1); q(2); dq(2); q(3); dq(3)]);
        index   = length(y);
        q(1)       = y(index,1);
        dq(1)      = y(index,2);
        q(2)       = y(index,3);
        dq(2)      = y(index,4);
        q(3)       = y(index,5);
        dq(3)      = y(index,6);
        simul_time(n)   = time;               % [sec]
        simul_q1(n)      = q(1);              % [rad]
        simul_q2(n)      = q(2);              % [rad]
        simul_q3(n)      = q(3);              % [rad]
        simul_dq1(n)     = dq(1);             % [rad/s]
        simul_dq2(n)     = dq(2);             % [rad/s]
        simul_dq3(n)     = dq(3);             % [rad/s]
        simul_X_x(n)     = X(1);              % [m]
        simul_X_y(n)     = X(2);              % [m]
        simul_dX_x(n)     = dX(1);            % [m/s]
        simul_dX_y(n)     = dX(2);            % [m/s]
        simul_X_d_x(n)     = X_d(1);          % [m]
        simul_X_d_y(n)     = X_d(2);          % [m]
        simul_dX_d_x(n)     = dX_d(1);        % [m/s]
        simul_dX_d_y(n)     = dX_d(2);        % [m/s]
        n = n+1
    end
end

%% 시뮬레이션 그래프

if(flag_Draw == 1)
    font_size_label = 20;
    font_size_title = 25;
    lineWidth_cur   = 3;
    lineWidth_tar   = 5;

    if(flag_Draw_Robot == 1)

        x1 = L1*cos(init_q1);                        % [m]   : Joint 1 X-axis Position
        y1 = L1*sin(init_q1);                        % [m]   : Joint 1 Y-axis Position
        x2 = L1*cos(init_q1+init_q2);                % [m]   : Joint 2 X-axis Position
        y2 = L1*sin(init_q1+init_q2);                % [m]   : Joint 2 Y-axis Position
        x3 = L1*cos(init_q1+init_q2+init_q3);        % [m]   : Joint 3 X-axis Position
        y3 = L1*sin(init_q1+init_q2+init_q3);        % [m]   : Joint 3 Y-axis Position



        FG1 = figure('Position', [500 0 700 700], 'Color',[1 1 1]);
        AX = axes('Parent',FG1);
        hold on

        Px1 = [0 x1];          Py1 = [0 y1];
        Px2 = [x1 x1+x2];      Py2 = [y1 y1+y2];
        Px3 = [x1+x2 x1+x2+x3];   Py3 = [y1+y2 y1+y2+y3];

        p1 = plot(Px1, Py1, '-or', 'Linewidth',3);
        p2 = plot(Px2, Py2, '-og', 'Linewidth',3);
        p3 = plot(Px3, Py3, '-ob', 'Linewidth',3);
        axis([-1.0 1.0 -1.6 0.4]);
        grid on

        xlabel('X-axis (m)','FontSize',font_size_label)
        ylabel('Y-axis (m)','FontSize',font_size_label)
        title('3-DOF Robot', 'FontSize',font_size_title)

        n = 1;
        for(time = start_t:delta_t:finish_t)
            cmd = sprintf('Time : %2.2f' , time);
            clc
            disp(cmd)
            q1 = simul_q1(n);     q2 = simul_q2(n);     q3 = simul_q3(n);
            x1 = L1*cos(q1);              % [m]   : Joint 1 X-axis Position
            y1 = L1*sin(q1);              % [m]   : Joint 1 Y-axis Position
            x2 = L1*cos(q1+q2);           % [m]   : Joint 2 X-axis Position
            y2 = L1*sin(q1+q2);           % [m]   : Joint 2 Y-axis Position
            x3 = L1*cos(q1+q2+q3);        % [m]   : Joint 3 X-axis Position
            y3 = L1*sin(q1+q2+q3);        % [m]   : Joint 3 Y-axis Position

            Px1 = [0 x1];             Py1 = [0 y1];
            Px2 = [x1 x1+x2];         Py2 = [y1 y1+y2];
            Px3 = [x1+x2 x1+x2+x3];   Py3 = [y1+y2 y1+y2+y3];
            set(p1,'XData', Px1, 'YData',Py1)
            set(p2,'XData', Px2, 'YData',Py2)
            set(p3,'XData', Px3, 'YData',Py3)
            drawnow
            n = n+1;
        end

    end
    if(flag_Draw_Graph == 1)

        FG2 = figure('Position', [900 300 600 500], 'Color',[1 1 1]);
        plot(simul_time,simul_X_d_x, ':b','LineWidth',lineWidth_tar); hold on;
        plot(simul_time,simul_X_d_y, ':r','LineWidth',lineWidth_tar); hold on;

        plot(simul_time,simul_X_x, 'b','LineWidth',lineWidth_cur); hold on;
        plot(simul_time,simul_X_y, 'r','LineWidth',lineWidth_cur); hold on;

        legend('Desired','Current')
        axis([start_t finish_t -1.25 1]);
        xticks([start_t:1:finish_t])
        yticks([-1:0.25:1])
        grid on
        legend({'tar_x','tar_y','cur_x','cur_y'},'location','best','Orientation','horizontal','FontSize',15);

        xlabel('time (s) ','FontSize',font_size_label)
        ylabel('Position (m) ','FontSize',font_size_label)
        title('Cartesian Space PID CTM Controller ','FontSize',font_size_title)



        FG3 = figure('Position', [900 300 600 500], 'Color',[1 1 1]);
        plot(simul_time,simul_dX_d_x, ':b','LineWidth',lineWidth_tar); hold on;
        plot(simul_time,simul_dX_d_y, ':r','LineWidth',lineWidth_tar); hold on;

        plot(simul_time,simul_dX_x, 'b','LineWidth',lineWidth_cur); hold on;
        plot(simul_time,simul_dX_y, 'r','LineWidth',lineWidth_cur); hold on;
        axis([start_t finish_t -1.25 1.25]);
        xticks([start_t:1:finish_t])
        yticks([-1.25:0.25:1.25])
        grid on

        legend({'tar_{dx}','tar_{dy}','cur_{dx}','cur_{dy}'},'location','best','Orientation','horizontal','FontSize',15);


        xlabel('time (s) ','FontSize',font_size_label)
        ylabel('Velocity (m/s) ','FontSize',font_size_label)
        title('Cartesian Space PD CTM Controller ','FontSize',font_size_title)
    end
end


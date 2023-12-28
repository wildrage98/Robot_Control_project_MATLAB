clc, clear all
close all

%% (3-DOF) Joint Space PID CTM 

% 19조_2019741054_이종우_2021741071_최나은

%% 시뮬레이션 설정

% 시뮬레이션 플래그
flag_Sim = 1;
flag_Draw = 1;
flag_Draw_Robot = 1;
flag_Draw_Graph = 1;

% 전역 변수
global Iz1 Iz2 Iz3 L1 L2 L3 g m1 m2 m3 r1 r2 r3 tau1 tau2 tau3

%Sim Param
delta_t     = 0.005;        % [sec] : sampling time
start_t     = 0.000;        % [sec] : start time
finish_t    = 5.000;        % [sec] : finish time

g           = 9.8148;       % [m/s^2] : Gravitational Acc

% 로봇 파라미터
DR = deg2rad(1);
RD = rad2deg(1);

m1  = 0.20;     m2 = 0.20;      m3 = 0.20;       % [kg]  : Link mass
L1  = 0.50;     L2 = 0.50;      L3 = 0.50;       % [m]   : Link Length
r1  = 0.10;     r2 = 0.10;      r3 = 0.10;       % [m]   : Centor of Mass
Iz1 = 0.05;     Iz2= 0.05;      Iz3= 0.05;       % [kgm^2] : Link Inertia
tau1= 0.00;     tau2=0.00;      tau3=0.00;       % [Nm]  : Control Torque

init_q1     = 0.00;   init_q2     = 0.00;       init_q3     = 0.00;        % [rad] : Init Joint Angle
init_dq1    = 0.00;   init_dq2    = 0.00;       init_dq3    = 0.00;        % [rad/s] : Init Angular Velocity
init_q      = [init_q1; init_q2; init_q3];                                 % [rad] : Init Joint Angles
q           = [init_q1; init_q2; init_q3];                                 % [rad] : Current Joint Angle
dq          = [init_dq1; init_dq2; init_dq3];                              % [rad/s] : Current Angular Velocity

% Target Position Param
q_d         = [init_q1; init_q2; init_q3];       % [rad] : Target Joint Angle
dq_d        = [init_dq1; init_dq2; init_dq3];                         % [rad/s] : Target Angular Velocity
ddq_d       = [0; 0; 0];                         % [rad/s^2] : Target Angular Acceleration

% 컨트롤러 게인
Wn          = 20;           % [rad/s] : Natural frequency
Kp          = Wn^2;         % [Nm/rad] : Propotional Gain
Kv          = 2 * Wn;       % [Nm*s/rad] : Derivative Gain
Ki          = 200.0;          % [] : Intergrator Gain

Error       = q_d - q;
Error_Sum   = 0;

%% Simulation

if(flag_Sim == 1)
    n = 1;
    for(time = start_t:delta_t:finish_t)
        if(time < 1)
            q_d     = init_q;
            dq_d    = [0;0;0];
            ddq_d   = [0;0;0];
        else
            if(q_d(1) < 45 * DR)
                q_d(1) = q_d(1) + (30*DR)*delta_t;
            else
                q_d(1) =  45 * DR;
            end
            if(q_d(2) < 70 * DR)
                q_d(2) = q_d(2) + (35*DR)*delta_t;
            else
                q_d(2) =  70 * DR;
            end
            if(q_d(3) < 90 * DR)
                q_d(3) = q_d(3) + (40*DR)*delta_t;
            else
                q_d(3) =  90 * DR;
            end
            dq_d(1)    = (q_d(1) - simul_q1_d(n-1))/delta_t;
            ddq_d(1)   = (dq_d(1) - simul_dq1_d(n-1))/delta_t;
            dq_d(2)    = (q_d(2) - simul_q2_d(n-1))/delta_t;
            ddq_d(2)   = (dq_d(2) - simul_dq2_d(n-1))/delta_t;
            dq_d(3)    = (q_d(3) - simul_q3_d(n-1))/delta_t;
            ddq_d(3)   = (dq_d(3) - simul_dq3_d(n-1))/delta_t;
        end

        D = get_Inertia3(q(2),q(3));
        H = get_Coriollis3(q(2),q(3),dq(1),dq(2),dq(3));
        C = get_Gravity3(q(1),q(2),q(3));

        Error = q_d - q;
        Error_Sum = Error_Sum + Error;

        u = ddq_d + Kv*(dq_d - dq) + Kp*(q_d - q) + Ki*Error_Sum*delta_t;
        tq_ctrl = D*u + H + C * 0.8;
        %         tq_ctrl = G * 0.8;

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

        simul_time(n)   = time;      % [sec]
        simul_q1(n)   = q(1);  simul_q2(n)   = q(2);  simul_q3(n)  = q(3);      % [rad]
        simul_dq1(n)   = dq(1); simul_dq2(n)   = dq(2);  simul_dq3(n)  = dq(3);     % [rad/s]
        simul_q1_d(n)    = q_d(1); simul_q2_d(n)    = q_d(2); simul_q3_d(n)    = q_d(3);    % [rad]
        simul_dq1_d(n)   = dq_d(1); simul_dq2_d(n)   = dq_d(2); simul_dq3_d(n)   = dq_d(3);  % [rad/s]
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

        q1 = init_q(1);
        q2 = init_q(2);
        q3 = init_q(3);


        x1 = L1*cos(q1);
        y1 = L1*sin(q1);
        Px1 = [0 x1];
        Py1 = [0 y1];

        x2 = L2*cos(q1+q2);
        y2 = L2*sin(q1+q2);
        Px2 = [x1 x1+x2];
        Py2 = [y1 y1+y2];

        x3 = L3*cos(q1+q2+q3);
        y3 = L3*sin(q1+q2+q3);
        Px3 = [x1+x2 x1+x2+x3];
        Py3 = [y1+y2 y1+y2+y3];

        FG1 = figure('Position', [200 300 700 700], 'Color',[1 1 1]);
        AX = axes('Parent',FG1);
        hold on

        p1 = plot(Px1, Py1, '-or', 'Linewidth',lineWidth_cur);
        p2 = plot(Px2, Py2, '-ob', 'Linewidth',lineWidth_cur);
        p3 = plot(Px3, Py3, '-og', 'Linewidth',lineWidth_cur);

        axis([-1.5 1.5 -1.5 1.5]);
        grid on

        xlabel('X-axis (m)','FontSize',font_size_label)
        ylabel('Y-axis (m)','FontSize',font_size_label)
        title('3-DOF Robot', 'FontSize',font_size_title)

        n = 1;
        for(time = start_t:delta_t:finish_t)
            cmd = sprintf('Time : %2.2f' , time);
            clc
            disp(cmd)
            q1 = simul_q1(n);
            q2 = simul_q2(n);
            q3 = simul_q3(n);

            x1 = L1*cos(q1);
            y1 = L1*sin(q1);
            Px1 = [0 x1];
            Py1 = [0 y1];

            x2 = L2*cos(q1+q2);
            y2 = L2*sin(q1+q2);
            Px2 = [x1 x1+x2];
            Py2 = [y1 y1+y2];

            x3 = L3*cos(q1+q2+q3);
            y3 = L3*sin(q1+q2+q3);
            Px3 = [x1+x2 x1+x2+x3];
            Py3 = [y1+y2 y1+y2+y3];

            set(p1,'XData', Px1, 'YData',Py1)
            set(p2,'XData', Px2, 'YData',Py2)
            set(p3,'XData', Px3, 'YData',Py3)
            drawnow
            n = n+1;
        end
    end
    if(flag_Draw_Graph == 1)

        FG2 = figure('Position', [900 400 600 500], 'Color',[1 1 1]);
        plot(simul_time,simul_q1_d * 180/pi, '--r','LineWidth',lineWidth_tar); hold on;
        plot(simul_time,simul_q1 * 180/pi, 'r','LineWidth',lineWidth_cur); hold on;
        
        plot(simul_time,simul_q2_d * 180/pi, '--g','LineWidth',lineWidth_tar); hold on;
        plot(simul_time,simul_q2 * 180/pi, 'g','LineWidth',lineWidth_cur); hold on;

        plot(simul_time,simul_q3_d * 180/pi, '--b','LineWidth',lineWidth_tar); hold on;
        plot(simul_time,simul_q3 * 180/pi, 'b','LineWidth',lineWidth_cur); hold on;

        legend({'tar_{q1}','cur_{q1}','tar_{q2}','cur_{q2}','tar_{q3}','cur_{q3}'},'location','best','Orientation','horizontal','FontSize',15);
        axis([start_t finish_t 0 120]);
        xticks([start_t:1:finish_t])
        yticks([0:45:90])
        grid on

        xlabel('time (s) ','FontSize',font_size_label)
        ylabel('Angle (deg) ','FontSize',font_size_label)
        title('Joint Space PID CTM Controller ','FontSize',font_size_title)



        FG3 = figure('Position', [900 400 600 500], 'Color',[1 1 1]);
        plot(simul_time,simul_dq1_d * 180/pi, '--r','LineWidth',lineWidth_tar); hold on;
        plot(simul_time,simul_dq1 * 180/pi, 'r','LineWidth',lineWidth_cur); hold on;

        plot(simul_time,simul_dq2_d * 180/pi, '--g','LineWidth',lineWidth_tar); hold on;
        plot(simul_time,simul_dq2 * 180/pi, 'g','LineWidth',lineWidth_cur); hold on;

        plot(simul_time,simul_dq3_d * 180/pi, '--b','LineWidth',lineWidth_tar); hold on;
        plot(simul_time,simul_dq3 * 180/pi, 'b','LineWidth',lineWidth_cur); hold on;

        legend({'tar_{dq1}','cur_{dq1}','tar_{dq2}','cur_{dq2}','tar_{dq3}','cur_{dq3}'},'location','best','Orientation','horizontal','FontSize',15);

        axis([start_t finish_t 0 120]);
        xticks([start_t:1:finish_t])
        yticks([-90:45:90])
        grid on

        xlabel('time (s) ','FontSize',font_size_label)
        ylabel('Anglular Velocity (deg/s) ','FontSize',font_size_label)
        title('Joint Space PID CTM Controller ','FontSize',font_size_title)
    end
end
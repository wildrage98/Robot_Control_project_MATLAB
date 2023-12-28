clc
clear all
close all

%% HW1: Using the Lagrangian function, derive 3 DOF robot and perform a free fall simulation. Discuss the results.

% 19조_2019741054_이종우_2021741071_최나은_HW1.m

% dydt를 얻는 부분 -> 'get_dydt.m' 스크립트 파일
% 계산량 및 실행 시간을 위해 사전에 dydt 수식을 얻어놓고, 'three_links.m' 파일에 저장
% 메인함수라 볼 수 있는 이 파일에서는 앞선 파일로부터 라그랑주 미분 방정식을 가져와 실시간 수치 적분, 시뮬레이션

%% free-fall Simulation

clear all
close all

% 전역변수, 링크 성분 초기화
global Iz1 Iz2 Iz3 L1 L2 L3 g m1 m2 m3 r1 r2 r3 tau1 tau2 tau3

L1 = 0.5; L2 = 0.5; L3 = 0.5;
r1 = 0.1; r2 = 0.1; r3 = 0.1;
m1 = 0.2; m2 = 0.2; m3 = 0.2;

Iz1 = 0.05; Iz2 = 0.05; Iz3 = 0.05;

g = 9.806; % 중력 가속도

dt = 0.02; ft = 5; % 타임 스텝 및 총 실행 시간

% 초기 관절 각도 및 각속도 설정
q1 = -pi/2; dq1 = 0;
q2 = pi/4; dq2 = 0;
q3 = pi/4; dq3 = 0;

data = []; % 실행 결과를 저장할 배열

n = 1;

% 시각화를 위한 도형, 축 생성
FG = figure('Position',[300 300 600 600],'Color',[1 1 1])
AX = axes('parent',FG);

hold on
grid on
axis([-2.0 2.0 -2.0 2.0]);

% GIF 파일 초기화
filename = 'simulation.gif';

% 초기 링크 위치 설정
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

p1 = plot(Px1,Py1,'-or','Linewidth',3);
p2 = plot(Px2,Py2,'-og','Linewidth',3);
p3 = plot(Px3,Py3,'-ob','Linewidth',3);

% 시뮬레이션 루프
for cnt=0:dt:ft
    % 관절 토크를 자유낙하 시뮬레이션을 위해 0으로 설정
    tau1 = 0.0;
    tau2 = 0.0;
    tau3 = 0.0;
    
    % 3-DOF 시스템에 대한 ode45를 사용한 수치 적분 수행
    [t, y] = ode45('three_links', [0 dt], [q1; dq1; q2; dq2; q3; dq3]);
    
    % 적분 결과에서 최종 상태 값을 추출
    index = length(y);
    q1 = y(index,1);
    dq1 = y(index,2);
    q2 = y(index,3);
    dq2 = y(index,4);
    q3 = y(index,5);
    dq3 = y(index,6);
    
    % 시각화를 위한 링크 위치 업데이트
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
    
    % 루프 카운터 증가
    n = n+1;
    
    % 시뮬레이션 실행 시간 표시
    cmd = sprintf('시간 : %2.2f',cnt);
    clc
    disp(cmd)
    
    % 부드러운 애니메이션을 위해 매 반복마다 플롯 업데이트
    if rem(n,2)==0
        set(p1,'XData',Px1,'YData',Py1)
        set(p2,'XData',Px2,'YData',Py2)
        set(p3,'XData',Px3,'YData',Py3)
        drawnow

        % 현재 프레임 캡쳐 및 GIF 파일에 저장
        frame = getframe(gcf);
        im = frame2im(frame);
        [imind,cm] = rgb2ind(im,256);
        if cnt == 0
            imwrite(imind,cm,filename,'gif', 'Loopcount',inf, 'DelayTime', 0.1);
        else
            imwrite(imind,cm,filename,'gif','WriteMode','append', 'DelayTime', 0.1);
        end
    end
end

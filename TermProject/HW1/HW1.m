clc
clear all
close all

%% HW1: Using the Lagrangian function, derive 3 DOF robot and perform a free fall simulation. Discuss the results.

% 19��_2019741054_������_2021741071_�ֳ���_HW1.m

% dydt�� ��� �κ� -> 'get_dydt.m' ��ũ��Ʈ ����
% ��귮 �� ���� �ð��� ���� ������ dydt ������ ������, 'three_links.m' ���Ͽ� ����
% �����Լ��� �� �� �ִ� �� ���Ͽ����� �ռ� ���Ϸκ��� ��׶��� �̺� �������� ������ �ǽð� ��ġ ����, �ùķ��̼�
 
%% free-fall Simulation

clear all
close all

% ��������, ��ũ ���� �ʱ�ȭ
global Iz1 Iz2 Iz3 L1 L2 L3 g m1 m2 m3 r1 r2 r3 tau1 tau2 tau3

L1 = 0.5;   L2 = 0.5;   L3 = 0.5;
r1 = 0.1;   r2 = 0.1;   r3 = 0.1;
m1 = 0.2;   m2 = 0.2;   m3 = 0.2;

Iz1 = 0.05; Iz2 = 0.05; Iz3 = 0.05;

g = 9.806; % �߷� ���ӵ�

dt = 0.02; ft = 5; % Ÿ�� ���� �� �� ���� �ð�

% �ʱ� ���� ���� �� ���ӵ� ����
q1 = -pi/2; dq1 = 0;
q2 = pi/4;  dq2 = 0;
q3 = pi/4;  dq3 = 0;

data = []; % ���� ����� ������ �迭

n = 1; 

% �ð�ȭ�� ���� ����, �� ����
FG = figure('Position',[300 300 600 600],'Color',[1 1 1])
AX = axes('parent',FG);

hold on
grid on
axis([-2.0 2.0 -2.0 2.0]);

% �ʱ� ��ũ ��ġ ����
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

% �ùķ��̼� ����
for cnt=0:dt:ft
    % ���� ��ũ�� �������� �ùķ��̼��� ���� 0���� ����
    tau1 = 0.0;
    tau2 = 0.0;
    tau3 = 0.0;
    
    % 3-DOF �ý��ۿ� ���� ode45�� ����� ��ġ ���� ����
    [t, y] = ode45('three_links', [0 dt], [q1; dq1; q2; dq2; q3; dq3]);
    
    % ���� ������� ���� ���� ���� ����
    index = length(y);
    q1  = y(index,1);
    dq1 = y(index,2);
    q2  = y(index,3);
    dq2 = y(index,4);
    q3  = y(index,5);
    dq3 = y(index,6);
    
    % �ð�ȭ�� ���� ��ũ ��ġ ������Ʈ
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
    
    % ���� ī���� ����
    n = n+1;
    
    % �ùķ��̼� ���� �ð� ǥ��
    cmd = sprintf('�ð� : %2.2f',cnt);
    clc
    disp(cmd) 
    
    % �ε巯�� �ִϸ��̼��� ���� �� �ݺ����� �÷� ������Ʈ
    if rem(n,2)==0
        set(p1,'XData',Px1,'YData',Py1)
        set(p2,'XData',Px2,'YData',Py2)
        set(p3,'XData',Px3,'YData',Py3)
        drawnow
    end
end

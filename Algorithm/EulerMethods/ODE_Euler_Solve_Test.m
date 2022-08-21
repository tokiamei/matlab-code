% 欧拉算法测试程序
% By ZFS@wust  2021
% 获取更多Matlab/Simulink原创资料和程序，清关注微信公众号：Matlab Fans

clear
clc
close all

%%  仿真步长 h=0.001 时
Hfun = @(t,x) [ x(2)*x(3); -x(1)*x(3); -0.51*x(1)*x(2)];   % 一阶微分方程导数表达式

% 参数
t = [0 12];     % 时间范围
h = 0.001;       % 时间步长
x0 = [0 1 1];   % 初始状态

% 显式欧拉法求解
[T1,X1] = ODE_ExplicitEuler( Hfun,t,h,x0 );
% 隐式欧拉法求解
[T2,X2] = ODE_ImplicitEuler( Hfun,t,h,x0 );
% 两步欧拉法求解
[T3,X3] = ODE_2OrderEuler( Hfun,t,h,x0 );
% 改进欧拉法求解
[T4,X4] = ODE_ImprovedEuler( Hfun,t,h,x0 );
% Matlab自带ode45求解
[T5,X5] = ode45( Hfun,t(1):h:t(2),x0 );

% 绘图对比
figure
subplot(311)
plot(T1,X1(:,1),T2,X2(:,1),T3,X3(:,1),T4,X4(:,1),T5,X5(:,1))
xlabel('Time(s)')
ylabel('x_1')
legend('显式欧拉法','隐式欧拉法','两步欧拉法','改进欧拉法','ode45')
subplot(312)
plot(T1,X1(:,2),T2,X2(:,2),T3,X3(:,2),T4,X4(:,2),T5,X5(:,2))
xlabel('Time(s)')
ylabel('x_2')
legend('显式欧拉法','隐式欧拉法','两步欧拉法','改进欧拉法','ode45')
subplot(313)
plot(T1,X1(:,3),T2,X2(:,3),T3,X3(:,3),T4,X4(:,3),T5,X5(:,3))
xlabel('Time(s)')
ylabel('x_3')
legend('显式欧拉法','隐式欧拉法','两步欧拉法','改进欧拉法','ode45')

%%  仿真步长 h=0.01 时

% 参数
h = 0.01;       % 时间步长

% 显式欧拉法求解
[T1,X1] = ODE_ExplicitEuler( Hfun,t,h,x0 );
% 隐式欧拉法求解
[T2,X2] = ODE_ImplicitEuler( Hfun,t,h,x0 );
% 两步欧拉法求解
[T3,X3] = ODE_2OrderEuler( Hfun,t,h,x0 );
% 改进欧拉法求解
[T4,X4] = ODE_ImprovedEuler( Hfun,t,h,x0 );
% Matlab自带ode45求解
[T5,X5] = ode45( Hfun,t(1):h:t(2),x0 );

% 绘图对比
figure
subplot(311)
plot(T1,X1(:,1),T2,X2(:,1),T3,X3(:,1),T4,X4(:,1),T5,X5(:,1))
xlabel('Time(s)')
ylabel('x_1')
legend('显式欧拉法','隐式欧拉法','两步欧拉法','改进欧拉法','ode45')
subplot(312)
plot(T1,X1(:,2),T2,X2(:,2),T3,X3(:,2),T4,X4(:,2),T5,X5(:,2))
xlabel('Time(s)')
ylabel('x_2')
legend('显式欧拉法','隐式欧拉法','两步欧拉法','改进欧拉法','ode45')
subplot(313)
plot(T1,X1(:,3),T2,X2(:,3),T3,X3(:,3),T4,X4(:,3),T5,X5(:,3))
xlabel('Time(s)')
ylabel('x_3')
legend('显式欧拉法','隐式欧拉法','两步欧拉法','改进欧拉法','ode45')


%%  仿真步长 h=0.1 时
% 参数
h = 0.1;       % 时间步长

% 显式欧拉法求解
[T1,X1] = ODE_ExplicitEuler( Hfun,t,h,x0 );
% 隐式欧拉法求解
[T2,X2] = ODE_ImplicitEuler( Hfun,t,h,x0 );
% 两步欧拉法求解
[T3,X3] = ODE_2OrderEuler( Hfun,t,h,x0 );
% 改进欧拉法求解
[T4,X4] = ODE_ImprovedEuler( Hfun,t,h,x0 );
% Matlab自带ode45求解
[T5,X5] = ode45( Hfun,t(1):h:t(2),x0 );


% 绘图对比
figure
subplot(311)
plot(T1,X1(:,1),T2,X2(:,1),T3,X3(:,1),T4,X4(:,1),T5,X5(:,1))
xlabel('Time(s)')
ylabel('x_1')
legend('显式欧拉法','隐式欧拉法','两步欧拉法','改进欧拉法','ode45')
subplot(312)
plot(T1,X1(:,2),T2,X2(:,2),T3,X3(:,2),T4,X4(:,2),T5,X5(:,2))
xlabel('Time(s)')
ylabel('x_2')
legend('显式欧拉法','隐式欧拉法','两步欧拉法','改进欧拉法','ode45')
subplot(313)
plot(T1,X1(:,3),T2,X2(:,3),T3,X3(:,3),T4,X4(:,3),T5,X5(:,3))
xlabel('Time(s)')
ylabel('x_3')
legend('显式欧拉法','隐式欧拉法','两步欧拉法','改进欧拉法','ode45')



%%  仿真步长 h=0.5 时
% 参数
h = 0.5;       % 时间步长

% 显式欧拉法求解
[T1,X1] = ODE_ExplicitEuler( Hfun,t,h,x0 );
% 隐式欧拉法求解
[T2,X2] = ODE_ImplicitEuler( Hfun,t,h,x0 );
% 两步欧拉法求解
[T3,X3] = ODE_2OrderEuler( Hfun,t,h,x0 );
% 改进欧拉法求解
[T4,X4] = ODE_ImprovedEuler( Hfun,t,h,x0 );
% Matlab自带ode45求解
[T5,X5] = ode45( Hfun,t(1):h:t(2),x0 );


% 绘图对比
figure
subplot(311)
plot(T3,X3(:,1),T4,X4(:,1),T5,X5(:,1))
xlabel('Time(s)')
ylabel('x_1')
legend('两步欧拉法','改进欧拉法','ode45')
subplot(312)
plot(T3,X3(:,2),T4,X4(:,2),T5,X5(:,2))
xlabel('Time(s)')
ylabel('x_2')
legend('两步欧拉法','改进欧拉法','ode45')
subplot(313)
plot(T3,X3(:,3),T4,X4(:,3),T5,X5(:,3))
xlabel('Time(s)')
ylabel('x_3')
legend('两步欧拉法','改进欧拉法','ode45')


%%  仿真步长 h=0.6 时
% 参数
h = 0.6;       % 时间步长

% 显式欧拉法求解
[T1,X1] = ODE_ExplicitEuler( Hfun,t,h,x0 );
% 隐式欧拉法求解
[T2,X2] = ODE_ImplicitEuler( Hfun,t,h,x0 );
% 两步欧拉法求解
[T3,X3] = ODE_2OrderEuler( Hfun,t,h,x0 );
% 改进欧拉法求解
[T4,X4] = ODE_ImprovedEuler( Hfun,t,h,x0 );
% Matlab自带ode45求解
[T5,X5] = ode45( Hfun,t(1):h:t(2),x0 );


% 绘图对比
figure
subplot(311)
plot(T4,X4(:,1),T5,X5(:,1))
xlabel('Time(s)')
ylabel('x_1')
legend('改进欧拉法','ode45')
subplot(312)
plot(T4,X4(:,2),T5,X5(:,2))
xlabel('Time(s)')
ylabel('x_2')
legend('改进欧拉法','ode45')
subplot(313)
plot(T4,X4(:,3),T5,X5(:,3))
xlabel('Time(s)')
ylabel('x_3')
legend('改进欧拉法','ode45')




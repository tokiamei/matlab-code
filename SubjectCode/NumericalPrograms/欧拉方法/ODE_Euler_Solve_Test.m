% ŷ���㷨���Գ���
% By ZFS@wust  2021
% ��ȡ����Matlab/Simulinkԭ�����Ϻͳ������ע΢�Ź��ںţ�Matlab Fans

clear
clc
close all

%%  ���沽�� h=0.001 ʱ
Hfun = @(t,x) [ x(2)*x(3); -x(1)*x(3); -0.51*x(1)*x(2)];   % һ��΢�ַ��̵������ʽ

% ����
t = [0 12];     % ʱ�䷶Χ
h = 0.001;       % ʱ�䲽��
x0 = [0 1 1];   % ��ʼ״̬

% ��ʽŷ�������
[T1,X1] = ODE_ExplicitEuler( Hfun,t,h,x0 );
% ��ʽŷ�������
[T2,X2] = ODE_ImplicitEuler( Hfun,t,h,x0 );
% ����ŷ�������
[T3,X3] = ODE_2OrderEuler( Hfun,t,h,x0 );
% �Ľ�ŷ�������
[T4,X4] = ODE_ImprovedEuler( Hfun,t,h,x0 );
% Matlab�Դ�ode45���
[T5,X5] = ode45( Hfun,t(1):h:t(2),x0 );

% ��ͼ�Ա�
figure
subplot(311)
plot(T1,X1(:,1),T2,X2(:,1),T3,X3(:,1),T4,X4(:,1),T5,X5(:,1))
xlabel('Time(s)')
ylabel('x_1')
legend('��ʽŷ����','��ʽŷ����','����ŷ����','�Ľ�ŷ����','ode45')
subplot(312)
plot(T1,X1(:,2),T2,X2(:,2),T3,X3(:,2),T4,X4(:,2),T5,X5(:,2))
xlabel('Time(s)')
ylabel('x_2')
legend('��ʽŷ����','��ʽŷ����','����ŷ����','�Ľ�ŷ����','ode45')
subplot(313)
plot(T1,X1(:,3),T2,X2(:,3),T3,X3(:,3),T4,X4(:,3),T5,X5(:,3))
xlabel('Time(s)')
ylabel('x_3')
legend('��ʽŷ����','��ʽŷ����','����ŷ����','�Ľ�ŷ����','ode45')

%%  ���沽�� h=0.01 ʱ

% ����
h = 0.01;       % ʱ�䲽��

% ��ʽŷ�������
[T1,X1] = ODE_ExplicitEuler( Hfun,t,h,x0 );
% ��ʽŷ�������
[T2,X2] = ODE_ImplicitEuler( Hfun,t,h,x0 );
% ����ŷ�������
[T3,X3] = ODE_2OrderEuler( Hfun,t,h,x0 );
% �Ľ�ŷ�������
[T4,X4] = ODE_ImprovedEuler( Hfun,t,h,x0 );
% Matlab�Դ�ode45���
[T5,X5] = ode45( Hfun,t(1):h:t(2),x0 );

% ��ͼ�Ա�
figure
subplot(311)
plot(T1,X1(:,1),T2,X2(:,1),T3,X3(:,1),T4,X4(:,1),T5,X5(:,1))
xlabel('Time(s)')
ylabel('x_1')
legend('��ʽŷ����','��ʽŷ����','����ŷ����','�Ľ�ŷ����','ode45')
subplot(312)
plot(T1,X1(:,2),T2,X2(:,2),T3,X3(:,2),T4,X4(:,2),T5,X5(:,2))
xlabel('Time(s)')
ylabel('x_2')
legend('��ʽŷ����','��ʽŷ����','����ŷ����','�Ľ�ŷ����','ode45')
subplot(313)
plot(T1,X1(:,3),T2,X2(:,3),T3,X3(:,3),T4,X4(:,3),T5,X5(:,3))
xlabel('Time(s)')
ylabel('x_3')
legend('��ʽŷ����','��ʽŷ����','����ŷ����','�Ľ�ŷ����','ode45')


%%  ���沽�� h=0.1 ʱ
% ����
h = 0.1;       % ʱ�䲽��

% ��ʽŷ�������
[T1,X1] = ODE_ExplicitEuler( Hfun,t,h,x0 );
% ��ʽŷ�������
[T2,X2] = ODE_ImplicitEuler( Hfun,t,h,x0 );
% ����ŷ�������
[T3,X3] = ODE_2OrderEuler( Hfun,t,h,x0 );
% �Ľ�ŷ�������
[T4,X4] = ODE_ImprovedEuler( Hfun,t,h,x0 );
% Matlab�Դ�ode45���
[T5,X5] = ode45( Hfun,t(1):h:t(2),x0 );


% ��ͼ�Ա�
figure
subplot(311)
plot(T1,X1(:,1),T2,X2(:,1),T3,X3(:,1),T4,X4(:,1),T5,X5(:,1))
xlabel('Time(s)')
ylabel('x_1')
legend('��ʽŷ����','��ʽŷ����','����ŷ����','�Ľ�ŷ����','ode45')
subplot(312)
plot(T1,X1(:,2),T2,X2(:,2),T3,X3(:,2),T4,X4(:,2),T5,X5(:,2))
xlabel('Time(s)')
ylabel('x_2')
legend('��ʽŷ����','��ʽŷ����','����ŷ����','�Ľ�ŷ����','ode45')
subplot(313)
plot(T1,X1(:,3),T2,X2(:,3),T3,X3(:,3),T4,X4(:,3),T5,X5(:,3))
xlabel('Time(s)')
ylabel('x_3')
legend('��ʽŷ����','��ʽŷ����','����ŷ����','�Ľ�ŷ����','ode45')



%%  ���沽�� h=0.5 ʱ
% ����
h = 0.5;       % ʱ�䲽��

% ��ʽŷ�������
[T1,X1] = ODE_ExplicitEuler( Hfun,t,h,x0 );
% ��ʽŷ�������
[T2,X2] = ODE_ImplicitEuler( Hfun,t,h,x0 );
% ����ŷ�������
[T3,X3] = ODE_2OrderEuler( Hfun,t,h,x0 );
% �Ľ�ŷ�������
[T4,X4] = ODE_ImprovedEuler( Hfun,t,h,x0 );
% Matlab�Դ�ode45���
[T5,X5] = ode45( Hfun,t(1):h:t(2),x0 );


% ��ͼ�Ա�
figure
subplot(311)
plot(T3,X3(:,1),T4,X4(:,1),T5,X5(:,1))
xlabel('Time(s)')
ylabel('x_1')
legend('����ŷ����','�Ľ�ŷ����','ode45')
subplot(312)
plot(T3,X3(:,2),T4,X4(:,2),T5,X5(:,2))
xlabel('Time(s)')
ylabel('x_2')
legend('����ŷ����','�Ľ�ŷ����','ode45')
subplot(313)
plot(T3,X3(:,3),T4,X4(:,3),T5,X5(:,3))
xlabel('Time(s)')
ylabel('x_3')
legend('����ŷ����','�Ľ�ŷ����','ode45')


%%  ���沽�� h=0.6 ʱ
% ����
h = 0.6;       % ʱ�䲽��

% ��ʽŷ�������
[T1,X1] = ODE_ExplicitEuler( Hfun,t,h,x0 );
% ��ʽŷ�������
[T2,X2] = ODE_ImplicitEuler( Hfun,t,h,x0 );
% ����ŷ�������
[T3,X3] = ODE_2OrderEuler( Hfun,t,h,x0 );
% �Ľ�ŷ�������
[T4,X4] = ODE_ImprovedEuler( Hfun,t,h,x0 );
% Matlab�Դ�ode45���
[T5,X5] = ode45( Hfun,t(1):h:t(2),x0 );


% ��ͼ�Ա�
figure
subplot(311)
plot(T4,X4(:,1),T5,X5(:,1))
xlabel('Time(s)')
ylabel('x_1')
legend('�Ľ�ŷ����','ode45')
subplot(312)
plot(T4,X4(:,2),T5,X5(:,2))
xlabel('Time(s)')
ylabel('x_2')
legend('�Ľ�ŷ����','ode45')
subplot(313)
plot(T4,X4(:,3),T5,X5(:,3))
xlabel('Time(s)')
ylabel('x_3')
legend('�Ľ�ŷ����','ode45')




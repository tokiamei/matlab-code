% ================== 初始化 =================================
clc;close all; clf;clear;
% ================== cputime ？为什么要写这个 ？=======================
cputime=0;
% ================== 开启定时器 =============================
tic;
% ================== 基本参数设定 ===========================
m = 1.443*1e-25; % 这是谁的质量 ？？？？
h_b = (6.626/(2*pi))*1e-34; % 约化普朗克常数
w = 2*pi*1e3; % 简谐势阱频率
dime = sqrt(h_b/(m*w)); % 坐标变量无量纲化系数

% ================= 坐标离散化 ==============================
x = -10:0.02:10-0.02;
y = -10:0.02:10-0.02;
[X,Y] = meshgrid(x,y) ;% 将 x 和 y 在二维平面内离散
delta_tau=0.01;%虚时间间隔
delta_x=0.02;%坐标间隔
delta_y=0.02;
N=1000;%x 或 y 方向离散点个数
kx=((2*pi)/(N*delta_x))*(-0.5*N:1:0.5*N-1);%kx 和 ky 分别是快速傅里叶变换后，动量空间内的坐标
ky=((2*pi)/(N*delta_y))*(-0.5*N:1:0.5*N-1);
kx_=fftshift(kx);
ky_=fftshift(ky);
[KX,KY]=meshgrid(kx_,ky_);%将 omega_x 和 omega_y 在二维空间内离散，用于快速傅里叶变换
K=sqrt(KX.^2+KY.^2);%后文用的主要是这个
%哈密顿算符中参数的设置。
delta=0;%单光子失谐
n1=1;
n2=-1;
m1=1;
m2=1;
omega0=1;
a_osc=0.34*1e-6;
beamwidth=5*a_osc;
I0=1;
I1=I0*((sqrt((dime*X).^2+(dime*Y).^2)/beamwidth).^(2*abs(n1))).*(exp(-((dime*X).^2+(dime*Y).^2)/beamwidth^2)).^2;%结合 m 和 n 的取值，拉盖尔多项式为 1，因此略去了拉盖尔多项式
I2=I0*((sqrt((dime*X).^2+(dime*Y).^2)/beamwidth).^(2*abs(n2))).*(exp(-((dime*X).^2+(dime*Y).^2)/beamwidth^2)).^2;
g=1;
g_up_down=1;
%设置初始波函数，并绘制出图像。
a1=2;
a2=0;
b1=0.01;
b2=0.01;
c1=0.01;
c2=0.01;
psi_former_up=a1*exp(-(b1*(X).^2+c1*(Y).^2));%自旋向上初始的波函数
psi_former_down=a2*exp(-(b2*(X).^2+c2*(Y).^2));%自旋向下的初始波函数
A=1/(sqrt(sum(sum(abs(psi_former_up).^2*delta_x*delta_y+abs(psi_former_down).^2*delta_x*delta_y))));
psi_former_up=A*psi_former_up;%自旋向上归一化的初始波函数
psi_former_down=A*psi_former_down;% 自 旋 向 下 归 一 化 的 初 始 波 函 数
% A2*psi_former_down;
psi_former_up(1,:)=zeros(1,1000);%边界条件
psi_former_up(1000,:)=zeros(1,1000);
psi_former_up(:,1)=zeros(1000,1);
psi_former_up(:,1000)=zeros(1000,1);
psi_former_down(1,:)=zeros(1,1000);%边界条件
psi_former_down(1000,:)=zeros(1,1000);
psi_former_down(:,1)=zeros(1000,1);
psi_former_down(:,1000)=zeros(1000,1);
figure%绘制初始自旋向上波函数的图像
subplot(2,2,1)
mesh(X,Y,abs(psi_former_up).^2);
xlabel('x'); ylabel('y'); zlabel('psi_former_up');
subplot(2,2,2)%绘制初始自旋向下波函数的图像

mesh(X,Y,abs(psi_former_down).^2);
xlabel('x'); ylabel('y'); zlabel('psi_former_down');
%虚时演化
gap=1;%虚时演化一个时间步长前后，用于存储波函数最大差值的模
n=1;
judge_gap=1e-8;%检验是否为基态波函数的判定值，用每次实际得到的 gap
while gap>judge_gap
spectrum_up=fft2(psi_former_up);%半个时间步长内，动能算符，自旋
spectrum_up=exp(-delta_tau*(1/4).*K.^2).*spectrum_up;
psi_1_up=ifft2(spectrum_up);
spectrum_down=fft2(psi_former_down);%半个时间步长内，动能算符，
spectrum_down=exp(-delta_tau*(1/4).*K.^2).*spectrum_down;
psi_1_down=ifft2(spectrum_down);
%下面的 psi_2_up_和 psi_2_down_是指未经过坐标变换的拉曼耦合项

% psi_2_up 和 psi_2_down。
%求解对角化的拉曼耦合项对应的薛定谔方程后，得到的态记做psi_3_up_和 psi_3_down_，进行坐标逆变换（左乘变换矩阵）后的态是psi_3_up 和 psi_3_down。
psi_2_up_=exp(-delta_tau*((1/2)*(g*(abs(psi_1_up).^2)+g_up_down*(abs(psi_1_down).^2))+(1/2)*((X).^2+(Y).^2))+delta/2).*psi_1_up;%算符劈裂后，N、V 算符和双光子失谐，波函数的求解
psi_2_down_=exp(-delta_tau*((1/2)*(g_up_down*(abs(psi_1_up).^2)+g*(abs(psi_1_down).^2))+(1/2)*((X).^2+(Y).^2))-delta/2).*psi_1_down;
psi_2_up=((sqrt(1+exp(-2i*(n1-n2)*angle(X+1i*Y))))./(2*exp(-1i*(n1-n2)*angle(X+1i*Y)))).*psi_2_up_+((sqrt(1+exp(-2i*(n1-n2)*angle(X+1i*Y))))./2).*psi_2_down_;%(坐标变换)乘逆的变换矩阵
psi_2_down=(-(sqrt(1+exp(-2i*(n1-n2)*angle(X+1i*Y))))./(2*exp(-1i*(n1-n2)*angle(X+1i*Y)))).*psi_2_up_+((sqrt(1+exp(-2i*(n1-n2)*angle(X+1i*Y))))./2).*psi_2_down_;

psi_3_up_=exp(-delta_tau*(omega0*sqrt(I1.*I2))).*psi_2_up;% 算符劈裂后，求解拉曼耦合项的求解
psi_3_down_=exp(-delta_tau*(-omega0*sqrt(I1.*I2))).*psi_2_down;
psi_3_up=((exp(-1i*(n1-n2)*angle(X+1i*Y)))./((sqrt(1+exp(-2i*(n1-n2)*angle(X+1i*Y)))))).*psi_3_up_+(-(exp(-1i*(n1-n2)*angle(X+1i*Y)))./((sqrt(1+exp(-2i*(n1-n2)*angle(X+1i*Y)))))).*psi_3_down_;%(坐标逆变换)乘变换矩阵
psi_3_down=(1./((sqrt(1+exp(-2i*(n1-n2)*angle(X+1i*Y)))))).*psi_3_up_+(1./((sqrt(1+exp(-2i*(n1-n2)*angle(X+1i*Y)))))).*psi_3_down_;
psi_4_up=exp(-delta_tau*((1/2)*((g*(abs(psi_3_up).^2)+g_up_down*(abs(psi_3_down).^2))/(1))+(1/2)*((X).^2+(Y).^2))+delta/2).*psi_3_up;%算符劈裂后，N、V 算符和双光子失谐，波函数的求解
psi_4_down=exp(-delta_tau*((1/2)*((g_up_down*(abs(psi_3_up).^2)+g*(abs(psi_3_down).^2))/(1))+(1/2)*((X).^2+(Y).^2))-delta/2).*psi_3_down;
spectrum_up=fft2(psi_4_up);%动能算符，自旋向上波函数的求解
spectrum_up=exp(-delta_tau*(1/4).*K.^2).*spectrum_up;
psi_5_up=ifft2(spectrum_up);
spectrum_down=fft2(psi_4_down);%动能算符，自旋向下波函数的求解
spectrum_down=exp(-delta_tau*(1/4).*K.^2).*spectrum_down;
psi_5_down=ifft2(spectrum_down);
coefficiency=(abs(psi_5_up).^2)*delta_x*delta_y+(abs(psi_5_down).^2)*delta_x*delta_y;
nor_coeffi=1/sqrt(sum(sum(coefficiency)));%归一化系数
psi_later_up=nor_coeffi*psi_5_up;%经过一个时间步长演化后的归一化

psi_later_down=nor_coeffi*psi_5_down;
Max_up=max(max(abs(psi_later_up-psi_former_up)));
Max_down=max(max(abs(psi_later_down-psi_former_down)));
gap=max(Max_up,Max_down);%波函数差值的模的最大值
psi_former_up=psi_later_up;

psi_former_down=psi_later_down;
n=n+1;
end
subplot(2,2,3)
mesh(X,Y,abs(psi_later_up).^2);
xlabel('x'); ylabel('y'); zlabel('psiup');
subplot(2,2,4)
mesh(X,Y,abs(psi_later_down).^2);
xlabel('x'); ylabel('y'); zlabel('psidown');
toc;
cputime=toc;
disp('CPU time:'), disp(cputime);
%================ 二维SOC凝聚体相图 =================%
clear;clc;clf;

%================ 参数区 ============================%
k0 = 10; % 激光波矢
delta = 0; % 失谐
omega = 1; % 拉曼耦合强度
w = 1; % 控制外势大小
% g0 = 9; % 控制相互作用大小
g11 = 10; g12 = 9;  g22 = 9;

%================ 坐标 与 时间 离散化 ===============%
len = 5; 
N = 100;
x = linspace(-len,len,N);
y = linspace(-len,len,N);
[xx,yy] = meshgrid(x,y);
dx = x(2)-x(1);
dt = 0.001; T = 10;
t = 0:dt:T;

%================ 势场 =============================%
% w1 = 0.1; w2 = 0.1;
% V = 1/2*(w1^2*xx.^2+w2^2*yy.^2);
V = 400*ones(N);
V(10:90,10:90) = 0;

%================ 设定初始波函数 ===================%
psiup = rand(N,N)+1i*rand(N,N);
psidn = rand(N,N)+1i*rand(N,N);
% psidn = 0.1*xx+0.1*yy;
psiup(1,:)= 0; psiup(end,:) = 0;
psiup(:,1) = 0; psiup(:,end) = 0;
psidn(1,:)= 0; psidn(end,:) = 0;
psidn(:,1) = 0; psidn(:,end) = 0;

%=============== 归一化 ============================%
norm = dx^2*sqrt(sum(abs(psiup).^2+abs(psidn).^2,'all'));
psiup = psiup/norm;
psidn = psidn/norm;

%================ 导数矩阵 =========================%
A = (-2*eye(N)+diag(ones(N-1,1),1)+diag(ones(N-1,1),-1))/dx^2;
B = (diag(ones(N-1,1),1)-diag(ones(N-1,1),-1))/2/dx;
HA = 1/2*A-1i*k0*B-delta/2*eye(N); % 简化一下迭代矩阵
HB = 1/2*A+1i*k0*B+delta/2*eye(N);

%================ 开始迭代 =========================%
for i = 1:length(t)-1
    psiup = psiup + dt*(HA*psiup+1/2*psiup*A-V.*psiup-omega/2*psidn-g11*abs(psiup).^2.*psiup-g12*abs(psiup).^2.*psidn);
    psidn = psidn + dt*(HB*psidn+1/2*psidn*A-V.*psidn-omega/2*psiup-g22*abs(psidn).^2.*psidn-g12*abs(psidn).^2.*psiup);
end

%=============== 归一化 ============================%
norm = dx^2*sqrt(sum(abs(psiup).^2+abs(psidn).^2,'all'));
psiup = psiup/norm; % 归一化
psidn = psidn/norm;

% ============== 绘图 ==============================%
figure(1);
subplot(2,2,1); surf(xx,yy,abs(psiup).^2); view(2); shading interp; colormap summer;colorbar;
subplot(2,2,2); surf(xx,yy,angle(psiup)); view(2); shading interp; colormap summer;colorbar;

subplot(2,2,3); surf(xx,yy,abs(psidn).^2); view(2); shading interp; colormap winter; colorbar;
subplot(2,2,4); surf(xx,yy,angle(psidn)); view(2); shading interp; colormap winter; colorbar;


















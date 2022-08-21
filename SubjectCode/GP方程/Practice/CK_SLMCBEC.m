%======================== 1D SOC BEC基态波函数 ===========
%======================== 向前欧拉方法 =============
clear;clc;
%======================== 参数区 ===========================
w = 1; % trapping frequency
k0 = 4; %激光波数
o = 0.015;
delta = 0; % 失谐
b11 = o; b12 = 0.8*o; b21 = 0.8*o; b22 = o;
% b11 = 0; b12 = 0; b21 = 0; b22 = 0;
omiga = 10; % 拉曼耦合强度

%======================== 时间和坐标的离散化 =============
L = 5;
a = -L;b = L;
dx = 0.1;
x = a:dx:b;
T = 80;
dt = 0.001;
t = 0:dt:T;

%======================= 初态 ============================
phi1 = rand(length(x), 1)+1i*rand(length(x), 1);
phi2 = rand(length(x), 1)+1i*rand(length(x), 1);

norm = sqrt(sum(abs(phi1.^2)+abs(phi2.^2)));
phi1 = phi1/norm;
phi2 = phi2/norm;

%======================= 算符矩阵 ========================
A = (-2*diag(ones(length(x),1)) + diag(ones(length(x)-1,1),1) + diag(ones(length(x)-1,1),-1))/dx^2;
%一阶导数算符矩阵
B = (diag(ones(length(x)-1,1),1) - diag(ones(length(x)-1,1),-1))/2/dx;

%======================= 势场 ============================
% v = 1/2*w*x.^2;
v = w*100*ones(1,length(x)); % box potential
v(10:90) = 0;
V = diag(v);

%======================= 迭代 =============================
for i = 1:length(t)-1
   phi1 = phi1 + dt*(1/2*A*phi1-V*phi1-1i*k0*B*phi1-delta/2*phi2-b11*abs(phi1.^2).*phi1+...
                      -b12*abs(phi2.^2).*phi1-omiga/2*phi2);

   phi2 = phi2 + dt*(1/2*A*phi2-V*phi2+1i*k0*B*phi2+delta/2*phi2-b21*abs(phi1.^2).*phi2+...
                      -b22*abs(phi2.^2).*phi2-omiga/2*phi1);
end

%===================== 归一化 ============================
norm = sqrt(0.01*sum(abs(phi1).^2+abs(phi2).^2));
phi1 = phi1/norm;
phi2 = phi2/norm;

%===================== 绘图 ==============================
f1 = figure(1);
plot(x,abs(phi1).^2,'--r',x,abs(phi2).^2,'g');
axis([a b 0 1]);
title('Zero mode phase');
xlabel('x'); ylabel('\rho');
legend('\rho_{\uparrow}','\rho_{\downarrow}');

% subplot(2,2,2); plot(x,angle(phi1(:,end)));
% subplot(2,2,3); plot(x,abs(phi2(:,end)).^2,'g');
% subplot(2,2,4); plot(x,angle(phi2(:,end)));

%===================== 保存图片 ==========================
% path = 'd:\matlab code\figure\SOC凝聚体\';
% str = strcat(path,'1D条纹相','.jpg');
% saveas(f1,str);











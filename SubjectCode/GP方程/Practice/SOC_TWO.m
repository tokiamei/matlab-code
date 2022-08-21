% ======================== 1D 两种类原子的自旋交换相互作用 SOC BEC 基态波函数 ===========
% ======================== 向前欧拉方法 =============
clear;clc;close;
% ======================== 参数区 ===========================
w = 1;
k0 = 4; % 激光波数
omegaA = 150; % A 原子的拉曼耦合强度
omegaB = 0; % B 原子拉曼耦合强度为 0，即 B 原子不产生自旋轨道耦合效应
g11 = 0.01; g22 = g11; g12 = 0.8*g11; % A 原子之间的相互作用强度
b11 = 0.6*g11; b22 = b11; b12 = 0.8*b11; % B 原子之间的相互作用强度
it = 0.006; % 两种原子种间相互作用
beta = 0.02; % 自旋交换相互作用强度
NA = 1000; % A 原子原子数
NB = 40; % B 原子原子数
% ======================== 时间和坐标的离散化 =============
N = 100;
L = 5;
a = -L;b = L;
x = linspace(a,b,N);
dx = x(2)-x(1);
T = 40;
dt = 0.001;
t = 0:dt:T;

% ======================= 初态 ============================
psiAup = rand(N,1)+1i*rand(N,1);
psiAdn = rand(N,1)+1i*rand(N,1);
psiBup = rand(N,1)+1i*rand(N,1);
psiBdn = rand(N,1)+1i*rand(N,1);

% ======================= 归一化 ==========================
normA = sqrt(sum(abs(psiAup.^2)+abs(psiAdn.^2)));
normB = sqrt(sum(abs(psiBup.^2)+abs(psiBdn.^2)));
psiAup = psiAup/normA;
psiAdn = psiAdn/normA;
psiBup = psiBup/normB;
psiBdn = psiBdn/normB;

% ======================= 算符矩阵 ========================
A = (-2*diag(ones(N,1)) + diag(ones(N-1,1),1) + diag(ones(N-1,1),-1))/dx^2;
B = (diag(ones(N-1,1),1) - diag(ones(N-1,1),-1))/2/dx;

% ======================= 势场 ============================

% v = 1/2*w*x.^2; % 简谐势场
v = w*100*ones(1,length(x)); % box potential
v(10:90) = 0;
V = diag(v);

% ======================= 迭代 =============================
for i = 1:length(t)-1
    psiAup = psiAup+dt*(1/2*A*psiAup-1i*k0*B*psiAup-V*psiAup-omegaA/2*psiAdn-NA*g11*abs(psiAup.^2).*psiAup+...
    -NA*g12*abs(psiAdn.^2).*psiAup-NB*it*psiAup.*(abs(psiBup).^2+abs(psiBdn).^2)-NB*beta*conj(psiBdn).*psiBup.*psiAdn);

    psiAdn = psiAdn+dt*(1/2*A*psiAdn+1i*k0*B*psiAdn-V*psiAdn-omegaA/2*psiAup-NA*g11*abs(psiAdn.^2).*psiAdn+...
    -NA*g12*abs(psiAup.^2).*psiAdn-NB*it*psiAdn.*(abs(psiBup).^2+abs(psiBdn).^2)-NB*beta*conj(psiBup).*psiBdn.*psiAup);
    
    
    psiBup = psiBup+dt*(1/2*A*psiBup-1i*k0*B*psiBup-V*psiBup-NB*b11*abs(psiBup.^2).*psiBup+...
    -NB*b12*abs(psiBdn.^2).*psiBup-NA*it*psiBup.*(abs(psiAup).^2+abs(psiAdn).^2)-NA*beta*conj(psiAdn).*psiBdn.*psiAup);

    psiBdn = psiBdn+dt*(1/2*A*psiBdn+1i*k0*B*psiBdn-V*psiBdn-NB*b11*abs(psiBdn.^2).*psiBdn+...
    -NB*b12*abs(psiBup.^2).*psiBdn-NA*it*psiBdn.*(abs(psiAup).^2+abs(psiAdn).^2)-NA*beta*conj(psiAup).*psiBup.*psiAdn);
end

% ===================== 归一化 ============================
normA = sqrt(dx^2*sum(abs(psiAup.^2)+abs(psiAdn.^2)));
normB = sqrt(dx^2*sum(abs(psiBup.^2)+abs(psiBdn.^2)));
psiAup = psiAup/normA;
psiAdn = psiAdn/normA;
psiBup = psiBup/normB;
psiBdn = psiBdn/normB;

% ===================== 绘图 ==============================
figure;
plot(x,abs(psiAup).^2,'--r',x,abs(psiAdn).^2,'r');
title('A 原子概率密度');
legend('\rho_{A,\uparrow}','\rho_{A,\downarrow}');
figure;
plot(x,abs(psiBup).^2,'--b',x,abs(psiBdn).^2,'b'); 
title('B 原子概率密度');
legend('\rho_{B,\uparrow}','\rho_{B,\downarrow}');



%===================== 保存图片 ==========================
% path = 'd:\matlab code\figure\SOC凝聚体\';
% str = strcat(path,'1D条纹相','.jpg');
% saveas(f1,str);

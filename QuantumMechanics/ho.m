clc
clear

X = 10;     % 长度 
N = 200;    % 格点数
dx = 2*X/N;   % 空间步长
x = linspace(-X,X,N);
n = 2;      % 次数n
V = 1/2*x'.^n;  % 势场

A = spdiags(V,0,N,N);  % 提取V的对角线并生成矩阵
% 计算哈密度算符
H = zeros(N,3);
H(1:N,1) = -0.5/dx^2;
H(1:N,2) = 1/dx^2;
H(1:N,3) = -0.5/dx^2;
B = spdiags(H,-1:1,N,N);
C = A+B;
E = 9;    % 能级个数

% 求特征值
[Vector, Value] = eigs(C,E,0);  % 求指定的几个特征值

% 画图
for i = 1:E
    psi = Vector(:,i);
    subplot(3,3,i)
    plot(x,psi.^2/sum(psi.^2*dx))
    ylim([0 0.8])
    xlabel('x');ylabel('P');
end
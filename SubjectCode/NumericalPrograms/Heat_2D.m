clear;
a = 1;
dx = 0.05;
Lx = 1;
x = 0:dx:Lx;
dy = 0.05;
Ly = 2;
y = 0:dy:Ly;
A1 = -2*eye(length(x)) + diag(ones(length(x)-1,1),1) + diag(ones(length(x)-1,1),-1);
A1(1,end) = 1; A1(end,1) = 1; %周期边界条件
A2 = -2*eye(length(y)) + diag(ones(length(y)-1,1),1) + diag(ones(length(y)-1,1),-1);
A2(1,end) = 1; A2(end,1) = 1; %周期边界条件
dt = 0.0001;
t = 0:dt:1;
[Y,X] = meshgrid(y,x); %x和y调换下位置就可以实现非方矩阵运算
U0 = exp(-10*((X-1/2).^2+(Y-1/2).^2));
U = U0;
% F = zeros(length(x),length(y));

for n = 1:length(t) - 1 
U = U + a^2*(1/dx^2*A1*U + 1/dy^2*U*A2)*dt;
if mod(n,10) == 1;
surf(X,Y,U);
axis([x(1) x(end) y(1) y(end) -1 1]);
getframe;
end
end

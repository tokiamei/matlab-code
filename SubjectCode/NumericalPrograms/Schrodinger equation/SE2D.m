clear;
hb = 1; m = 1;
dx = 0.2;
Lx = 10;
x = -Lx:dx:Lx;
dy = 0.1;
Ly = 10;
y = -Ly:dy:Ly; 
A1 = -2*eye(length(x)) + diag(ones(length(x)-1,1),1) + diag(ones(length(x)-1,1),-1);
A1(1,end) = 1; A1(end,1) = 1; %周期边界条件
A2 = -2*eye(length(y)) + diag(ones(length(y)-1,1),1) + diag(ones(length(y)-1,1),-1);
A2(1,end) = 1; A2(end,1) = 1; %周期边界条件
% m1 = x; m2 = x'.^2; m3 = y; m4 = 2 - y;
dt = 0.00001;
t = 0:dt:2;


[Y,X] = meshgrid(y,x); %x和y调换下位置就可以实现非方矩阵运算
wx = 1;wy = 2;
V = 1/2*m*wx^2*X.^2 + 1/2*m*wy^2*Y.^2;
D = 3;
U0 = exp(-((X).^2+(Y-3).^2)/D^2).*exp(1i*3*X);
U = U0;
% F = zeros(length(x),length(y));

for n = 1:length(t) - 1 
U = U + 1i*hb/2/m*(1/dx^2*A1*U + 1/dy^2*U*A2)*dt + 1/1i/hb*V.*U*dt;
% U(:,1) = m1; U(:,end) = m2; U(1,:) = m3; U(end,:) = m4;
if mod(n,5000) == 1;
surf(X,Y,abs(U.^2));
shading interp;
axis([x(1) x(end) y(1) y(end) 0 3]);

getframe;
end
end

% surf(X,Y,abs(U.^2));
% shading interp;
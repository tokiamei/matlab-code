% ============= 周期结构
% 时间：2022年4月3日10:48:51
hb = 1;
m = 1;
w = 1;
N = 1000;
c = 1;
L = 20*c;
dx = L/N;
% x = -L/2:dx:L/2;
x = 0:dx:L;
% v = 1/2*m*w^2*x.^2;
v0 = 10;
v = v0*cos(2*pi*x/c);
V = diag(v(2:end-1));
A = 1/dx^2*(-2*eye(N-1) + diag(ones(N-2,1),1) + diag(ones(N-2,1),-1));
H = -hb^2/2/m*A + V;
E = eig(H);
plot(E,'.'); %横坐标是n，就是说第n个本征值

% [P,D] = eig(H);
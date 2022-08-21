%%%%%%%%%%%%向后欧拉方法求解GP方程
%%%%%%%%%%%%向后欧拉方法复杂度太高
clear;clc;
g = 1.5; %相互作用强度
dt = 0.001; %time step
T = 10;
t = 0:dt:T;
a = -5; b = 5;
dx = 0.1;
x = a:dx:b;
w = 0;
v = 1/2*w*x.^2;
V = diag(x);
psi = zeros(length(x),length(t));
psi(:,1) = 1/2*x;
% psi(:,1) = exp(-x.^2);
% norm=conj(psi(:,1)')*psi(:,1);
% psi(:,1)=psi(:,1)/norm;

A = -2*diag(ones(length(x),1)) + diag(ones(length(x)-1,1),1) + diag(ones(length(x)-1,1),-1);


for i = 1:length(t)-1
   for j = 1:length(x)
       fh = @(t)dt*(1/2*dx^2*A*t-V*t-g*abs(psi(j,i))^2*t)-psi(j,i)-t;
      psi(j,i+1) = fzero(fh,psi(j,i)); 
   end
end

plot(x,abs(psi(:,i+1)).^2);











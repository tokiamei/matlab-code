%改进的欧拉方法求解一维GP方程
clear;clc;
g = 1.5; %相互作用强度
dt = 0.001; %time step
T = 10;
t = 0:dt:T;
a = -4; b = 5;
dx = 0.1;
x = a:dx:b;
w = 1;
v = 1/2*w*x.^2;
V = diag(x);
psi = zeros(length(x),length(t));
psi(:,1) = 1/2*x;
% psi(:,1) = exp(-x.^2);
% norm=conj(psi(:,1)')*psi(:,1);
% psi(:,1)=psi(:,1)/norm;

A = -2*diag(ones(length(x),1)) + diag(ones(length(x)-1,1),1) + diag(ones(length(x)-1,1),-1);

for n=1:length(t)-1
%    psi(:,n+1) = psi(:,n) + dt*(1/2/dx^2*A*psi(:,n)-V*psi(:,n)-g*abs(psi(:,n).^2.*psi(:,n))); 
    k1 = psi(:,n)+dt*(1/2/dx^2*A*psi(:,n)-V*psi(:,n)-g*abs(psi(:,n).^2).*psi(:,n));
    k2 = psi(:,n)+dt*(1/2/dx^2*A*k1-V*k1-g*abs(k1.^2).*k1);
    psi(:,n+1) = 1/2*(k1+k2);
   
   if mod(n,20)==1
       plot(x,abs(psi(:,n+1)).^2);
       axis([x(1) x(end) 0 5]);
      getframe; 
   end
end
% plot(x,abs(psi(:,n+1)).^2);
% plot(x,abs(psi(:,n+1)).^2);
% hold on;


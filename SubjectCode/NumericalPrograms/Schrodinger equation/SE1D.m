hb = 1;
m = 1;
dt = 0.0005;
t = 0:dt:10;
L = 10;
dx = 0.1;
x = -L:dx:L;
psi = zeros(length(x),length(t));
x0 = 0; Dx = 2;
k0 = 5; 
psi0 = exp(-x'.^2/Dx^2).*exp(1i*k0*x');
psi(:,1) = psi0;
A = -2*eye(length(x)) + diag(ones(length(x)-1,1),1) + diag(ones(length(x)-1,1),-1);
% A(1,end) = 1; %周期边界条件实现的方法
% A(end,1) = 1;
w = 1;
% v = 1/2*m*w^2*x.^2;
v = x;
V = diag(v);
H = -hb^2/(2*m)*A/dx^2 + V;

for n = 1:length(t)-1
    k1 = 1/(1i*hb)*H*psi(:,n);
    k2 = 1/(1i*hb)*H*(psi(:,n)+k1*dt/2);
    k3 = 1/(1i*hb)*H*(psi(:,n)+k2*dt/2);
    k4 = 1/(1i*hb)*H*(psi(:,n)+k3*dt);
    psi(:,n+1) = psi(:,n) + 1/6*dt*(k1+2*k2+2*k3+k4);
%     psi(:,n+1) = psi(:,n) + 1/(1i*hb)*H*psi(:,n)*dt;
%     psi(1,n+1) = 0; %这样就是默认1和end前后还有一个0
%     psi(end,n+1) = 0;
    if mod(n,50) == 1
    plot(x, real(psi(:,n+1)), x, imag(psi(:,n+1)));
%     plot(x,abs(psi(:,n+1).^2))
    axis([x(1) x(end) -1 1.5]);
    getframe;
    end
end


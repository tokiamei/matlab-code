clear;
a = 1;
dx = 0.02;
x = 0:dx:1;
dt = 0.0001;
t = 0:dt:1;
u = zeros(length(x),length(t));
u(:,1) = exp(-20*(x-1/2).^2);
m1 = 0+0*sin(t);
m2 = 1-0*sin(t);
A = -2*eye(length(x)) + diag(ones(length(x)-1,1),1) + diag(ones(length(x)-1,1),-1);

for n = 1:length(t)-1
    u(:,n+1) = u(:,n) + a^2*dt/dx^2*A*u(:,n);
    u(1,n+1) = m1(n+1);
    u(end,n+1) = m2(n+1);
    plot(x,u(:,n),'r');
    axis([x(1) x(end) 0 1]);
    getframe;
end

% [T,X] = meshgrid(t,x); %三维热流
% surf(X,T,u);
% shading interp;

% plot(x,u(:,end));
% hold on;




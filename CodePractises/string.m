dt = 0.01;
t = 0:dt:2;
m = pi;
L = 2;
dx = 0.01;
x = 0:dx:L;
a = 1;

for n = 1:length(t)
   u = sin(m*pi*x/L)*sin(m*pi*a*t(n)/L) + sin(2*m*pi*x/L)*sin(2*m*pi*a*t(n)/L); 
   plot(x,u);
%    axis equal;
   axis([0 L, -2 2]);
   getframe;
end

% plot(x,u);
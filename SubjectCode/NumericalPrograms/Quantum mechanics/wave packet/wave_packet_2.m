clear;
hb = 1;
m = 1;
Lx = 20;
dx = 0.01;
x = -Lx:dx:Lx;
k0 = 3*pi/4;
dk = 0.1;
Lk = 10;
k = -Lk:dk:Lk;
Dk = 0.2;
phik = (pi/2)^(-1/4)*Dk^(-1/2)*exp(-(k-k0).^2/Dk^2);
T = Lx/(hb*k0/m);
dt = 0.05;
t = 0:dt:T;
% w = hb*k.^2/2/m; %É«É¢¹ØÏµ
% w = 5*k;
w = sin(k);

for i = 1:length(t)
   psix = 0;
   for j = 1:length(k)
      psix = psix +  phik(j)*exp(1i*(k(j)*x-w(j)*t(i)))*dk;
      
   end
   
   subplot(2,1,1);
   plot(x, real(psix), x, imag(psix));
   axis([x(1) x(end) -2 2]);
   subplot(2,1,2);
   plot(x,abs(psix.^2));
   axis([x(1) x(end) 0 3]);
   getframe;
end
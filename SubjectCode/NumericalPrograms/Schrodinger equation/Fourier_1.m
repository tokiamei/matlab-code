%%%%%%%%%%%%用波包的构造来理解傅里叶变换以及测不准关系
clear;
Dk = 0.5;
k0 = 5;
dk = 0.02;
Lk = 10;
k = -Lk:dk:Lk;
dx = 0.01;
Lx = 10;
% x = -Lx:dx:Lx;
x = -10;  %%%x较大时，就不是相长干涉，是相消干涉
phik = exp(-(k-k0).^2/Dk^2);
p = phik.*exp(1i*k.*x)*dk;

ps = cumsum(p);
plot(real(ps),imag(ps));
axis equal;
axis([-1 1 -1 1]);
hold on;

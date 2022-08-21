dx = 0.01;
x = 0:dx:10;
% y = @(x)exp(-x.^2);
y = @(x)x;
a = Euler01(y,x,dx);
plot(x,a,'r');
hold on;
plot(x,x,'b')

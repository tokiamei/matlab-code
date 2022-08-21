syms x;
y1 = exp(x^2)*cos(x);
t1 = taylor(y1,x,8,'order',8);
subplot(1,2,1);ezplot(t1,[-5,5]);grid;
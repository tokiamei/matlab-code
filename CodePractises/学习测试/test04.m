syms x;
disp('原函数f');
f = exp(-2*abs(x));
disp('象函数F');
F = fourier(f);
%作图
subplot(1,2,1);ezplot(f);grid;
subplot(1,2,2);ezplot(F);grid;
syms x;
f1 = (1+sin(x))/(1+cos(x));
df1 = diff(f1,x);
figure (1);
subplot(1,2,1);ezplot(f1);grid on;
subplot(1,2,2);ezplot(df1);grid on;

syms x;
y1 = asin(x);
y2 = atan(x);
y3 = atan(x);
y4 = acot(x);

y = [y1,y2;y3,y4];
dy = diff(y);
figure (2);
sublpot(2,2,1);ezplot(y1);grid;
% subplot(2,2,)





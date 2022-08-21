x = 1:0.1:5;
f1 = x.^2;
f2 = x.^3;
normf =  sqrt(f1*f1.'+f2*f2.');
f1 = f1/normf;
f2 = f2/normf;
plot(x,f1,x,f2);
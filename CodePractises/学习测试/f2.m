function F = f2(x0)
%%%%%%%fsolve的用法
x = x0(1);
y = x0(2);
z = x0(3);
f1 = 3*x-cos(y*z)-1/2;
f2 = x^2-81*(y+0.1)^2+sin(z)+1.06;
f3 = exp(-x*y)+20*z+(10*pi-3)/3;
F = [f1,f2,f3];

end
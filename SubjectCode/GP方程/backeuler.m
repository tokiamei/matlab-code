function [T,Y] = backeuler(fun,t,y0,h,dx)

T = t;
M = length(T);
Y(:,1) = y0;
N = length(y0);
A = -2*diag(ones(length(x),1)) + diag(ones(length(x)-1,1),1) + diag(ones(length(x)-1,1),-1);
v = 1/2*w*x.^2;
V = diag(x);

end
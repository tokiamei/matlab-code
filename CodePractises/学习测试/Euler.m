function [t,x] = Euler(fun,t0,tt,x0,N)
h = (tt-t0)/N;
t = t0:h:tt;
t = t';
x(1,:) = x0';
for k = 1:N
    f = feval(fun,t(k),x(k,:));
%     f = f';
    x(k+1,:) = x(k,:) + h*f;
end

end
function gen = Erfen(f,a,b,tol)
% 二分法计算非线性方程
% f 是方程f(x) = 0中的f(x)
%如果输入变量缺省 则默认误差为1E-3
%半独立制作的二分法求解程序

if (nargin == 3)
   tol = 1.0e-3; 
end

% gen = compute_gen(f,a,b,tol);
% 
% 
% function r = compute_gen(f,a,b,tol)
%计算左端点函数值
fa = feval(f,a);
%右端函数值
fb = feval(f,b);
while (fa*fb<0 && abs(a-b)>tol)
    t = (a+b)/2;
    ft = feval(f,t);
    if (fa*ft<=0)
        b = t;
    else
        a = t;   
    end
end
    gen = t;
end



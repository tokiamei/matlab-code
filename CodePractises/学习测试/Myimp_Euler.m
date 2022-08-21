function [t,y] = Myimp_Euler(fun,t,h,y0)
 % 隐式欧拉方法
 % 再跟GF发火，你就是不折不扣的小丑!!!
 T = t;
 N = length(T);
 y(:,1) = y0'; % 初值设为行向量
 m = length(y0);
 
 for i = 1:N-1
     for j = 1:m
         fh = @(x)h*fun(T(i+1),x)-x+y(j,i);
         y(j,i+1) = fzero(fh,y(j,i));

%     y(j,i+1) = fzero(y(j,i+1)-y(j,i)-h*feval(fun,t(j+1),y(:,i+1)),y(j,i));
     end
 end
end
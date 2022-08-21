function [t,x]=GjEuler(fun,t0,tt,x0,N)
h=(tt-t0)/N;
t=t0+[0:N]'*h; % 将时间离散化
x(1,:)=x0';
for i=1:N
    f1=h*feval(fun,t(i),x(i,:));
    f1=f1';
    f2=h*feval(fun,t(i+1),x(i,:)+f1);
    f2=f2';
    x(i+1,:)=x(i,:)+1/2*(f1+f2);
end




end
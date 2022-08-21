function[x,y]=imp_euler(f,a,b,y0,h)
%隐式欧拉格式
%f是待求函数的一阶导形式
%a,b分别是积分上下限
%y0是初始条件y(0)
%h是步长
s = (b-a)/h;
X = zeros(1,s+1);
Y = zeros(1,s+1);
X = a:h:b;
Y(1) = y0;

for i = 1:s
    Y(i+1) = Y(i)+h*feval(f,X(i+1));
end

x = X';
y = Y';
end











function s = tixing(f,a,b,n)
 % 梯形公式求积分
 % a,b为积分区间
 % n为区间划分细度
h = (b-a)/n; 
s = 0;

for k = 1:(n-1)
   x = h*k;
   s = s+feval(f,x);
end
s = h*(feval(f,a)+feval(f,b))/2 + h*s;

end
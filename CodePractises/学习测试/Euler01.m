function s = Euler01(y,x,h,s(0))

for i=1:length(x)-1
   s(i+1) = s(i) +  h*feval(y,x(i+1));
end
end
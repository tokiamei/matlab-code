dx = 0.1;
x = 0:dx:10;
y = zeros(1,length(x));

for i = 1:length(x)-1
   y(i+1) = y(i) + dx*y(i+1);
end
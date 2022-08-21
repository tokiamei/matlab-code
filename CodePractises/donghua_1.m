clear;
t = 0:0.01:10;
x = cos(t);
y = cos(3*t + pi/2);
plot(x, y);

for n = 1:length(t)
   scatter(x(n), y(n));
   hold on;
   scatter(x(n),0);
   scatter(0,y(n));
   plot([-2 2], [0 0]);
   plot([0 0], [-2 2]);
   axis equal;
   axis([-2 2 -2 2]);
  
   M(n) = getframe;
%    hold off;
end

movie(M);
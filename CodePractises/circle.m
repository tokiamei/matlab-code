dt = 0.01;
t = 0:dt:100;
w = 1;
% R = 1 + 4*cos(3*t);
R = 1;  
x = R.*cos(w*t) ;
y = R.*sin(w*t);
plot(x, y);

for n = 1:length(t)
%    x(n) = R(n).*cos(w*t(n));
%    y(n) = R(n).*sin(w*t(n));
   scatter(x(n), y(n));
   hold on;
   axis equal;
axis([-5 5 -5 5]);
M(n) = getframe;

end

% plot(x,y);
% axis equal;

% ================= polarization =============
clear;clc;close;
w = 5;
t = 0:0.1:10;
x = cos(w*t);
y = sin(w*t);
% plot(x,y);

figure(1);
for n = 1:length(t)
   scatter3(x(n),y(n),t(n));
%    plot3(x(n),y(n),t(n));
   hold on;
   axis equal;
   axis([ -1 1  -1 1 0 10]);
%    axis([ -10 10  -10 10 0 10]);
   view(2);
   frame = getframe(gcf); % 制作gif方法
    img =  frame2im(frame);
    [img,cmap] = rgb2ind(img,256);
    if n == 1
        imwrite(img,cmap,'lpolarization.gif','gif','LoopCount',Inf,'DelayTime',0.1);
    else
        imwrite(img,cmap,'lpolarization.gif','gif','WriteMode','append','DelayTime',0.1);
    end
end



% plot3(x,y,t);








clear;
t = 0:0.1:10;

for n = 1:length(t)
    c = 1 + 0.3*cos(2*t(n));
th = 0:pi/100:2*pi;
phi = 0:pi/100:2*pi;
[Th, Phi] = meshgrid(th, phi); 
R1 = 2;
R2 = 1+cos(10*(Phi-Th-t(n)));
X = (R1+R2.*cos(Th)).*cos(Phi)*c;
Y = (R1+R2.*cos(Th)).*sin(Phi)*c;
Z = R2.*sin(Th)*c;
surf(X,Y,Z);
axis equal;
axis([-5 5 -5 5 -5 5]);
shading interp;

frame = getframe(gcf); %%%%%%%%%%%制作gif方法
    img =  frame2im(frame);
    [img,cmap] = rgb2ind(img,256);
    if n == 1
        imwrite(img,cmap,'animation.gif','gif','LoopCount',Inf,'DelayTime',1);
    else
        imwrite(img,cmap,'animation.gif','gif','WriteMode','append','DelayTime',1);
    end
% grid off;
end





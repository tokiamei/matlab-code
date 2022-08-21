clear; %%%%%%%%%%%制作avi动画
t = 0:0.1:10;
aviobj=VideoWriter('donuts.avi');%新建叫example.avi的文件
open(aviobj); %打开
for n = 1:length(t)
    
th = 0:pi/50:2*pi;
phi = 0:pi/50:2*pi;
[Th, Phi] = meshgrid(th, phi); 
R1 = 2;
R2 = 1+0.2*cos(10*(Phi- Th -t(n)));
X = (R1+R2.*cos(Th)).*cos(Phi);
Y = (R1+R2.*cos(Th)).*sin(Phi);
Z = R2.*sin(Th);
surf(X,Y,Z);
axis equal;
axis([-3.5 3.5 -3.5 3.5 -1.5 1.5]);
shading interp;
currFrame = getframe;
writeVideo(aviobj,currFrame);
% grid off;
end

close(aviobj); %关闭
% movie2avi(M, 'new.avi');

clear; %%%%%%%%%%%����avi����
t = 0:0.1:10;
aviobj=VideoWriter('donuts.avi');%�½���example.avi���ļ�
open(aviobj); %��
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

close(aviobj); %�ر�
% movie2avi(M, 'new.avi');

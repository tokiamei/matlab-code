clear;
hb = 1;
m = 1;
wx = 3;
wy = 3;
g = 3;
dx = 0.1;
Lx = 10;
x = -Lx:dx:Lx;
dy = 0.1;
Ly = 10;
y = -Ly:dy:Ly;
[Y,X] = meshgrid(y,x);

dt = 0.001;
t = 0:dt:4;
v0 = 1;
V = v0*1/2*m*(wx^2*X.^2 + wy^2*Y.^2);

Ax = -2*eye(length(x)) + diag(ones(length(x)-1,1),1) + diag(ones(length(x)-1,1),-1);
% Ax(1,end) = 1; Ax(end,1) = 1;
Ay = -2*eye(length(y)) + diag(ones(length(y)-1,1),1) + diag(ones(length(y)-1,1),-1);
% Ay(1,end) = 1; Ay(end,1) = 1;
D = 1;
% psi0 = exp(-(X.^2 + (Y-3).^2)/D^2).*exp(1i*3*X);
psi0 = exp(-(X.^2 + (Y-3).^2)/D^2);
psi = psi0;

for n = 1:length(t)-1
   psi = psi + 1i*hb/2/m*(1/dx^2*Ax*psi + 1/dy^2*psi*Ay)*dt + 1/1i/hb*V.*psi*dt + 1/1i/hb*g*abs(psi.^2).*psi*dt;

% k1 = 1i*hb/2/m*(1/dx^2*Ax*psi + 1/dy^2*psi*Ay) + 1/1i/hb*V.*psi + 1/1i/hb*g*abs(psi.^2).*psi;
% k2 = 1i*hb/2/m*(1/dx^2*Ax*(psi+1/2*k1*dt) + 1/dy^2*(psi+1/2*k1*dt)*Ay) + 1/1i/hb*V.*(psi+1/2*k1*dt) +...
%                                            1/1i/hb*g*abs((psi+1/2*k1*dt).^2).*(psi+1/2*k1*dt);
% k3 = 1i*hb/2/m*(1/dx^2*Ax*(psi+1/2*k2*dt) + 1/dy^2*(psi+1/2*k2*dt)*Ay) + 1/1i/hb*V.*(psi+1/2*k2*dt) +...
%                                            1/1i/hb*g*abs((psi+1/2*k2*dt).^2).*(psi+1/2*k2*dt);
% k4 = 1i*hb/2/m*(1/dx^2*Ax*(psi+k3*dt) + 1/dy^2*(psi+k3*dt)*Ay) + 1/1i/hb*V.*(psi+k3*dt) + 1/1i/hb*g*abs((psi+k3*dt).^2).*(psi+k3*dt);
%    psi = psi + 1/6*(k1 + 2*k2 + 2*k3 + k4)*dt;
   if mod(n,100) == 1
      surf(X,Y,abs(psi.^2));
      shading interp;
      colormap winter;
      colorbar;
      axis([x(1) x(end) y(1) y(end) 0 3]);
      frame = getframe(gcf); %%%%%%%%%%%制作gif方法
%     img =  frame2im(frame);
%     [img,cmap] = rgb2ind(img,256);
%     if n == 1
%         imwrite(img,cmap,'animation.gif','gif','LoopCount',Inf,'DelayTime',1);
%     else
%         imwrite(img,cmap,'animation.gif','gif','WriteMode','append','DelayTime',1);
%     end
   end
end
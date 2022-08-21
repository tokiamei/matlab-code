hb = 1;
m = 1;
wx = 1;
wy = 1;
dx = 0.1;
Lx = 10;
x = -Lx:dx:Lx;
dy = 0.1;
Ly = 10;
y = -Ly:dy:Ly;
[Y,X] = meshgrid(y,x);

dt = 0.0001;
t = 0:dt:2;
V = 1/2*m*(wx^2*X.^2 + wy^2*Y.^2);

Ax = -2*eye(length(x)) + diag(ones(length(x)-1,1),1) + diag(ones(length(x)-1,1),-1);
Ay = -2*eye(length(y)) + diag(ones(length(y)-1,1),1) + diag(ones(length(y)-1,1),-1);
D = 1;
psi0 = exp(-(X.^2 + Y.^2)/D^2).*exp(1i*3*X);
psi = psi0;

for n = 1:length(t)-1
   psi = psi + 1i*hb/2/m*(1/dx^2*Ax*psi + 1/dy^2*psi*Ay)*dt + 1/1i/hb*V.*psi*dt;
   if mod(n,200) == 1
      surf(X,Y,abs(psi.^2));
      shading interp;
      axis([x(1) x(end) y(1) y(end) 0 1]);
      frame = getframe(gcf); %%%%%%%%%%%制作gif方法
    img =  frame2im(frame);
    [img,cmap] = rgb2ind(img,256);
    if n == 1
        imwrite(img,cmap,'animation.gif','gif','LoopCount',Inf,'DelayTime',1);
    else
        imwrite(img,cmap,'animation.gif','gif','WriteMode','append','DelayTime',1);
    end
   end
end
clear;
hb = 1;
m = 1;
wx = 1;
wy = 1;
Lx = 10;
dx = 0.1;
x = -Lx:dx:Lx;
Ly = 10;
dy = 0.1;
y = -Ly:dy:Ly;
[Y,X] = meshgrid(y,x);
dt = 0.01;
t = 0:dt:1;
Ax = -2*eye(length(x)) + diag(ones(length(x)-1,1),1) + diag(ones(length(x)-1,1),-1);
Ay = -2*eye(length(y)) + diag(ones(length(y)-1,1),1) + diag(ones(length(y)-1,1),-1);

V = 1/2*m*(wx^2*X.^2 + wy^2*Y.^2);
% U = zeros(Nx+1,Ny+1);
D = 3;
k0 = 0;
U0 = exp(-(X.^2+Y.^2)/D^2).*exp(1i*k0*X);
U = U0;

for n = 1:dt:length(t)-1
    U = U + dt/1i/hb*(-hb^2/2/m*(1/dx^2*Ax*U + 1/dy^2*U*Ay));
    
    if mod(n,200) == 1
        surf(X,Y,abs(U.^2));
        shading interp;
        axis([x(1) x(end) y(1) y(end) 0 1]);
        getframe;
    end
end




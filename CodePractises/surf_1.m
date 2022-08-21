clear;
x = -4:0.1:4;
y = -5:0.1:5;
[X,Y] = meshgrid(x,y);
% Z = X.^2 + Y.^2;
Z = exp(-X.^2 - Y.^2).*sin(X).*sin(Y);
surf(X,Y,Z);
shading interp;
colormap spring;
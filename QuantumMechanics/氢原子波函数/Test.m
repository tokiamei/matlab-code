[X,Z]=meshgrid(linspace(-30,30,120));
Y=zeros(size(X));

psi=HydWave(4,2,0,X,Y,Z);
surf(X,Z,psi,'EdgeColor','none')
view(2);caxis([0,2])
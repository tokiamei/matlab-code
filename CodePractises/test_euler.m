dx = 0.1;
x = 1:dx:10;
dt = 0.01;
t = 0:dt:2;
psi = zeros(length(x),length(t));
psi(:,1) = x;

for j = 1:length(x)
for i = 1:length(t)-1
   psi(j,i+1) = psi(j,i) + dt*psi(j,i+1); 
end
end










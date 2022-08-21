dt = 0.01;
t = 0:dt:5;
f = @(x)x;
x = 1:10;
psi = zeros(length(x),length(t));
psi(:,1) = x;

for j = 1:length(x)
    for i=1:length(t)-1
   fzero(psi(j,i+1) = psi(j,i) + dt*psi(j,i+1)); 
    end
end
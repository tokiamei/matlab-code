clear;clc;
k0 = 3; %激光波数
L = 10;
a = -L;b = L;
dx = 0.1;
x = a:dx:b;
T = 40;
dt = 0.005;
t = 0:dt:T;

o = 0.02;
delta = 0; %失谐
b11 = o; b12 = 0.8*o; b21 = 0.8*o; b22 = o;
% b11 = 0; b12 = 0; b21 = 0; b22 = 0;
omiga = 200; %拉曼耦合强度
phi1 = zeros(length(x), length(t));
phi2 = zeros(length(x), length(t));

phi1(:,1) = rand(length(x),1)+1i*rand(length(x),1);
phi2(:,1) = rand(length(x),1)+1i*rand(length(x),1);

norm = sum(abs(phi1(:,1).^2)+abs(phi2(:,1).^2));
phi1(:,1) = phi1(:,1)/norm;
phi2(:,1) = phi2(:,1)/norm;

%拉普拉斯算符矩阵
A = [-2*diag(ones(length(x),1)) + diag(ones(length(x)-1,1),1) + diag(ones(length(x)-1,1),-1)]/dx^2;
%一阶导数算符矩阵
B = [diag(ones(length(x)-1,1),1) - diag(ones(length(x)-1,1),-1)]/2/dx;
%势场
ita = 1;
% v = 1/2*ita*x.^2;
v = ita*400*ones(1,length(x));
v(20:180) = 0;
V = diag(v);
for i = 1:length(t)-1
    % 龙格库塔法
    m1 = [1/2*A*phi1(:,i)-V*phi1(:,i)-1i*k0*B*phi1(:,i)-delta/2*phi2(:,i)-b11*abs(phi1(:,i).^2).*phi1(:,i)+...
                      -b12*abs(phi2(:,i).^2).*phi1(:,i)-omiga/2*phi2(:,i)];
    m2 = [1/2*A*(phi1(:,i)+1/2*m1*dt)-V*(phi1(:,i)+1/2*m1*dt)-1i*k0*B*(phi1(:,i)+1/2*m1*dt)-delta/2*phi2(:,i)-b11*abs((phi1(:,i)+1/2*m1*dt).^2).*(phi1(:,i)+1/2*m1*dt)+...
                      -b12*abs(phi2(:,i).^2).*(phi1(:,i)+1/2*m1*dt)-omiga/2*phi2(:,i)];
    m3 = [1/2*A*(phi1(:,i)+1/2*m2*dt)-V*(phi1(:,i)+1/2*m2*dt)-1i*k0*B*(phi1(:,i)+1/2*m2*dt)-delta/2*phi2(:,i)-b11*abs((phi1(:,i)+1/2*m2*dt).^2).*(phi1(:,i)+1/2*m2*dt)+...
                      -b12*abs(phi2(:,i).^2).*(phi1(:,i)+1/2*m2*dt)-omiga/2*phi2(:,i)];
    m4 = [1/2*A*(phi1(:,i)+m3*dt)-V*(phi1(:,i)+m3*dt)-1i*k0*B*(phi1(:,i)+m3*dt)-delta/2*phi2(:,i)-b11*abs((phi1(:,i)+m3*dt).^2).*(phi1(:,i)+m3*dt)+...
                      -b12*abs(phi2(:,i).^2).*(phi1(:,i)+m3*dt)-omiga/2*phi2(:,i)];
   

    n1 = [1/2*A*phi2(:,i)-V*phi2(:,i)+1i*k0*B*phi2(:,i)+delta/2*phi2(:,i)-b21*abs(phi1(:,i).^2).*phi2(:,i)+...
                      -b22*abs(phi2(:,i).^2).*phi2(:,i)-omiga/2*phi1(:,i)];
    n2 = [1/2*A*(phi2(:,i)+1/2*n1*dt)-V*(phi2(:,i)+1/2*n1*dt)+1i*k0*B*(phi2(:,i)+1/2*n1*dt)+delta/2*(phi2(:,i)+1/2*n1*dt)-b21*abs(phi1(:,i).^2).*(phi2(:,i)+1/2*n1*dt)+...
                      -b22*abs((phi2(:,i)+1/2*n1*dt).^2).*(phi2(:,i)+1/2*n1*dt)-omiga/2*phi1(:,i)];
    n3 = [1/2*A*(phi2(:,i)+1/2*n2*dt)-V*(phi2(:,i)+1/2*n2*dt)+1i*k0*B*(phi2(:,i)+1/2*n2*dt)+delta/2*(phi2(:,i)+1/2*n2*dt)-b21*abs(phi1(:,i).^2).*(phi2(:,i)+1/2*n2*dt)+...
                      -b22*abs((phi2(:,i)+1/2*n2*dt).^2).*(phi2(:,i)+1/2*n2*dt)-omiga/2*phi1(:,i)];
    n4 = [1/2*A*(phi2(:,i)+1/2*n1*dt)-V*(phi2(:,i)+n3*dt)+1i*k0*B*(phi2(:,i)+n3*dt)+delta/2*(phi2(:,i)+n3*dt)-b21*abs(phi1(:,i).^2).*(phi2(:,i)+n3*dt)+...
                      -b22*abs((phi2(:,i)+n3*dt).^2).*(phi2(:,i)+n3*dt)-omiga/2*phi1(:,i)];
  
   phi1(:,i+1) = phi1(:,i) + 1/6*(m1+m4+2*m2+2*m3)*dt;
   phi2(:,i+1) = phi2(:,i) + 1/6*(n1+n4+2*n2+2*n3)*dt;
   
    
end

norm = sqrt(sum(abs(phi1(:,end)).^2+abs(phi2(:,end)).^2));
phi1(:,end) = phi1(:,end)/norm;
phi2(:,end) = phi2(:,end)/norm;
plot(x,abs(phi1(:,end)).^2,'r',x,abs(phi2(:,end)).^2,'g');
% axis([a b 0 2e-5]);
% axis([a b 0 1]);












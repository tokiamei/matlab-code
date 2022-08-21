%%%%%%%%%%%%%%%GP������ֵ���
% close;
m = 1;
hb = 1;
a = -2;
L = 5;
dx = 0.09;
x = -L:dx:L;
dt = 0.00001;
t = 0:dt:1;
psi = zeros(length(x),length(t));
nrho = zeros(length(x),length(t));
x0 = 0; %��˹����λ��x = 0;
w = 1; %���
k0 = 10;
% psi(:,1) = exp(-(x'-x0).^2/w^2).*exp(1i*k0*x');
psi(:,1) = exp(-(x'-x0).^2/w^2);
nrho(:,1) = abs(psi(:,1).*psi(:,1))/(psi(:,1)'*psi(:,1));
% psi(:,1) = sin(1/2*x);

A = -2*eye(length(x)) + diag(ones(length(x)-1,1),1) + diag(ones(length(x)-1,1),-1);
om = 2;
v =1/2*m*om^2*x.^2;
% v = 5;
V = diag(v);
H = -hb^2/(2*m)*A/dx^2 + V;

for n = 1:length(t)-1
%      k1 = 1/(1i*hb)*H*psi(:,n);
%     k2 = 1/(1i*hb)*H*(psi(:,n)+k1*dt/2);
%     k3 = 1/(1i*hb)*H*(psi(:,n)+k2*dt/2);
%     k4 = 1/(1i*hb)*H*(psi(:,n)+k3*dt);
%     psi(:,n+1) = psi(:,n) + 1/6*dt*(k1+2*k2+2*k3+k4);

psi(:,n+1) = psi(:,n) + dt/(1i*hb)*(-hb^2/(2*m)*A/dx^2*psi(:,n)+V*psi(:,n)+a*abs(psi(:,n).^2).*psi(:,n));
nrho(:,n+1) = abs(psi(:,n+1).*psi(:,n+1))/(psi(:,n+1)'*psi(:,n+1));
%     psi(1,n+1) = 0;
%     psi(end,n+1) = 0;
%     if mod(n,200) == 0
% %     plot(x,real(psi(:,n+1)), 'r',x,imag(psi(:,n+1)), 'b');
%     plot(x, abs(psi(:,n+1).*psi(:,n+1)));
% %     axis equal;
%     axis([x(1) x(end) -2 4]);
%     getframe;
%     end
end
% 



%%����ά����
f1 = figure(1);
[T,X] = meshgrid(t,x);
% surf(X,T,nrho);
surf(X,T,real(psi));
title('BEC��̬��������ʱ����ݻ�');
xlabel('����','color', 'r'); ylabel('ʱ��', 'color', 'b'); zlabel('̬�ܶ�', 'color','g');
shading interp;

%������άͼ��
path = 'D:\matlab figure\GPE\��̬������\';
str = strcat(path,'������̬�ܶ���άͼ��','.jpg');
% saveas(f1, str);

%%%%a>0ʱ���������������Ӽ��໥����ǿ�ȵı仯
f2 = figure(2);
plot(x, real(psi(:,n+1)));
% plot(x,nrho(:,n+1))
title('��̬������̬�ܶ�');
xlabel('����'); ylabel('̬�ܶ�');
axis([0 4 0 0.5]);
hold on;

































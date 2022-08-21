%%%%%%%%%%%%�ò����Ĺ�������⸵��Ҷ�任�Լ��ⲻ׼��ϵ
clear;
Dk = 0.5;
k0 = 5;
dk = 0.02;
Lk = 10;
k = -Lk:dk:Lk;
dx = 0.01;
Lx = 10;
% x = -Lx:dx:Lx;
x = -10;  %%%x�ϴ�ʱ���Ͳ����೤���棬����������
phik = exp(-(k-k0).^2/Dk^2);
p = phik.*exp(1i*k.*x)*dk;

ps = cumsum(p);
plot(real(ps),imag(ps));
axis equal;
axis([-1 1 -1 1]);
hold on;

%%%%%%%%%%%%%%构造波包，利用具有给定动量的波函数进行线性组合，得到局域在一定范围里的波包
          %%%也就是把平面波进行线性组合得到一个局域的波函数
clear;
Dk = 5; %%%%%%%波包的宽度，胖的程度，波包窄说明动量是确定的，波包宽说明动量是不确定的，波矢的“不确定度”
        %%%%%%%%%%动量空间宽度越宽，坐标空间函数越窄，坐标和动量之间的不确定度成反比例
dk = 0.001;
Lk = 10;
k = -Lk:dk:Lk;
k0 = 5;
% phik = exp(-k.^2/Dk^2); %%%%%表示不同的平面波在psi中的权重，0时候占的权重最大
phik = exp(-(k-k0).^2/Dk^2);
% plot(k, phik, k, phik1);
% hold on;

dx = 0.01;
x = -5:dx:5;
psix = 0;
for n = 1:length(k)
   psix = psix + 1/(2*pi)^1/2*phik(n)*exp(1i*k(n)*x)*dk; 
end

subplot(3,1,1);
plot(x, real(psix), x, imag(psix));
subplot(3,1,2);
plot(x, abs(psix.*psix));
subplot(3,1,3);
plot(k, abs(phik.^2));

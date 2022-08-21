dx=pi/60;
col=0:dx:pi;
az=0:dx:2*pi;
[phi,theta]=meshgrid(az,col);
% ============== 计算 l=3 的网格上的 P_{l}^{m}(\cos \theta) ===========
l=3;
Plm=legendre(l,cos(theta));
% ============== 由于 legendre 为 m 的所有值计算答案，因此 Plm 会包含一些额外的函数值。 ====
% ============== 提取 m=2 的值并丢弃其余值。 ==========================
% ============== 使用 reshape 函数将结果定向为与 phi 和 theta 具有相同大小的矩阵。
m=2;
if l~=0 % 这里写条件判断的原因是为了以后方便改参数
    Plm=reshape(Plm(m+1,:,:),size(phi)); % 写size（theta）也可以
end
% ============== Plm = Plm(m+1,:,:); % 不能这么写，这样写Plm是一个三维矩阵，所以必须用resahpe
% 计算Y_3^2的球谐函数值
a=(2*l+1)*factorial(l-m);
b=4*pi*factorial(l+m);
C=sqrt(a/b);
Ylm=C.*Plm.*exp(1i*m*phi);
% ============== 球面坐标转换为笛卡尔坐标并绘图 ==============================
% surf(phi,theta,abs(real(Ylm)));
[Xm,Ym,Zm]=sph2cart(phi,pi/2-theta,abs(real(Ylm)));
surf(Xm,Ym,Zm)
title('$Y_3^2$ spherical harmonic','interpreter','latex')

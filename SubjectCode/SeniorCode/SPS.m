%% 参数设定
% 失谐为0时，色散曲线随拉曼耦合的变化
ky=0; % [momentum along y direction]
kz=0; % [momentum along z direction]
k0=2; % [wave vector of laser]
omiga1=1; % [SOC strength]
         % 以omiga = 4Er = 8为分界线；omiga < 4Er = 8时有简并基态，omiga > 4Er = 8只有一个基态
omiga2=2;
omiga3=4;
omiga4=5;
omiga5=10;

m=1; % [atom mass]
Er=(k0)^2/2; % [recoil energy]   Er=2 
%% 绘制图像
kx=linspace(-3,3); 
EU1=(kx.*kx+(k0)^2)/(2*m) + sqrt(kx.*kx*(k0)^2/m^2+omiga1^2/4);  % helicity +
EU2=(kx.*kx+(k0)^2)/(2*m) + sqrt(kx.*kx*(k0)^2/m^2+omiga2^2/4);
EU3=(kx.*kx+(k0)^2)/(2*m) + sqrt(kx.*kx*(k0)^2/m^2+omiga3^2/4);
EU4=(kx.*kx+(k0)^2)/(2*m) + sqrt(kx.*kx*(k0)^2/m^2+omiga4^2/4);
EU5=(kx.*kx+(k0)^2)/(2*m) + sqrt(kx.*kx*(k0)^2/m^2+omiga5^2/4);

ED1=(kx.*kx+(k0)^2)/(2*m) - sqrt(kx.*kx*(k0)^2/m^2+omiga1^2/4);  % helicity _
ED2=(kx.*kx+(k0)^2)/(2*m) - sqrt(kx.*kx*(k0)^2/m^2+omiga2^2/4);
ED3=(kx.*kx+(k0)^2)/(2*m) - sqrt(kx.*kx*(k0)^2/m^2+omiga3^2/4);
ED4=(kx.*kx+(k0)^2)/(2*m) - sqrt(kx.*kx*(k0)^2/m^2+omiga4^2/4);
ED5=(kx.*kx+(k0)^2)/(2*m) - sqrt(kx.*kx*(k0)^2/m^2+omiga5^2/4);
% plot(kx,EU1,kx,ED1,'DisplayName','omiga=1')
plot(kx,EU1,kx,EU2,kx,EU3,kx,EU4,kx,EU5,kx,ED1,kx,ED2,kx,ED3,kx,ED4,kx,ED5,'linewidth',1);
set(gca,'Fontsize',16,'Fontname','Times');
xlabel('kx','Interpreter','latex');
ylabel('E','Interpreter','latex');
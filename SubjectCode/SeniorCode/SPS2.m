%% 参数设定
% 在拉曼耦合强度较小时，能量随着失谐变换的曲线
ky=0; % [momentum along y direction]
kz=0; % [momentum along z direction]
k0=2; % [wave vector of laser]
omiga=1; % [SOC strength]
         % 以omiga = 4Er = 8为分界线；omiga < 4Er = 8时有简并基态，omiga > 4Er = 8只有一个基态
delta1=1; % detuning
delta2=2;
delta3=3;
delta4=4;
m=1; % [atom mass]
Er=(k0)^2/2; % [recoil energy]   Er=2 
%% 绘制图像
kx=linspace(-5,5); 
E0=(kx.*kx+(k0)^2)/(2*m) + sqrt(kx.*kx*(k0)^2/m^2+kx*k0*delta1/m+delta1^2/4+omiga^2/4);  % helicity +
E1=(kx.*kx+(k0)^2)/(2*m) - sqrt(kx.*kx*(k0)^2/m^2+kx*k0*delta1/m+delta1^2/4+omiga^2/4);  % helicity _
E2=(kx.*kx+(k0)^2)/(2*m) - sqrt(kx.*kx*(k0)^2/m^2+kx*k0*delta2/m+delta2^2/4+omiga^2/4);  % helicity _
E3=(kx.*kx+(k0)^2)/(2*m) - sqrt(kx.*kx*(k0)^2/m^2+kx*k0*delta3/m+delta3^2/4+omiga^2/4);  % helicity _
E4=(kx.*kx+(k0)^2)/(2*m) - sqrt(kx.*kx*(k0)^2/m^2+kx*k0*delta4/m+delta4^2/4+omiga^2/4);  % helicity _
plot(kx,E1,kx,E2,kx,E3,kx,E4,'linewidth',1);
set(gca,'Fontsize',16,'Fontname','Times');
xlabel('kx');
ylabel('E');
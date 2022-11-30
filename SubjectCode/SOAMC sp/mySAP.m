 ================== 离散化 =========================
LMin = 0; LMax = 100; % rho 的范围
step = 200; % 200 个离散点
rho = linspace(LMin, LMax, step); % rho 的取值
dRho = rho(2) - rho(1); % 相邻 两个离散点的间距

n1 = 1; n2 = -1; % phase winding



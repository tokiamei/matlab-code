w = 1; % trapping frequency
amu = 1.66e-27;
mA = 87*amu; mB = 23*amu; % 87Rb 23Na
hbar = 1.054e-34;
xs = sqrt(hbar/mA/w);
x = 500e-6/xs; % 量纲化的坐标大小

aA = 3.4e-7;
gA = 4*pi*aA/xs;
lamda = 532e-9;
k = 2*pi/lamda;








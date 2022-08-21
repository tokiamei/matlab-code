%傅里叶逆变换
syms w;
f = 1/2*(dirac(w-2)+dirac(w+2));
F = ifourier(f);
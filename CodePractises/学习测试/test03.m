syms x y;
z = exp(x^2);
int1 = int(z,y,x^3,x);
s = int(int1,x,0,1);
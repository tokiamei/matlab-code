function y = if_test(x)
%测试if
if (x>2)
    y = exp(x);
elseif (-2<x && x<2)
        y = x^2+2;
else
    y = -x^2;


end
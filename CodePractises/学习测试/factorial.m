function f = factorial(n)
if (n <= 0)
    f =  1;
else
    f = factorial(n - 1) * n;

end
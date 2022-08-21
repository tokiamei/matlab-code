function f = Fibonacci(n)
if n == 1 || n == 2
    f = 1;
else
    f = Fibonacci(n - 2) + Fibonacci(n - 1);
end
end
function y = recur01(func, n)
% ================== 测试递归 =================

if (n == 1)
   y = 2; 
else
    y = func(n) + recur01(func, n-1);
end

end
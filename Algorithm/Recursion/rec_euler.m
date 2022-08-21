function y = rec_euler(func, y0, t)
% =========================================
% ======== 本程序试图利用递归来计算欧拉公式
% h = t(2) - t(1); % 步长
% n = length(t);
% y = zeros(n, 1);
h = 1;
y = 0;
if (t == 1)
    y = y0;
else
    y = y + h*rec_euler(func, y0, t-1);
end
end
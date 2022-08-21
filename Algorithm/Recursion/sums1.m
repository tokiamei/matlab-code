function s = sums1(n)
if (n == 1)
   s = 1; 
else
    s = sums1(n-1) + n;
end
end
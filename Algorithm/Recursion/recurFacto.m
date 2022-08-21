function cnt = recurFacto(n)
cnt = 0;
if (n <= 0) 
    cnt = 1;
end
for i = 1:n
    cnt = cnt + recurFacto(n -1);
end
end
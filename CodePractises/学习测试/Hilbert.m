function a = Hilbert(n)
for i = 1:n
   for j = 1:n
      a(i,j) = 1/(i+j-1); 
   end
end
function M = map_cell_types(S)

% S is the number of genotypes in the system.
% The output M is an array with first row giving the row number of the array. The second and third rows give information on which genotype each row of x (each ODE) contains. If a row of x represents a cell with genotype 1 fused with genotype 2, this corresponding row of M will contain a 1 and 2 in its 2nd and 3rd columns.   


k=1;
M=zeros(1,3);

for i=1:S
  
  M(i,1)=i;
  M(i,2)=i;
  
end


for i=1:S
   for j=1:S

       if j<i
           continue
       else
       end

   M(k+S,1)=k+S;
   M(k+S,2)=i;
   M(k+S,3)=j;
   k=k+1;
   end
end

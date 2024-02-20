function f = Survival_function(mgenotypes,alphagenotypes,x,S,beta,Cz,Cp)

% Simulates the frequency after each successive generation.
% x is the output of the function "Fertilisation_Kinetics" which represents the number of cells of each type e.g. number of unfused cells of genotype 1, number of genotype 1 cells that have fused with genotype 3 etc.
% \beta is the resistance to survival and C is the fusion cost as defined in the text.
% mgenotypes and alphagenotypes are the mass and fusion rates of each genotype.

Mv=map_cell_types(S);  % writes the array Mv which is used to determine what each column of x (each ODE) represents. e.g. when S=2, x(4) represents F^{1}_{2} the number of genotype 1 cells that've fused with genotype 2 cells.
f=zeros(S,1);

for k=1:S          % This for loop with index k calculates the fitness of each genotype in the system i.e. w_k for k in [1,S] as defined in the text.
   wk=0; 
   
        %------------------------------------------------------------------------------------------------------------------------------------------------------------------
        %%%This block of code computes the absolute fitness of the k-th genotype.
   
   
        Mk = Mv(any(Mv==k,2),:);        % Mk outputs an array showing only the columns of x (which ODE) that contains genotype k cells.     
           for ii=1:length(Mk(:,1))     % This for loop with index ii calculates the absolute fitness of the k-th genotype.

              if Mk(ii,3)==0 || Mk(ii,2)==0       % checks if the ii-th column of Mk represents a zygote, if so we multiply by (1-C).
              wk=wk+x(length(x(:,1)),Mk(ii,1))*exp(-beta/( mgenotypes(Mk(ii,2)+1) + mgenotypes(Mk(ii,3)+1) ) )*(1-Cp);    % any rows with k in their 2nd or 3rd columns. the +1 is here because the first entry of mgenotypes is set to 0 for a purpose. If an entry in Mk is 0, then it adds mgenotypes(1).  
              else
              wk=wk+x(length(x(:,1)),Mk(ii,1))*exp(-beta/( mgenotypes(Mk(ii,2)+1) + mgenotypes(Mk(ii,3)+1) ) )*(1-Cz);
              end
       
           end 
         %---------------------------------------------------------------------------------------------------------------------------------------------------------------
      
      w(k)=wk;                   % writes the absolute fitness of k-th genotype into a vector w. 
end
w(S+1:end)=[];


for i=1:length(w)                % frequency of each genotype in the subsequent generation
f(i)=w(i)/sum(w);   
end

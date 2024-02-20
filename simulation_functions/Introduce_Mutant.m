function [mgenotypes,alphagenotypes,f]=Introduce_Mutant(mgenotypes,alphagenotypes,f,f0,delta,alphamax)

% This code simulates the coevolutionary dynamics between mass m and fusion
% rate \alpha over 1 mutations. The difference between this and "Evol_Dynamics_Mass_Alpha" is that the evolutionary dynamics can be started from any number of genotypes. 
% This code is used for simulating
% "Evol_Dynamics_Mass_Alpha_Switching_Gillespie".

% Parameters: mgenotypes and alphagenotypes - the mass and fusion rate of all the genotypes.
%             f - initial frequency of all the genotypes.
%             f0 - initial frequency of rare mutant.

% This block of code is for initialisation.
%-------------------------------------------------------------

if length(mgenotypes)~=length(alphagenotypes)
error('length of mgenotypes not equal to length of alphagenotypes')
end 

if length(f)~=length(alphagenotypes)-1
error('length of f not equal to length of alphagenotypes-1')
end


%-------------------------------------------------------------




    
       keepgenotypes=0;     % this bit is the part that picks who mutates
       r=rand;
       prob = f;
       if any(f<f0)
             index_to_exclude = f<f0; % 
             logical_index = true(size(f)); % 
             logical_index(index_to_exclude) = false; % 
             new_vector = f(logical_index); % 
             q= sum(r >= cumsum([0;new_vector/sum(new_vector)])); %fine here, no error
       else
       q= sum(r >= cumsum([0;prob/sum(f)]));                                % q is the index of 
       end

       
       if randn<0                                                           % if randn<0 mutate in mass, else mutate in \alpha.

                mgenotypes(length(mgenotypes)+1)=mgenotypes(q+1)+delta*sign(randn);
         
                if mgenotypes(end)<delta
                mgenotypes(end)=delta;    
                end 
                alphagenotypes(length(alphagenotypes)+1)=alphagenotypes(q+1);              % IMPORTANT NOTE! the alpha value of the mutant should be the alpha value of the genotype that underwent mutation in mass! Hence the "alphagenotypes(q+1)"
                S=length(mgenotypes)-1;
   
       else
                           
                alphagenotypes(length(alphagenotypes)+1)=alphagenotypes(q+1)+delta*sign(randn);

                if alphagenotypes(end)<0
                alphagenotypes(end)=0;    
                end

                if alphagenotypes(end)>=alphamax
                alphagenotypes(end)=alphamax;    
                end

                mgenotypes(length(mgenotypes)+1)=mgenotypes(q+1);                         % IMPORTANT NOTE! the m value of the mutant should be the m value of the genotype that underwent mutation in alpha! Hence the "mgenotypes(q+1)"
                S=length(alphagenotypes)-1;   

       end
 
           
             % This block of code checks if there's any duplicate genotypes.
             %-------------------------------------------------------------
             aa=[mgenotypes' alphagenotypes'];   
             [u,I,J] = unique(aa, 'rows', 'first');
             hasDuplicates = size(u,1) < size(aa,1);
             ixDupRows = setdiff(1:size(aa,1), I);
             dupRowValues = aa(ixDupRows,:);
             duprow=0;

             for i=1:length(aa(:,1))
             duprow(i)=isequal(aa(i,:),dupRowValues);
             end
             
             if hasDuplicates>0
             f(find(duprow,1, 'first')-1)=f(find(duprow,1, 'first')-1) + f0;
             f(q)=f(q)-f0;
             mgenotypes(max(ixDupRows))=[];
             alphagenotypes(max(ixDupRows))=[];
             S=length(mgenotypes)-1;
             else             
             f(end+1)=f0;
             f(q)=f(q)-f0;
             end
             if f(q)<10^(-3)   
                 f(q)=[];      
                 mgenotypes(q+1)=[];   
                 alphagenotypes(q+1)=[];   
                 S=length(mgenotypes)-1;   
             end
             %-------------------------------------------------------------

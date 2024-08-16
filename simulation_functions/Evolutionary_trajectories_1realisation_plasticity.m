function [m,alpha,m2,alpha2]=Evolutionary_trajectories_1realisation_plasticity(m0,alpha0,A,M,T,Cz,Cp,beta1,beta2,lambda21,lambda12,mu,NEVOL,f0,delta,alphamax)

% This code simulates the coevolutionary dynamics between mass m and fusion
% rate \alpha over NEVOL mutations for a system undergoing bet-hedging in a
% stochastically switching environment. Here, the outputs genotypeDATA_m and genotypeDATA_alpha gives the set of all genotypes in one simulation of the coevoultionary dynamics. This enables us to
% visualise evolutionary branching in both mass and \alpha.

% As is the case for the code "Evol_Dynamics_Mass_Alpha_EvolBranching", branching is visualised by plotting up the set of all genotypes using the
% functions "Evol_Branching_plots" or "Evol_Branching_plots_Coevolution".


% Parameters: m0 - initial mass
%             \alpha_0 - initial fusion rate
%             f0 - initial frequency of rare mutant
%             mu - mutation rate.
%             \delta - mutational stepsize
%             beta1 - resistance to survival in Env1 environment
%             beta2 - resistance to survival in Env2 environment
%             lambda21 - rate of switching from Env2 to Env1 environment
%             lambda12 - rate of switching from Env1 to Env2 environment


% This block of code is for initialisation.
%-------------------------------------------------------------
Nmut=1;
Nev=1;

if randn<0
mgenotypes=[0,m0,m0+delta*sign(randn)];
alphagenotypes=[0,alpha0,alpha0];
mgenotypes2=[0,m0,m0+delta*sign(randn)];
alphagenotypes2=[0,alpha0,alpha0];
else
mgenotypes=[0,m0,m0];
mgenotypes2=[0,m0,m0];
if alpha0<delta
alphagenotypes=[0,alpha0,alpha0+delta]; 
alphagenotypes2=[0,alpha0,alpha0+delta];
else
alphagenotypes=[0,alpha0,alpha0+delta*sign(randn)];
alphagenotypes2=[0,alpha0,alpha0+delta*sign(randn)];
    if alphagenotypes(end)>alphamax
    alphagenotypes(end)=alphamax;
    end
    if alphagenotypes2(end)>alphamax
    alphagenotypes2(end)=alphamax;
    end

end
end


f=[1-f0,f0];
f2=[1-f0,f0];
S=2;
S2=2;
mcur=0;
alphacur=0;
mcur2=0;
alphacur2=0;


random2=geornd(   (mu+lambda21)*ones(1,7500000)      )+ones(1,7500000) ;
random1=geornd(   (mu+lambda12)*ones(1,7500000)      )+ones(1,7500000) ;

m=zeros(1,NEVOL);
alpha=zeros(1,NEVOL);
m2=zeros(1,NEVOL);
alpha2=zeros(1,NEVOL);

Env2=1; Env1=0; beta=beta2; tnext=random2(Nev); 


%-------------------------------------------------------------



while Nmut<=NEVOL
 keepgenotypes=0;
 keepgenotypes2=0;

 if Env2==1
     lambda=lambda21;
 else
     lambda=lambda12;
 end

   if rand<mu/(mu+lambda)


       if Env1==1
       f=f(:);       
       [mgenotypes,alphagenotypes,f]=Introduce_Mutant(mgenotypes,alphagenotypes,f,f0,delta,alphamax);
       S=length(mgenotypes)-1; 
       elseif Env2==1
       f2=f2(:);       
       [mgenotypes2,alphagenotypes2,f2]=Introduce_Mutant(mgenotypes2,alphagenotypes2,f2,f0,delta,alphamax);
       S2=length(mgenotypes2)-1; 
       end


       
       mutation=1; switching=0;
       Nmut=Nmut+1;       
 
       fprintf('Processing %d...',Nmut);
   else
   
       if Env2==1
           Env2=0; Env1=1; tnext=random1(Nev); beta=beta1;
       else
           Env2=1; Env1=0; tnext=random2(Nev); beta=beta2;
       end

       mutation=0; switching=1;
   end


   if Env1==1
   f = Invasion_dynamics(mgenotypes,alphagenotypes,S,beta,Cz,Cp,tnext,T,A,M,f);
   elseif Env2==1
   f2 = Invasion_dynamics(mgenotypes2,alphagenotypes2,S2,beta2,Cz,Cp,tnext,T,A,M,f2);
   end


             if Env1==1

             % Block of code that Discards genotypes with frequency below
             % 10^(-3) in environment 1
             %-------------------------------------------------------------
             for i=1:length(f)
                if f(i)<10^(-3)
                keepgenotypes(i,:)=0;
                else
                keepgenotypes(i,:)=1;    
                end
             end

             f(keepgenotypes==0)=[];
             f=f/sum(f);
             mgenotypes([1; keepgenotypes]==0)=[];
             alphagenotypes([1; keepgenotypes]==0)=[];
             S=length(mgenotypes)-1;
            %--------------------------------------------------------------


             elseif Env2==1
            % Block of code that Discards genotypes with frequency below
            % 10^(-3) in environment 2
            %-------------------------------------------------------------

             for i=1:length(f2)
                if f2(i)<10^(-3)
                keepgenotypes2(i,:)=0;
                else
                keepgenotypes2(i,:)=1;    
                end
             end

             if length(keepgenotypes2)~=length(f2)
                 keyboard
             end

             f2(keepgenotypes2==0)=[];
             f2=f2/sum(f2);
             mgenotypes2([1; keepgenotypes2]==0)=[];
             alphagenotypes2([1; keepgenotypes2]==0)=[];
             S2=length(mgenotypes2)-1;
             %-------------------------------------------------------------
             end

             if mutation==1

                if Env1==1
                % This block of code calculates the mean mass and alpha of the
                % population in environment 1 before the next mutation occurs.
                %-------------------------------------------------------------
                for i=1:length(f)        
                mcur=mcur+(f(i)/sum(f))*mgenotypes(i+1);
                alphacur=alphacur+(f(i)/sum(f))*alphagenotypes(i+1);
                end       
                m(Nmut-1)=mcur;
                alpha(Nmut-1)=alphacur;
                   if Nmut>2
                   m2(Nmut-1)=m2(Nmut-2);
                   alpha2(Nmut-1)=alpha2(Nmut-2);
                   end

                mcur=0;
                alphacur=0;
                %-------------------------------------------------------------
                elseif Env2==1

                % This block of code calculates the mean mass and alpha of the
                % population in environment 2 before the next mutation occurs.
                %-------------------------------------------------------------
                for i=1:length(f2)        
                mcur2=mcur2+(f2(i)/sum(f2))*mgenotypes2(i+1);
                alphacur2=alphacur2+(f2(i)/sum(f2))*alphagenotypes2(i+1);
                end       
                m2(Nmut-1)=mcur2;
                alpha2(Nmut-1)=alphacur2;
                      if Nmut>2
                      m(Nmut-1)=m(Nmut-2);
                      alpha(Nmut-1)=alpha(Nmut-2);
                      end
                mcur2=0;
                alphacur2=0;
                %-------------------------------------------------------------
                end

              end
  
Nev=Nev+1;

end

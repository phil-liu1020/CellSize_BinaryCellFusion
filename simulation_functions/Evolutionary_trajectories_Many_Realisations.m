function [m,alpha]=Evolutionary_trajectories_Many_Realisations(number_of_realisations,m0,alpha0,A,M,T,Cz,Cp,beta1,beta2,lambda12,lambda21,mu,NEVOL,f0,delta,alphamax, return_genotypes )
                                                         
% Obtains multiple realisations of the simulation for the coevolutionary dynamics between mass and \alpha in a system undergoing bet-hedging in a stochastically switching environment.
% Parameters: Nrealz - number of realisations
%             f0 - initial frequency of rare mutant

m=zeros(number_of_realisations,0);
alpha=zeros(number_of_realisations,0);


    parfor i=1:number_of_realisations
    
    [~,~,m(i,:),alpha(i,:)]=Evolutionary_trajectories_1realisation_BH(m0,alpha0,A,M,T,Cz,Cp,beta1,beta2,lambda21,lambda12,mu,NEVOL,f0,delta,alphamax,return_genotypes);
                                                                     
    end

assignin('base','m',m)
assignin('base','alpha',alpha)

function [genotypeDATA_m,genotypeDATA_alpha,m,alpha]=Evolutionary_trajectories(number_of_realisations,m0,alpha0,A,M,T,Cz,Cp,beta1,beta2,lambda12,lambda21,mu,NEVOL,f0,delta,alphamax, switching_environments, return_genotypes )

if number_of_realisations>1

[m,alpha]=Evolutionary_trajectories_Many_Realisations(number_of_realisations,m0,alpha0,A,M,T,Cz,Cp,beta1,beta2,lambda12,lambda21,mu,NEVOL,f0,delta,alphamax, switching_environments,  return_genotypes );    

genotypeDATA_m=0;
genotypeDATA_alpha=0;

else

[genotypeDATA_m,genotypeDATA_alpha,m,alpha]=Evolutionary_trajectories_1realisation(m0,alpha0,A,M,T,Cz,Cp,beta1,beta2,lambda21,lambda12,mu,NEVOL,f0,delta,alphamax,return_genotypes);

end

% This is the daddy function.
% fig 3a (m(0),alpha(0))=(0.02,0)
% initialising parameters

A=100; M=1; T=0.025; Cz=0.3; Cp=0; beta1=1; beta2=1; lambda21=0; lambda12=0; mu=0.001; delta=0.01; NEVOL=1500; f0=0.002; m0=0.02; alpha0=0; 
return_genotypes=0; number_of_realisations=25; alphamax=1000;

cd ..

addpath(genpath('simulation_functions'))

cd data_generation_scripts


[~,~,m,alpha]=Evolutionary_trajectories(number_of_realisations,m0,alpha0,A,M,T,Cz,Cp,beta1,beta2,lambda12,lambda21,mu,NEVOL,f0,delta,alphamax,return_genotypes);


cd ..

save('Data_files\Fig_3\panel_a\m_fig3a_002.mat','m');
save('Data_files\Fig_3\panel_a\alpha_fig3a_002.mat','alpha');

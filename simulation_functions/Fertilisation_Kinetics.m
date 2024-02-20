function x = Fertilisation_Kinetics(mgenotypes,alphagenotypes,f,S,A,M,T)

% Simulates the fertilisation kinetics for a fertilisation period of T.
% Input parameters: mgenotypes and alphagenotypes - the mass and fusion rate of all the genotypes.
%                   f - frequency of each genotype.
%                   S - number of distinct genotypes.

if length(mgenotypes)~=S+1 || length(alphagenotypes)~=S+1
    error('length of mgenotypes or alphagenotypes incorrect')
end

if length(f)~=S
    error('length of f incorrect')
end


for i=1:S

N0(i)=( A*M*f(i) )/mgenotypes(i+1);          % N0(i) is the number of unfused cells of genotype i.

end

[t,x] = ode45(@(t,x)unisexual_ODE_Multigenotype(t,x,alphagenotypes,S),[0,T],[N0,zeros(1,nchoosek(S+1 ,2))]);

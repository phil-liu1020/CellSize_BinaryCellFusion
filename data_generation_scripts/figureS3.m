beta1=(0.08:0.08:4);

fusion_env2beta1=zeros(1,length(beta1));
fusion_env1beta1=zeros(1,length(beta1));
fusion_env2betaBH=zeros(1,length(beta1));
fusion_env1betaBH=zeros(1,length(beta1));
figS15output=zeros(1,length(beta1));

addpath(genpath('simulation_functions'))

for i=1:length(beta1)

m0=beta1(i);

if beta1(i)==beta2
    continue
end

[m,alpha,m2,alpha2]=Evolutionary_trajectories_1realisation_plasticity(m0,alpha0,A,M,T,C,beta1(i),beta2,lambda21,lambda12,mu,NEVOL,f0,delta,alphamax);

if (abs(m2(end)-beta2/4)<2*delta && alpha2(end)>4*delta) || alpha2(end)>4*delta
fusion_env2beta1(i)=1;
end

if (abs(m(end)-beta1(i)/4)<2*delta && alpha(end)>4*delta) || alpha(end)>4*delta
fusion_env1beta1(i)=1;
end

        if beta1(i)>beta2
        P1=lambda12/(lambda12+lambda21);
        P2=lambda21/(lambda12+lambda21);    
        else
        P1=lambda21/(lambda12+lambda21);
        P2=lambda12/(lambda12+lambda21);
        end

m0=P1*beta1(i)+P2*beta2;

[m,alpha,m2,alpha2]=Evolutionary_trajectories_1realisation_plasticity(m0,alpha0,A,M,T,C,beta1(i),beta2,lambda21,lambda12,mu,NEVOL,f0,delta,alphamax);

if (abs(m2(end)-beta2/4)<2*delta && alpha2(end)>4*delta) || alpha2(end)>4*delta
fusion_env2betaBH(i)=1;
end

if (abs(m(end)-beta1(i)/4)<2*delta && alpha(end)>4*delta) || alpha(end)>4*delta
fusion_env1betaBH(i)=1;
end




figS15output(i)=fusion_env2beta1(i)+fusion_env1beta1(i)+fusion_env2betaBH(i)+fusion_env1betaBH(i);


end

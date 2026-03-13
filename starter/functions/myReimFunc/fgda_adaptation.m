function [W,ManifoldMappingMutot] = fgda_adaptation(Ctrain,Ytrain,C,ManifoldMappingMutot,METHOD_COV,ARG_COV)
% This was used to adapt the mutotal in an exponential update fashion with
% changing reference
% Ctrain: training set
% Ytrain: Labels
% C: KarcherMean (Barycentre , either learnt through an exponential geodesic update or the one proposed in paper)
% ManifolMappingMutot was just used to track the changes on manifold (to avoid the transportation from tangent space to tangent space)
% Output
% W: Geodesic filters
% A minor tweaking can be done on the expense of computational time : In
% place of euclidean mean of within class covariances, riemannian mean can
% be used


labels = unique(Ytrain);
Nclass = length(labels);
% C=weighted_riemann_mean(Ctrain);

Strain = Tangent_space(Ctrain,C);
Nelec = size(Strain,1);

mu = zeros(Nelec,Nclass);
Covclass = zeros(Nelec,Nelec,Nclass);
for i=1:Nclass
    mu(:,i) = mean(Strain(:,Ytrain==labels(i)),2);
    Covclass(:,:,i) = covariances(Strain(:,Ytrain==labels(i)),METHOD_COV,ARG_COV);
end
  


%% Exponential Update: Adding the new values to the mean
% 
mutot=0.991*Tangent_space(ManifoldMappingMutot,C)+0.009*Strain(:,end);
ManifoldMappingMutot=UnTangent_space(mutot,C);

    
Sb = zeros(Nelec,Nelec);    
for i=1:Nclass
    Sb = Sb+(mu(:,i) - mutot)*(mu(:,i)-mutot)';
end
    
S = mean(Covclass,3);

% S=covariances(Strain,METHOD_COV,ARG_COV);

[W Lambda] = eig(Sb,S);
[~, Index] = sort(diag(Lambda),'descend');

W = W(:,Index);

function [W,mutot] = fgda_adaptationSameCg(trial,C,S,mutot,mu,Nclass,etaValue)
% In this scheme for a fixed reference covariance matrix 'C' we can project
% an incoming trial to the tangent space and just update the total mean 
% Inputs:
%  trial: An incoming Trial
%      C: reference Covariance matrix
%      S: Within class covarinace : (Use the tweaking of riemannian mean on the expense of computational time)
%  mutot: Total mean which is updated in a sequential online exponential
%         window fashion
% NClass: Number of motor imagery classes
%     mu: individual means of all the classes
%etaValue: Learning Parameter
% C=weighted_riemann_mean(Ctrain);

Strain = Tangent_space(trial,C);
Nelec = size(Strain,1);


%% Adding the new values to the mean
% 
mutot=(1-etaValue)*mutot+etaValue*Strain;
Sb = zeros(Nelec,Nelec);    
for i=1:Nclass
    Sb = Sb+(mu(:,i) - mutot)*(mu(:,i)-mutot)';
end
    

[W Lambda] = eig(Sb,S);
[~, Index] = sort(diag(Lambda),'descend');

W = W(:,Index);

function [W,Data] = fgda_adaptation_covariance_too(trial,C,Sb,Data)

%Fluke : Rubbish function: Probably for experimentation


% C=weighted_riemann_mean(Ctrain);

Strain = Tangent_space(trial,C);
Data=[Data,Strain];
Nelec = size(Strain,1);


%% Adding the new values to the mean
% 
% mutot=(1-etaValue)*mutot+etaValue*Strain;
% Sb = zeros(Nelec,Nelec);    
% for i=1:Nclass
%     Sb = Sb+(mu(:,i) - mutot)*(mu(:,i)-mutot)';
% end



S=cov1para(Data');


[W Lambda] = eig(Sb,S);
[~, Index] = sort(diag(Lambda),'descend');

W = W(:,Index);

function [W C] = fgda_play(Ctrain,Ytrain,METHOD_MEAN,ARG_MEAN,METHOD_COV,ARG_COV)

labels = unique(Ytrain);
Nclass = length(labels);

C = mean_covariances(Ctrain,METHOD_MEAN,ARG_MEAN);

% C=weighted_riemann_mean(Ctrain);

Strain = Tangent_space(Ctrain,C);
Nelec = size(Strain,1);

mu = zeros(Nelec,Nclass);
Covclass = zeros(Nelec,Nelec,Nclass);

for i=1:Nclass
    mu(:,i) = mean(Strain(:,Ytrain==labels(i)),2);
    Covclass(:,:,i) = covariances(Strain(:,Ytrain==labels(i)),METHOD_COV,ARG_COV);
%     Covclass(:,:,i) = cov1para(Strain(:,Ytrain==labels(i))');

end
% Cov_total=zeros(Nelec,Nelec);
% for i=1:size(Strain,2)
%     Cov_total=Cov_total+cov1para(Strain(:,i)');
% end
% Cov_total=Cov_total/288;


Cov_total=cov1para(Strain');


% 
% for i=1:size(Ctrain,3)
% %     mu(:,i) = mean(Strain(:,Ytrain==labels(i)),2);
%     Covclass_temp(:,:,i) = cov1para((Strain(:,i)-temp)');
% end
  
% 
% Covclass_temp=cov1para((Strain)');



% mutot = mean(mu,2);
mutot=zeros(Nelec,1);
    
Sb = zeros(Nelec,Nelec);    
for i=1:Nclass
    
%     Sb = Sb+.25*(mu(:,i) - mutot)*(mu(:,i)-mutot)';
    Sb=Sb+cov1para(mu(:,i)');
end
    

% d=rand(1,size(Covclass,2))*1e-5;

% temp_diagonal=diag(d);
% for i=1:size(Covclass,3)
%     Covclass(:,:,i)=Covclass(:,:,i)+temp_diagonal;
% end
% S=riemann_mean(Covclass);

S=mean(Covclass,3);
% S=Stotal-Sb;
% Covclass=covariances((Strain-repmat(mutot,1,size(Ctrain,3))),METHOD_COV,ARG_COV);
% S = Covclass_temp-Sb;
% S=S+Sb;


% S=72*(Cov_total-Sb);

[W Lambda] = eig(Sb,S);
[~, Index] = sort(diag(Lambda),'descend');

W = W(:,Index);

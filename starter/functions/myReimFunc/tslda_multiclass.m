function Ytest = tslda_multiclass(Ctest,Ctrain,Ytrain,varargin)
%% tslda_multiclass(Ctest,Ctrain,Ytrain,varargin)
% Input in this format
    labels = unique(Ytrain);
    Nclass = length(labels);
    
    % Tangent space mapping
    C = riemann_mean(Ctrain);
    Strain = Tangent_space(Ctrain,C);
    Nelec = size(Strain,1);

    % Regularized LDA
    mu = zeros(Nelec,Nclass);
    Covclass = zeros(Nelec,Nelec,Nclass);

    for i=1:Nclass
        mu(:,i) = mean(Strain(:,Ytrain==labels(i)),2);
        Covclass(:,:,i) = covariances(Strain(:,Ytrain==labels(i)),'shcov');
    end

%     mutot = zeros(size(mu(:,i)));

    mutot=mean(mu,2);

    Sb = zeros(Nelec,Nelec);    
    for i=1:Nclass
        Sb = Sb+(mu(:,i) - mutot)*(mu(:,i)-mutot)';
    end

    S = mean(Covclass,3);

    [W Lambda] = eig(Sb,S);
    [~, Index] = sort(diag(Lambda),'descend');
    
    W=W(:,Index);
    
    W = W(:,1:Nclass-1);
    
    % classification
%     if update == 1
%         C = mean_covariances(Ctest,method_mean);
%     end
    Stest = Tangent_space(Ctest,C);
    
    for trialIndex=1:size(Stest,2)
        
        for class=1:Nclass
            temp=W'*(Stest(:,trialIndex)-mu(:,class));
            distance(class)=temp'*temp;
         end
        [~,Ytest(trialIndex)]= min(distance);
    end
    
end
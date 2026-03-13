function [LDAParams] = MulticlassShrinkLDA(Ctrain,Ytrain)
%% Adaptive Regularised lda_multiclass(Ctest,Ctrain,Ytrain,varargin)
% Input in this format
    
    labels = unique(Ytrain);
    Nclass = length(labels);
    Nfeatures=size(Ctrain,1);
    
    
    % Regularized LDA
    mu = zeros(Nfeatures,Nclass);
    Covclass = zeros(Nfeatures,Nfeatures,Nclass);

    for i=1:Nclass
        mu(:,i) = mean(Ctrain(:,Ytrain==labels(i)),2);
        Covclass(:,:,i) = covariances(Ctrain(:,Ytrain==labels(i)),'shcov');
    end

%     mutot = zeros(size(mu(:,i)));

    mutot=mean(mu,2);

    Sb = zeros(Nfeatures,Nfeatures);    
    for i=1:Nclass
        Sb = Sb+(mu(:,i) - mutot)*(mu(:,i)-mutot)';
    end

    S = mean(Covclass,3);

    [W Lambda] = eig(Sb,S);
    [~, Index] = sort(diag(Lambda),'descend');
    
    
    W = W(:,Index);
    LDAParams.W=W(:,1:Nclass-1);
    LDAParams.bias=mu;
    
end



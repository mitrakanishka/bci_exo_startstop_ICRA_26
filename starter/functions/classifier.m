function [est_label, posterior, cm, auc, model, varargout] = classifier(test, test_label, train, train_label, classifier_type, gamma)

idx = [];
switch (classifier_type)
    case 'SVM'
        model = fitcsvm(train, train_label, 'Standardize', false, 'KernelFunction', 'gauss', 'KernelScale', 'Auto', 'Prior', 'uniform');
        [est_label, posterior] = predict(model, test);
        
    case 'LinearSVM'
        model = fitclinear(train, train_label, 'Prior', 'uniform');
        [est_label, posterior] = predict(model, test);
        
    case 'LDA'
        model = fitcdiscr(train, train_label, 'Prior', 'uniform', 'DiscrimType', 'linear');
        [est_label, posterior] = predict(model, test);
    
    case 'diagLDA'
        model = fitcdiscr(train, train_label, 'Prior', 'uniform', 'DiscrimType', 'diaglinear');
        [est_label, posterior] = predict(model, test);
        
    case 'diagQDA'
        model = fitcdiscr(train, train_label, 'Prior', 'uniform', 'DiscrimType', 'diagquadratic');
        [est_label, posterior] = predict(model, test);
        
    case 'RandomForest'
        model = TreeBagger(1000, train, train_label, 'OOBPrediction', 'On', 'Method', 'classification');
        [est_label, posterior] = predict(model, test);
        est_label = str2num(cell2mat(est_label));
        
    case 'SLR_VAR'
        [idx, ~, ~, ~, ~, ~, ~, posterior, ~, est_label, model] = biclsfy_slrvar(train, train_label, test, test_label, 'displaytext', 0, 'nlearn', 300);

    case 'L1_SLR'
        [idx, ~, ~, ~, ~, ~, posterior, ~, est_label, model] = biclsfy_l1slrc(train, train_label, test, test_label, 5, 'displaytext', 0);
end

if not(isempty(idx))
    varargout = {idx};
else 
    varargout = {[]};
end

if not(isempty(test_label))
    cm = confusionmat(est_label, test_label);
    [~, ~, ~, auc] = perfcurve(test_label, posterior(:,2), 1);
else
    cm = [];
    auc = [];
end

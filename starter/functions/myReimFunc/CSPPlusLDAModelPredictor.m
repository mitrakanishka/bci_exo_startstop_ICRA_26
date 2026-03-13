function [LDA_Model,spati] = CSPPlusLDAModelPredictor(training_set,numSpatialFilters)
%% This function is for training an LDAModelgiven Data and its labels

    unique_labels=unique(training_set.y);
    Nclass=size(unique_labels,2);
    spati=OVRAnydataset_onlyTrainingSet(training_set,numSpatialFilters);

    %% LDA Training
    for i=1:size(training_set.x,3)
        temp=spati*training_set.x(:,:,i)';
        var_temp=var(temp');
        feature_mat(:,i)=log(var_temp./sum(var_temp));
    end

    LDA_Model=MulticlassShrinkLDA(feature_mat,training_set.y);  

end


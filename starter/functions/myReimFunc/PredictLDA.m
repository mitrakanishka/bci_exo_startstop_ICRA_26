function [Predicted_label] = PredictLDA(testFeature,Model)
%% Prediction Regularised lda_multiclass(Ctest,Ctrain,Ytrain,varargin)
% Input in this format
    
    Nclass=size(Model.bias,2);
    
    for class=1:Nclass
        distance_to_class=Model.W'*(testFeature-Model.bias(:,class));
        distance(class)=distance_to_class'*distance_to_class;
    end
    [~,Predicted_label]=min(distance);
    
end



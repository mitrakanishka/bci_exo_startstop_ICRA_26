function [meanValueCovMat,varCovMat] = GeoPotatoFromTrainingMat(training_set,labels)

    Nclass=unique(labels);
    
    for i=1:size(training_set,3)
        for class=1:size(Nclass,2)
            cov{i}=training_set(:,:,labels==class);
            [meanValueCovMat(class),varCovMat(class)]=Geo_potato_parameter(cov{i});
        end
        
    end
    
    
end


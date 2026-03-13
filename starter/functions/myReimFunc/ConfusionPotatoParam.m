function [ConfusionData,ConfusionPotatoCenter,ConfusionPotatoMean,ConfusionPotatoVar] = ConfusionPotatoParam(data,Threshold);

    unique_labels=unique(data.labels);
    Nclass=size(unique_labels,2);
    ConfusionData=[];
    
    for i=1:Nclass
        cov{i}=data.data(:,:,data.labels(data.idxTraining)==unique_labels(i));
        C{i} = riemann_mean(data.data(:,:,data.labels(data.idxTraining)==unique_labels(i)));
        [meanValueCovMat(i) varCovMat(i)] = Geo_potato_parameter(cov{i});

        for j=1:size(cov{i},3)
            if(GeoPotatoZScore(distance_riemann(cov{i}(:,:,j),C{i}),meanValueCovMat(i),varCovMat(i))>Threshold(i))
                ConfusionData=cat(3,ConfusionData,cov{i}(:,:,j));
            end
        end
    
    end
    
    ConfusionPotatoCenter=riemann_mean(ConfusionData);
    [ConfusionPotatoMean,ConfusionPotatoVar]=Geo_potato_parameter(ConfusionData);

    



end


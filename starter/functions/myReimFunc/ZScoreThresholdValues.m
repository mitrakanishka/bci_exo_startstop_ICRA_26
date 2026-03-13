function [ThresholdValues] = ZScoreThresholdValues(data,ntop)
    % Normal Potato Parameter

    unique_labels=unique(data.labels);
    Nclass=size(unique_labels,2);

    for i=1:Nclass
        cov{i}=data.data(:,:,data.labels(data.idxTraining)==unique_labels(i));
        C{i} = riemann_mean(data.data(:,:,data.labels(data.idxTraining)==unique_labels(i)));
        [meanValueCovMat(i) varCovMat(i)] = potato_parameter(cov{i});

        for j=1:size(cov{i},3)
            distance(i,j)=distance_riemann(cov{i}(:,:,j),C{i});
        end
        
        for j=1:size(distance(i,:),2)
            zscorevalues(j)=PotatoZScore(distance(i,j),meanValueCovMat(i),varCovMat(i));
        end
        zscorevalues=sort(zscorevalues,'descend');
        
        if ~exist('ntop')
             n_temp_top=OptimalJumpIndex(zscorevalues);
             ThresholdValues(i)=zscorevalues(n_temp_top);
        else
            if(ntop==0)
                ThresholdValues(i)= 10;
            else        
            ThresholdValues(i)=zscorevalues(ntop);
            end
        end
    end

end


function [cov,zscorevaluesUnsorted] = ZscoreThresholdFGMDMSpecific(training_set,training_labels,ntop)
    
    data.data=training_set;
    data.labels=training_labels;
    data.idxTraining=[1:size(training_labels,2)];
    
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
%         zscorevalues=sort(abs(zscorevalues),'descend');
        zscorevaluesUnsorted{i}=zscorevalues;
        zscorevalues=sort(zscorevalues,'descend');
        
        plot(zscorevalues)
%         hold on
        
%         if ~exist('ntop')
%              n_temp_top=OptimalJumpIndex(zscorevalues);
%              ThresholdValues(i)=zscorevalues(i,n_temp_top);
%         else
%             if(ntop==0)
%                 ThresholdValues(i)= 10;
%             else        
%             ThresholdValues(i)=zscorevalues(ntop(i));
%             end
%         end
        clear zscorevalues distance
    end
end


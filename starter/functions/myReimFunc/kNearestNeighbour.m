function [Predicted_Label] = kNearestNeighbour(Trial,Cov,K)

    distance=[];
    for i=1:size(Cov,2)
        temp=Cov{i};
        for j=1:size(Cov{i},3)
            distance=[distance,distance_riemann(Trial,temp(:,:,j))];
        end
    end
    
    labels=[];
    for i=1:size(Cov,2)
        labels=[labels,i*ones(1,size(Cov{i},3))];
    end
    
    unique_labels=unique(labels);
    
    [~,index]=sort(distance,'ascend');
    
    labels=labels(index);
    
    labels=labels([1:K]);
    
    for i=1:size(unique_labels,2)
        count(i)=numel(find(labels==i));
    end

    [maxi_count,Predicted_Label]=max(count);
    

    
    if(numel(find(count==maxi_count))>1)
        Predicted_Label=0;
    end
    
    if(Predicted_Label~=0)
        count_temp=count;
        count_temp(count_temp==maxi_count)=[];
        [maxi_count_temp,~]=max(count_temp);
        if(maxi_count-maxi_count_temp<K/10)
            Predicted_Label=0;
        end
    end
            
        
        
    
    
    
    
end


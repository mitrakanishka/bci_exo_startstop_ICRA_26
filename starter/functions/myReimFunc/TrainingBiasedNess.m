function [acc,Output,conf] = TrainingBiasedNess(data)

unique_labels=unique(data.labels);
Nclass=size(unique_labels,2);

for i=1:Nclass
    cov{i}=data.data(:,:,data.labels(data.idxTraining)==unique_labels(i));
%     C{i} = mean_covariances(data.data(:,:,data.labels(data.idxTraining)==unique_labels(i)),'riemann');
    C{i} = riemann_mean(data.data(:,:,data.labels(data.idxTraining)==unique_labels(i)));
end

for i=1:size(data.idxTraining,2)
    
    trial=data.data(:,:,data.idxTraining(i));
    for class=1:Nclass
        distance(class) = distance_riemann(trial,C{class});
    end
    [~,detected_trial(i)]=min(distance);
        % 
    time_taken(i)=toc;
    
end

trueYtest  = data.labels(data.idxTraining);
acc=100*numel(find(trueYtest-detected_trial==0))/size(data.idxTest,2);

conf=confusionmat(trueYtest,detected_trial);

Output=sum(conf);



    
    
    
    
    











end


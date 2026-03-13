function [EEGsignals] = train_test_generator_POTIOC(subjectID,testingSessionID)
    load(strcat('D:\Dataset\EEGDataCamille\Matlab\subject_[',num2str(subjectID),']','8-30Hz_post.mat'))
    training_set.x=[];
    training_set.y=[];
    for run=1:5
        training_set.x=cat(3,training_set.x,subjectsRuns{1,run}.x);
        training_set.y=cat(2,training_set.y,subjectsRuns{1,run}.y);
    end
    
    training_set.s=subjectsRuns{1,1}.s;
    
    %% testing set
    
    testing_set.x=[];
    testing_set.y=[];
    for run=1:5
        testing_set.x=cat(3,testing_set.x,subjectsRuns{testingSessionID,run}.x);
        testing_set.y=cat(2,testing_set.y,subjectsRuns{testingSessionID,run}.y);
    end
    
    testing_set.s=subjectsRuns{1,1}.s;
    
    
    %% concatenation
    
    EEGsignals.x=cat(3,training_set.x,testing_set.x);
    EEGsignals.y=cat(1,training_set.y',testing_set.y')';
    EEGsignals.s=training_set.s;
    EEGsignals.idxTraining=[1:size(training_set.x,3)];
    EEGsignals.idxTest=[size(training_set.x,3)+1:size(training_set.x,3)+size(testing_set.x,3)];
    


end


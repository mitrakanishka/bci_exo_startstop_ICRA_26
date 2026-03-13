% One V/s rest strategy for feature extraction
% numSpatialFilters: num of spatial filters to be extacted for each one vs
% rest class
% Compatible with any dataset

function spati= OVRAnydataset(EEGsignals,numSpatialFilters)

    labels=unique(EEGsignals.y);
    
    Nclass=size(labels,2);
    training_labels=EEGsignals.y(EEGsignals.idxTraining);
    % Generating the spatial filters
    spati=[];
    for class=1:Nclass
        temp_training.x=EEGsignals.x(:,:,EEGsignals.idxTraining);
        temp_training.y(training_labels==class)='o';
        temp_training.y(training_labels~=class)='r';
        [spatialFilters,~]=learnAutoShrinkCSP(temp_training);
        spati_temp=spatialFilters([1:numSpatialFilters/2,size(spatialFilters,1)-numSpatialFilters/2+1:size(spatialFilters,1)],:);
        spati=[spati;spati_temp];
    end
        

end
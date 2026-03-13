function [ data ] = Cov_barachant( EEGSignals )
% Extraction of covariance matrices in barachant format ( To adjust accordingly with that of covariance toolbox)

nbChannels = size(EEGSignals.x,2);
nbTrials = size(EEGSignals.x,3);
classLabels = unique(EEGSignals.y);
nbClasses = length(classLabels);
% if nbClasses ~= 2
%     disp('ERROR! CSP can only be used for two classes');
%     return;
% end
covMatrices = cell(nbClasses,1); %the covariance matrices for each class

%computing the normalized covariance matrices for each trial
trialCov = zeros(nbChannels,nbChannels,nbTrials);
for t=1:nbTrials
    E = EEGSignals.x(:,:,t)';
    EE = E * E';
    trialCov(:,:,t) = cov1para(EEGSignals.x(:,:,t)./ sqrt(trace(EE))); %estimation of trial covariances using automatic covariance matrix shrinking
%     trialCov(:,:,t) = cov1para(EEGSignals.x(:,:,t)); %estimation of trial covariances using automatic covariance matrix shrinking
%     trialCov(:,:,t)=EE;
    %trialCov(:,:,t) = EE ./ trace(EE);
end

data.data=trialCov;
data.idxTest=EEGSignals.idxTest;
data.idxTraining=EEGSignals.idxTraining;
data.labels=EEGSignals.y;


end


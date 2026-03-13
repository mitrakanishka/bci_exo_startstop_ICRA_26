function [modelParams] = AdaptiveLDAParameterUpdate(modelParams,IncomingTrial,IncomingTrialLabel,learningRate)
    %% only adaptive mean -- Unsupervised
%         modelParams.bias(:,IncomingTrialLabel)=(1-learningRate)*modelParams.bias(:,IncomingTrialLabel)+learningRate*IncomingTrial;
        modelParams.bias(:,IncomingTrialLabel)=euclidean_geodesic(modelParams.bias(:,IncomingTrialLabel),IncomingTrial,learningRate);


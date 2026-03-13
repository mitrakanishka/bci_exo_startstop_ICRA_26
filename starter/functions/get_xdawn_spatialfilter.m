function spatialFilter = get_xdawn_spatialfilter(eeg, labels)

[nSamples, nChannels, nTrials] = size(eeg);
eeg = permute(eeg, [2 1 3]);
eeg = reshape(eeg, [nChannels nSamples*nTrials]);

index.blockLength = nSamples;
index.indexStimulus = 1:nSamples:nSamples*nTrials;

index.indexStimulus(labels == 0) = [];

[enhancedResponse, ~, ~] = mxDAWN(eeg', index, 0);
spatialFilter = enhancedResponse.spatialFilter;

end
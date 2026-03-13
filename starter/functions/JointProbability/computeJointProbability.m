function jp = computeJointProbability(data, params)

[~, nChannels, nTrials] = size(data);
jp = nan(nChannels, nTrials);

for iTrial = 1:nTrials
    for iChannel = 1:nChannels
        jp(iChannel, iTrial) = compute_jp(data(:, iChannel, iTrial), params.data2idx{iChannel}, params.distribution(:, iChannel));
    end
end
end

function jp = compute_jp(data, data2idx, distribution)

idx = data2idx(data);
idx = min(idx, length(distribution));
idx = max(idx, 1);
probability = distribution(idx);
jp = -sum(probability.*log(probability));
end
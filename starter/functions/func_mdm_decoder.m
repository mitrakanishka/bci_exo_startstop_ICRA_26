function [pred_class, prob] = func_mdm_decoder(decoder, eeg_cov)
dist_class_0 = riemann_distance(decoder.prototype_ctrl, eeg_cov);
dist_class_1 = riemann_distance(decoder.prototype_test, eeg_cov);

%% Prediction of the classLabel and probabilities
prob = 1-softmax([dist_class_0^2; dist_class_1^2]);
if (prob(1) >= prob(2))
    pred_class = 0;
else
    pred_class = 1;
end

end

%% helper functions
function a = riemann_distance(A,B)

a = sqrt(sum(log(eig(A,B)).^2));

end


function [Confidence] = ConfidenceCheck(distanceValue,Confidence_Value)
% This fnction check the confidence of the unsupervised clustered trial

[sortedArray,~]=sort(distanceValue,'descend');
Difference=(sortedArray(1)-sortedArray(2))/sortedArray(1);
if(Difference*100>Confidence_Value)
    Confidence=1;
else
    Confidence=0;
end



end


function MCC = computeMCC(confusionMat)

MCC = nan(size(confusionMat,3), 1);
for iMatrix = 1:size(confusionMat, 3)
    TN = confusionMat(1, 1, iMatrix);
    TP = confusionMat(2, 2, iMatrix);
    FP = confusionMat(1, 2, iMatrix);
    FN = confusionMat(2, 1, iMatrix);

    denominator = (TP*TN) - (FP*FN);

    if denominator == 0.0
        MCC(iMatrix) = 0;
    else
        MCC(iMatrix) = denominator/sqrt((TP+FP)*(TP+FN)*(TN+FP)*(TN+FN));
    end
end
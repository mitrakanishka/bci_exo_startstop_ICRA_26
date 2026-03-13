function [Zscore] = PotatoZScore(distance,meanValueCovMat,varCovMat)

    Zscore=((distance-meanValueCovMat)/varCovMat);

end


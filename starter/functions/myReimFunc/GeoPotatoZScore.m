function [Zscore] = GeoPotatoZScore(distance,meanValueCovMat,varCovMat)

    Zscore=log(distance/meanValueCovMat)/log(varCovMat);
%     Zscore=Zscore);
end


function [meanValueCovMat,varCovMat] = Geo_potato_parameter(covMats)
    %Computing the mean covariance matrix
    meanCovMat = riemann_mean(covMats);

    %computing the variance of the covariance matrices
    nbMat = size(covMats,3); %the number of covariance matrices
    meanValueCovMat = 0;
    for t=1:nbMat
        meanValueCovMat = meanValueCovMat + log(distance_riemann(covMats(:,:,t),meanCovMat));
    end
    meanValueCovMat = exp(meanValueCovMat / nbMat);
    
    %% Variance
    
    varCovMat=0;
    
    for t=1:nbMat
        varCovMat = varCovMat + (log(distance_riemann(covMats(:,:,t),meanCovMat)/meanValueCovMat))^2;
    end
    varCovMat = exp(sqrt(varCovMat / nbMat));
end


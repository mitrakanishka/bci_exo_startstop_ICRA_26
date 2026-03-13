function [meanValueCovMat,varCovMat] = potato_parameter(covMats,meanCovMat)
    %Computing the mean covariance matrix
%     meanCovMat = riemann_mean(covMats);

    %computing the variance of the covariance matrices
    nbMat = size(covMats,3); %the number of covariance matrices
    meanValueCovMat = 0;
    for t=1:nbMat
        meanValueCovMat = meanValueCovMat + distance_riemann(covMats(:,:,t),meanCovMat);
    end
    meanValueCovMat = meanValueCovMat / nbMat;
    
    %% Variance
    
    varCovMat=0;
    
    for t=1:nbMat
        varCovMat = varCovMat + (distance_riemann(covMats(:,:,t),meanCovMat)-meanValueCovMat)^2;
    end
    varCovMat = sqrt(varCovMat / nbMat);
end


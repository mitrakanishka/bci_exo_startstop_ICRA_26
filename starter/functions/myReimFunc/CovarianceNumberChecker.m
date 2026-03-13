function [cov] = CovarianceNumberChecker(cov,AverageTrialsPerClass)
    for i=1:size(cov,2)
        temp=cov{i};
        temp=temp(:,:,[size(temp,3)-AverageTrialsPerClass+1:size(temp,3)]);
        cov{i}=temp;
    end
end


function [restDis,classDis,restStab,classStab,betweenClass,classLabels] = PerformanceMetricCovMats(data_rest,data)
%
%Output:
%restDis: a vector with Ncl elements, Ncl being the number of classes in the data. Each element is the distance between the class mean and the rest mean
%classDis: a value indicating the distinctiveness between all classes (whether for 2 or more classes).
%restStab: the stability (inverse of the variance) of the covariance matrices for the rest class (a single value)
%classStab: a vector with each element being the the stability (inverse of the variance) covariance matrices for the corresponding class
%classLabels: the classes labels of the Ncl classes (as indexed in restDis, classDis and classStab)
%RestDisPerBand: (optional) only for multiband EEG data sets and spatioSpectral patterns: provide the list of restDist values for each
%   input band
%classDisPerBand: (optional) only for multiband EEG data sets and spatioSpectral patterns: provide the list of classDist values for each
%   input band
%restStabPerBand: (optional) only for multiband EEG data sets and spatioSpectral patterns: provide the list of rest stab values for each
%   input band
%classStabPerBand: (optional) only for multiband EEG data sets and spatioSpectral patterns: provide the list of class stab values for each
%   input band
%
%by Fabien Lotte (fabien.lotte@inria.fr)
%created: 19/01/2017
%last revised: 04/10/2017



classLabels=unique(data.labels);
nbClasses = length(classLabels);
restDis = zeros(1,nbClasses);

%classDis = zeros(nbClasses);


%Rest
[restMeanCovMat restVarCovMat] = riemann_mean_var(data_rest.data(:,:,data_rest.idxTraining));    
restStab = 1/(1+restVarCovMat);
dimCov = size(restMeanCovMat,1);


% Task
classMeanCovMat = zeros(dimCov,dimCov,nbClasses);
classVarCovMat = zeros(1,nbClasses);

for c=1:nbClasses
    curClassIdx = (data.labels(data.idxTraining) == classLabels(c));
    curClassCovMat = data.data(:,:,curClassIdx);
    [classMeanCovMat(:,:,c),classVarCovMat(c)] = riemann_mean_var(curClassCovMat);
end
classStab = 1./(1+classVarCovMat);


%computing the distances between each class and the rest covariance matrices
for c1=1:nbClasses

    restDis(c1) = distance_riemann(restMeanCovMat,classMeanCovMat(:,:,c1)) / (0.5 * (restVarCovMat + classVarCovMat(c1)));

end

%computing the distinctiveness between all classes
if nbClasses == 2 %2 class case: distance between the two classes mean, divided by their average variance

    classDis = distance_riemann(classMeanCovMat(:,:,1), classMeanCovMat(:,:,2)) / (0.5 * (classVarCovMat(1) + classVarCovMat(2)));
    
elseif nbClasses > 2 %multiclass case: between-class variance divided by the within-class variance (estimated using Riemannian distance)
    %computing the average of the mean covariance matrices for each class (grand average) and the variance around this
    %grand average (i.e., the between class variance)

    [grandMeanCov betweenClassVar] = riemann_mean_var(classMeanCovMat);
    %computing the within class variance, i.e., the mean variance over class
    withinClassVar = 0;
    for c=1:nbClasses
        withinClassVar = withinClassVar + classVarCovMat(c);
    end
    withinClassVar = withinClassVar / nbClasses;
    classDis = betweenClassVar / withinClassVar;

else
    disp('! ERROR! You need at least two classes to compute the class distinctiveness!');
    return;
end

%% PairWiseDistances
betweenClass=[];
for c1=1:nbClasses
    for c2=c1+1:nbClasses
        betweenClass =[betweenClass,distance_riemann(classMeanCovMat(:,:,c2),classMeanCovMat(:,:,c1)) / (0.5 * (classVarCovMat(c2) + classVarCovMat(c1)))];
    end
end
end




function [C] = weighted_riemann_mean(data)
%UNTITLED20 Summary of this function goes here
%   Detailed explanation goes here
convCriteria=10;

while convCriteria>1e-4
    if(convCriteria==10)
        C = riemann_mean(data);
    end
    for j=1:size(data,3)
        weight(j)=distance_riemann(data(:,:,j),C);
    end
    weight=1./weight;
    weight=weight./sum(weight);
    norm_before=norm(C);
    C=karcherBarycenter(data,weight);
    norm_after=norm(C);
    convCriteria=(abs(norm_before-norm_after))/norm_before;
    fprintf('percent change in Norm: %d\n',100*(abs(norm_before-norm_after))/norm_before)


end





end


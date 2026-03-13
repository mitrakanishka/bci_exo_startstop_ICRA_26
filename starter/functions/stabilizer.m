function [transformed_data] = stabilizer(data,reference)
% A function to shift covariance matrices on Riemannian manifold
% For detailed visualisation check Figure 5:
%           Yger F, Sugiyama M., Supervised logeuclidean metric learning for symmetric positive definite matrices. 
%           arXiv preprint arXiv:1502.03505. 2015 Feb 12. 
% Input:
%     data: The covariance matrices of each EEG trial as [Nc*NC*Nt]
%                Nc: Number of Channels Nt: Number of trials
%   refernce: The Affine transformation of size [Nc*Nc]
% Output:
%     transformed_data: The covariance matrices transfored according to affine transform 'reference'
%                      (transformed_data)=reference*data*reference';     
        
    transformed_data=inv(sqrtm(reference))*data*inv(sqrtm(reference))';
        

end

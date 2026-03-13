function [TransformedData] = MapCovariances(data,Center)

% First map the data to identity

    reference=riemann_mean(data);
    
    data_temp=Affine_transformation(data,reference);
% Map to a new Center    
    for i=1:size(data_temp,3)
        TransformedData(:,:,i)=sqrtm(Center)*data_temp(:,:,i)*sqrtm(Center);
    end
end


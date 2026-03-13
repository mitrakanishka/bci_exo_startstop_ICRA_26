function [NearestNeighborMean] = DistanceSorter(chunk,trial,WindowSize)

    for n=1:size(chunk,3)
        dist(n)=distance_riemann(chunk(:,:,n),trial);
    end
    
    [~,index]=sort(dist,'ascend');
    
    chunk=chunk(:,:,index);
    
    NearestNeighborMean=riemann_mean(chunk(:,:,1:WindowSize));
    
    


end


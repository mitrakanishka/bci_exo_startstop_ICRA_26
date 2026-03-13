function [Stability] = StabilityWRTReference(DataChunk,reference)
    
    for i=1:size(DataChunk,3)
        distance(i)=distance_riemann(DataChunk(:,:,i),reference);
    end
    
    Stability=1/(1+distance);


end


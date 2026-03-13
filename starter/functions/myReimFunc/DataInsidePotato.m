function [OutputData] = DataInsidePotato(mean_cov_data,PotatoData,CompleteData)

cov_data_ref=PotatoData;
[meanValue,varianceValue]=potato_parameter(cov_data_ref);
C_temp=[];
for i=1:size(CompleteData,3)
    distance=distance_riemann(CompleteData(:,:,i),mean_cov_data);
    if(PotatoZScore(distance,meanValue,varianceValue)<2)
        C_temp=cat(3,C_temp,CompleteData(:,:,i));
    end
end

OutputData=C_temp;


end


function [outputArg1,outputArg2] = tSNERepersentation(data)



    colors={'r','g','b','k','c'}; % Different Colors maximum 4 diffrent classes : 1. 4 in BCIComp 2. 3 in POTIOC Dataset
    Markers={'o','*'} ;   % 1.Train 2.Test

%     if(Session=='Train')
%         mrk=Markers{1}
%     else
%         mrk=Markers{2}
%     end


    unique_labels=unique(data.labels);
    Nclass=size(unique_labels,2);
    
    VisualiseData=data.data;
    
    for i=1:Nclass
        cov{i}=data.data(:,:,data.labels(data.idxTraining)==unique_labels(i));
        %     C{i} = mean_covariances(data.data(:,:,data.labels(data.idxTraining)==unique_labels(i)),'riemann');
        C{i} = riemann_mean(data.data(:,:,data.labels(data.idxTraining)==unique_labels(i)));
        VisualiseData=cat(3,VisualiseData,C{i});
    end



    
    
    
    T=Tangent_space(VisualiseData,riemann_mean(data.data(:,:,data.idxTraining)));
    tSNE_Mapping=tsne(T');
    
    unique_class=unique(data.labels);
    
    for i=1:size(data.idxTraining,2)
        
        for class =1:Nclass
            if(data.labels(data.idxTraining(i))==unique_class(class))
                plot(tSNE_Mapping(i,1),tSNE_Mapping(i,2),strcat('o',colors{class}))
            end
        end

        hold on
        
        
    end
    
%     i=size(data.idxTest,2)+size(data.idxTraining,2)+1;
%     for class=1:Nclass
%         plot(tSNE_Mapping(i+class-1,1),tSNE_Mapping(i+class-1,2),strcat('p',colors{class}),'MarkerSize',20,'MarkerFaceColor',[1,0.87,0])
%     end
   %% testing
    
        
%     for i=1:size(data.idxTest,2)
%         for class =1:Nclass
%             if(data.labels(data.idxTest(i))==unique_class(class))
%                 plot(tSNE_Mapping(data.idxTest(i),1),tSNE_Mapping(data.idxTest(i),2),strcat('*',colors{class}))
%             end
%         end
% 
%         hold on
%     end

    
    
end






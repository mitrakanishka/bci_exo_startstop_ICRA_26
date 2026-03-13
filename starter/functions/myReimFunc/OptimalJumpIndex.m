function [jumpIndex] = OptimalJumpIndex(Array)
    x=diff(Array);
    
    temp_x=abs(x);
    [~,ind]=max(temp_x);
    
    jumpIndex=ind+1;
end


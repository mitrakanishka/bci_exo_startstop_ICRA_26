function [ acc ] = func_accuracy(conf_matrix)

acc = trace(conf_matrix)/sum(conf_matrix(:));

end
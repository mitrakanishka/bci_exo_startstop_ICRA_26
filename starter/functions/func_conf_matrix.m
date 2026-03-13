function [ conf_matrix ] = func_conf_matrix(y_true, y_pred, class_labels)

if nargin == 2
    class_labels = unique(y_true);
end
n_class = length(class_labels);
conf_matrix = zeros(n_class, n_class);

for i = 1:n_class
    for j = 1:n_class
        y_i = y_true==class_labels(i);
        y_j = y_pred==class_labels(j);
        conf_matrix(i,j) = y_i'*y_j;
    end
end


end
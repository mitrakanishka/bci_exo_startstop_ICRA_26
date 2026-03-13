function [ p ] = func_binomial_chance_level(n, c, alpha)

if c > 1
    p = binoinv(1-alpha, n, 1/c)*100/n;
else
    p = binoinv(1-alpha, n, c)*100/n;
end

end
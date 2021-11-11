function [A,prior] = gen_CF_factors(I,F)

N     = length(I);
A     = cell(N,1);

prior = rand(F, 1); 
prior = prior / sum(prior);

for n = 1 : N
    for f = 1:F
        A{n} = 1*rand(I(n),F) + 1i*1*rand(I(n),F);
    end
    A{n}(ceil(size(A{n},1)/2),:) = ones(1,F);
end

end
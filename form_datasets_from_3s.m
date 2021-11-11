function Y = form_datasets_from_3s(Dataset,marg,I,v)

Y  = cell(size(marg,1),1);

%%% Create the frequency grid G %%%
C = repmat({v},1,3);
D = cell(size(C));
[D{:}] = ndgrid(C{:});

v_combs = [];
v_combs_sub = [];

for i=1:size(D,2)
    v_combs = [v_combs D{i}(:)];
end
   
for i=1:size(marg,1)
    
    Y{i} = zeros(I(marg{i}));  
    Y_vec = mean(exp(1i*v_combs*Dataset(:,marg{i})'),2);
    Y{i} = reshape(Y_vec,I(marg{i}));

end


   
end

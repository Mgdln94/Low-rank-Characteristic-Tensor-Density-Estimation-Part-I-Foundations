function T = N_cpdgen(A,lambda)

if length(A)>1
    
    % mode-1 matrix unfolding in terms of CPD factors   
    T_matricized_debug = (kr(A{end:-1:2}))*(diag(lambda))*((A{1}).');
    % Tensorization
    T_sz = cellfun('size',A(:).',1);
    T = mat2tens(T_matricized_debug.',T_sz,1); %mode_row = 1
   
% %    % Nikos Method
%     T_matricized = A{1}*diag(lambda)*kr(A{end:-1:2}).';
%     T_N = reshape(T_matricized, T_sz);
   
else
    T = A{1}*lambda;
end
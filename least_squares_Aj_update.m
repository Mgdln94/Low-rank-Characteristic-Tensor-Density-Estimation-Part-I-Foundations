function A_j = least_squares_Aj_update(marg,rows,cols,lambda,A,Y_1,n)
       
        F = length(lambda); 
          
        X_bar = [];
        Y_bar = [];
                    
        for i=1:length(rows)
            
                all_idx = marg{rows(i)};
                current_idx = all_idx(cols(i));
                all_idx(cols(i)) = [];
                 
%                 matr_idx = [current_idx sort(all_idx,'ascend')];
%                 kr_idx = sort(all_idx,'descend');
                                
                X_bar = [X_bar ;kr(A{sort(all_idx,'descend')})];
                Y_bar = [Y_bar ;tens2mat(Y_1{i},cols(i)).'];
        
        end

        X_bar = X_bar*diag(lambda);
        
        middle_idx = ceil(size(A{n},1)/2);
        Y_bar(:,middle_idx) = [];
 
        A_j = lsqminnorm(X_bar,Y_bar).';
        A_j = [A_j(1:middle_idx-1,:);ones(1,F);A_j(middle_idx:end,:)];
             
         
end
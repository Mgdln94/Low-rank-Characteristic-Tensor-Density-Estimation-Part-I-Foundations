function [v,F] = K_cross_validation(feature_matrix_X,k,W_sel,F,marg,opts,minVal_vec, maxVal_vec)


  %% Choosing window - K cross-validation
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    training_sz = size(feature_matrix_X,1);
    N = size(feature_matrix_X,2);
    split = ceil(training_sz/k); % Split the dataset into k equal partitions
   
    log_lik = [];
    MSE_proposed_approach= [];
 
    params = combvec(F,W_sel)';
    for s_i=1:size(params,1)
       

        v = 2*pi*(-params(s_i,2):1:params(s_i,2));             
        I = ones(1,N)*(2*params(s_i,2)+1);
        
        [opts.A0,~] = gen_PMF_factors_demo(I,params(s_i,1));
        opts.l0 = 1/params(s_i,1)*ones(params(s_i,1),1);

        log_lik_local = 0;
        MSE_proposed_approach_local  = 0;
        
        for k_indx = 1:k
            
            val_data = feature_matrix_X(1:split,:); % Use first fold as validation data
            train_data = feature_matrix_X(split+1:end,:); % Use the union of other folds as training data and calculate avg-LL
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %% Form Tensor-Triples from the factors %%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            Y = form_datasets_from_3s(train_data,marg,I,v);
            [A_3,l_3,~] = C3_CPD(Y,I,params(s_i,1),opts);
                        
            pdf_val_set = zeros(size(val_data,1),1);
            for s = 1:size(val_data,1)
                pdf_val_set(s) = PDF_point_eval(A_3,l_3,val_data(s,:),params(s_i,2));
            end
        
            % Assemble the samples of validation set and calculate their probability.                         
            log_lik_local = log_lik_local + mean(log(pdf_val_set)); %Average val-set log-likelihood per datapoint
            feature_matrix_X = vertcat(train_data,val_data);
            
            % Predict
            estimated_value = zeros(size(val_data,1),1); 
            for s = 1:size(val_data,1)
                estimated_value(s) = PDF_predict(A_3,l_3,val_data(s,:),params(s_i,2),minVal_vec, maxVal_vec);
            end
            
            val_original_dom =  minVal_vec(N) + val_data(:,N).*(maxVal_vec(N) - minVal_vec(N));
            MSE_proposed_approach_local = MSE_proposed_approach_local +(1/length(val_data))*sum(abs(estimated_value-val_original_dom));
               
        end
        params(s_i,:)
        log_lik = [log_lik log_lik_local/k] % Take the average of these validation avg-LL as the avg-LL of the sample.
        MSE_proposed_approach = [MSE_proposed_approach MSE_proposed_approach_local/k]
    end
    
    
    [cutoff_value,v_cutoff_idx] = max(log_lik); % Choose the window width that minimizes the avg-LL of the sample.
    v = 2*pi*(-params(v_cutoff_idx,2):1:params(v_cutoff_idx,2)); 
    F = params(v_cutoff_idx,1);
    

end



function [norm_data, minVal_vec, maxVal_vec]  = normalize_in_hypercube(feature_matrix_X)

    N = size(feature_matrix_X,2);
    norm_data = feature_matrix_X;
    
    %% Data normalization    
    % Min-Max Feature scaling: Feature scaling is used to bring all values into the range [0,1].
    % This is also called unity-based normalization. This can be generalized to restrict the range of values in the dataset between any arbitrary points a and b. 
    minVal_vec = []; 
    maxVal_vec = [];
    for i = 1:N
        % Min-Max as linear transformation: x' = (b-a)/(tmax-tmin)*x + (a-((b-a)/(tmax-tmin)*tmin)) 
        % Keep the below formulas to go back to original domain
        
        minVal = min(feature_matrix_X(:,i));
        maxVal = max(feature_matrix_X(:,i));
        norm_data(:,i) = (feature_matrix_X(:,i) - minVal) / ( maxVal - minVal );
        
        %% To go back to the original domain
        % your_original_data = minVal + norm_data.*(maxVal - minVal);
        minVal_vec = [minVal_vec minVal];
        maxVal_vec = [maxVal_vec maxVal];
        
    end



end



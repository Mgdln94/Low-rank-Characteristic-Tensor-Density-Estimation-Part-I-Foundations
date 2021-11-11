
function expected_val_original_dom = PDF_predict(A_3,l_3,data_point,K,minVal_vec, maxVal_vec)
   
            F = size(A_3{1},2);
            N = size(A_3,1);
            
            v = 2*pi*(-K:1:K); 
            
            %% Calculating c_1
            p1 = ones(1,F);
            for n = 1:N-1
                p1 = p1.*(abs(real((A_3{n}).'*exp(-1i*v'*data_point(n)))))';
            end
            c1 = p1*l_3;
            
            %% Calculating c_2
            alpha = -1i*v';
            c_2 = exp(alpha)./alpha - exp(alpha)./alpha.^2 + 1./alpha.^2;
            c_2(isinf(c_2)|isnan(c_2)) = 1/2; % Replace NaNs and infinite values with  1/2
            
          
            %% Calculating expectation
            p2 = ones(1,F);
            for n = 1:N-1
                p2 = p2.*(abs(real((A_3{n}).'*exp(-1i*v'*data_point(n)))))';
            end
            p2 = p2.*(abs(real((A_3{N}).'*c_2)))';
            expected_val = p2*l_3;
            

            expected_val_original_dom = minVal_vec(N) + (expected_val/c1).*(maxVal_vec(N) - minVal_vec(N));
            
        
   
end
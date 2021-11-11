function PDF_point_eval = PDF_point_eval(k_vec,A_3,l_3,data_point)

            F = size(A_3{1},2);
            N = size(A_3{1},1);
            
            p1 = ones(1,F);
            for n = 1:N
                p1 = p1.*A_3{n}(k_vec(n),:);
            end
            
            CF_point_eval = p1*l_3*exp(-1i*2*pi*data_point(n));
   
end
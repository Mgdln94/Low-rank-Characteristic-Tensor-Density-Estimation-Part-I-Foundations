
function PDF_point_eval = PDF_point_eval(A_3,l_3,data_point,K)
   
            F = size(A_3{1},2);
            N = size(A_3,1);
            
            v = 2*pi*(-K:1:K); 
            p1 = ones(1,F);
            for n = 1:N
                p1 = p1.*(abs(real((A_3{n}).'*exp(-1i*v'*data_point(n)))))';
            end
            
            PDF_point_eval = p1*l_3;
   
end
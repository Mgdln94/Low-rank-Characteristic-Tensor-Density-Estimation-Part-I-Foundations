function [lambda,U] = admm_lambda_update(G,V,lambda,U,n,opts,rho,max_iter)
[ ~, k ] = size(lambda);
    

%     flag = 1;
%     rho = rho/10;
%     while flag ~= 0 && rho < 10
%          rho = 10*rho;
         [L,flag] = chol(G + rho*eye(k), 'lower');                 
%     end
    if flag ~= 0
        tol = 1e-3;
        for itr = 1:max_iter
            lambda_0 = lambda;
            lambda_t = L'\ ( L\ ( V + rho*(lambda+U)') );
            lambda  = proxr(real(lambda_t'-U), opts, n); % Fix this U is complex and thus the next lambda update becomes complex
            U  = U + lambda - lambda_t';
            r  = lambda - lambda_t';
            s  = (lambda - lambda_0);
            if  norm(r(:)) < tol  && norm(s(:)) < tol
                break
            end
        end
    end
end

function lambda = proxr(Ab,opts,n)
switch opts.constraint{n}
    case 'simplex'
        lambda = reshape(ProjectOntoSimplex(Ab(:),1),size(Ab));         
end
end
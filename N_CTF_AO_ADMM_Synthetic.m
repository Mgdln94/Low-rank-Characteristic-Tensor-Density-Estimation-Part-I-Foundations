function [A,lambda,Out] = N_CTF_AO_ADMM_Synthetic(Y,I,F,opts)
% input:
%       y   : measurement vector
%       I   : number of rows of each factor
%       F   : rank
%       opts.
%           loss        = 'ls'
%           constraint  = 'nonnegative'
%                       = 'simplex_col'
%                       = 'simplex'
%
%           A0          : initial tensor factors, default: random
%           max_iter    : max number of iterations, default: 1000
%           tol_impr    : min cost improvement to terminate the algorithm, default: 1e-7
%           X_true      : N-way tensor (optional)
%           A_true      : true factors (optional)
% Output:
%       A   : cell array containing the factors
%       Out.
%           iter            : number of iterations
%           hist_cost       : history of cost
%           hist_rel_cost   : history of relative cost
%           hist_ten_cost   : history of tensor relative cost
%           time_instants   : time at the end of each iteration

%% Parameters
N  = length(I);
U  = cell(N+1,1); % scaled dual variables
GG = cell(N,1);

if ~isfield(opts,'loss'),       opts.loss       = 'ls';          end
if ~isfield(opts,'constraint'), opts.constraint = 'nonnegative'; end
if ~isfield(opts,'max_iter'),   opts.max_iter   = 2000;          end
if ~isfield(opts,'tol_impr');   opts.tol_impr   = 1e-7;          end

rho    = opts.rho;
A      = opts.A0;          % initial tensor factors
lambda = opts.lambda0;     % Prior of hidden variable

for n = 1:N
    GG{n} = A{n}'*A{n};    % A^T*A cache
end

for n = 1:N
    U{n} = zeros(size(A{n}));
end
U{N+1}   = zeros(size(lambda'));

% cost and relative cost
Out.hist_cost         = zeros(opts.max_iter,1);
Out.hist_rel_cost     = zeros(opts.max_iter,1);

% tensor and factors relative cost
Out.hist_tensor_cost  = zeros(opts.max_iter,1);
Out.hist_factor_cost  = cell(N,1);

Out.time_instants     = zeros(opts.max_iter,1);

for n = 1:N
    Out.hist_fac_cost{n}   = zeros(opts.max_iter,1);
end
opts.F = F;

%%% Iterations of block-coordinate update
iter = 1;

%%% Matricize each smaller tensor and concatenate
Ym = matricize(Y,N,opts);

%%% Precompute to save time
prod  = cell(N,1);
for n = 1:N
    [row_n,~] = find(opts.marg == n);
    row_n     = sort(row_n);
    for r = 1:length(row_n)
        prod{n} = [prod{n}; setdiff(opts.marg(row_n(r),:),n)];
    end
end

tic
while(1)
    %%%  Solve each subproblem with ADMM
    %%%  marginal index is increasing
    for n = 1:N
        V = zeros(F,F);
        %%% Compute KhatriRao products
        W     = zeros(size(Ym{n},1),F);
        ind_n = 0;
        for i = 1 : size(prod{n},1)
            G = ones(F,F);
            for j = prod{n}(i,:)
                G = G .* GG{j};
            end
            V       = V + G;
            tmp     = khatrirao(A{prod{n}(i,:)},'r');
            tmp_sz  = size(tmp,1);
            W(ind_n+1:ind_n+tmp_sz,:) = tmp;
            ind_n   = ind_n + tmp_sz;
        end
        V = diag(lambda)*V*diag(lambda);
        W = W*diag(lambda);
        max_iter = 1000;
        [A{n}, U{n}] = N_CTF_AO_ADMM_sub(Ym{n},W,V,A{n},U{n},n,opts,rho(1),max_iter);
        
%         cvx_begin quiet
%         variable A_cvx(I(n),F)
%         minimize(sum(sum_square(Ym{n}-W*A_cvx')))
%         subject to
%         A_cvx>=0;
%         ones(1,I(n))*A_cvx == ones(1,F);
%         cvx_end
%         A{n} = A_cvx;
        
        GG{n} = A{n}'*A{n};
    end
        
    V = zeros(F,F);
    W = zeros(size(Ym{N+1},1),F);
    ind_n = 0;
    
    for i = 1 : size(opts.marg,1)
        prodd = opts.marg(i,:);
        G = ones(F,F);
        for j = prodd
            G = G .* GG{j};
        end
        V       = V + G;
        tmp     = khatrirao(A{prodd}, 'r');
        tmp_sz  = size(tmp,1);
        W(ind_n+1:ind_n+tmp_sz,:) = tmp;
        ind_n   = ind_n + tmp_sz;
    end
    max_iter = 1000;
    [lambda, U{N+1}] = N_CTF_AO_ADMM_sub(Ym{N+1},W,V,lambda',U{N+1},N+1,opts,rho(2),max_iter);
    lambda = lambda';
    
%     cvx_begin quiet
%     variable A_cvx(F,1)
%     minimize sum_square(Ym{n+1}-W*A_cvx)
%     subject to
%     A_cvx>=0;
%     ones(1,F)*A_cvx == 1;
%     cvx_end
%     lambda = A_cvx;

    [Out.hist_cost(iter),Out.hist_rel_cost(iter)]  = Loss_Coupled(Y, A, opts,lambda);
    Out.time_instants(iter) = toc;
    
    if iter>1
        if (iter == opts.max_iter ||  abs(Out.hist_rel_cost(iter) - Out.hist_rel_cost(iter-1)) < opts.tol_impr )
            Out.iter = iter;
            Out.time_instants(iter+1:end) = [];
            Out.hist_rel_cost(iter+1:end) = [];
            Out.hist_cost(iter+1:end)     = [];
            break;
        end
    end
    if mod(iter,10) == 0, fprintf('Iteration : %d rel cost : %d rel cost diff : %d \n', iter, Out.hist_rel_cost(iter), abs(Out.hist_rel_cost(iter) - Out.hist_rel_cost(iter-1))); end;
    iter = iter + 1;
end
end

function Ym = matricize( Y, N, opts )
Ym  = cell(N+1,1);
r   = size(opts.marg,1);
for i = 1:r
    Y_ = tensor(Y{i});
    for d = 1:ndims(Y_)
        temp  = tenmat(Y_,d);
        Ym{opts.marg(i,d)} = [Ym{opts.marg(i,d)}; temp.data'];
    end
    Ym{N+1} = [Ym{N+1}; Y_(:)];
end
end

function [err,rel_err] = Loss_Coupled(Y,A,opts,lambda)
rel_cost = zeros(size(opts.marg,1),1);
c = zeros(size(opts.marg,1),1);
nrm     = 0;
err     = 0;
for i = 1 : size(opts.marg,1)
    Y_ = tenmat(tensor(Y{i}),1)';
    c(i) = norm(Y_.data - khatrirao(A{opts.marg(i,[2:end])} , 'r') * diag(lambda) * A{opts.marg(i,1)}' ,'fro');
    rel_cost(i) = c(i) / norm(Y{i}(:));
    err     = err + c(i)^2;
    nrm = nrm + norm(Y{i}(:))^2;
end
rel_err = err/nrm;
err     = 1/2 * err;

end
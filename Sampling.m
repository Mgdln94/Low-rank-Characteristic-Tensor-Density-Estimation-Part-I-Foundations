% clc
% close all
% clear all


rng(0) % For reproducibility

folder = fileparts(which(mfilename));
addpath(genpath(folder)); % Add current folder plus all subfolders to the path.

% https://web.stanford.edu/~hastie/StatLearnSparsity_files/DATA/zipcode.html
data1 = csvread('zip_train.csv');
data2 = csvread('zip_test.csv');
data = [data1;data2];

first_Column = data(:,1);
data = data((first_Column == 0),:); % Change first_Column == 1 if you want to model 1s or choose any other number to model
data = data(:,2:end);

M = size(data,1);
N = size(data,2);

test_sample_size =  floor(0.2*M);
W_sel = 15;
F = 8;
num_samples = 20;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Data normalization %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

a = 0;
b = 1;
[data_norm, a_bar, b_bar] = normalize_data_in_hypercube(data,a,b);


stp = 0.5;
t = a:stp:b;

marg_third = [];
marg_third = combnk(1:N,3);
marg_third = marg_third(randperm(size(marg_third, 1), 5000), :); % If you need fast results, you can replace 5000 with 2000 or something smaller

for i=1:2:254
    marg_third = [marg_third;[i i+1 i+2]];
end
marg_third = [marg_third;[254 255 256]];

marg_third = unique(marg_third,'rows');

marg = num2cell(marg_third,2);
opts.marg   = marg;

sims = 1;
log_lik_proposed_approach = zeros(1,sims);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



for sim = 1:sims
    
    sim
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Parameters %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    for n = 1:N
        opts.constraint{n} = 'none';
    end
    opts.constraint{N+1}   = 'simplex';
    
    opts.max_iter = 35;
    opts.rho      = 0.01;
    opts.tol_impr = 1e-3;
    
    train_data = data_norm;
    v = pi*(-W_sel:1:W_sel);
    smooth_param = (length(v)-1)/2;
    I = ones(1,N)*length(v);
    
    [opts.A0,~] = gen_CF_factors(I,F);
    opts.l0 = 1/F*ones(F,1);
    
    %Y = form_datasets_from_3s(train_data,marg,I,v);
    [A_3,l_3,~] = C3_CPD(Y,I,F,opts);
    pdfs_coupled = cell(N,1);
    for i = 1:N
        pdfs_coupled{i} = (abs(real((A_3{i}).'*exp(-1i*v'*t))))';
    end
    
    figure()
    %% Sample from the distribution %%  
    rand_lamda_index = randsample(1:F,num_samples,true,l_3);
    samples = ones(num_samples,N);
    for i = 1:num_samples
        for n = 1:N
            samples(i,n) = randsample(t, 1, true, pdfs_coupled{n}(:,rand_lamda_index(i)));
        end
        subplot(2,num_samples,i), imshow(mat2gray(vec2mat(samples(i,:),16)))
        subplot(2,num_samples,i+20), imshow(mat2gray(vec2mat(data_norm(i,:),16)));
    end
    
     
end


function [feature_matrix_X, a_bar, b_bar]  = normalize_data_in_hypercube(feature_matrix_X,a,b)

    N = size(feature_matrix_X,2);
    
    %% Data normalization    
    % Min-Max Feature scaling: Feature scaling is used to bring all values into the range [0,1].
    % This is also called unity-based normalization. This can be generalized to restrict the range of values in the dataset between any arbitrary points a and b. 
    a_bar = []; 
    b_bar = [];
    for i = 1:N
        % Min-Max as linear transformation: x' = (b-a)/(tmax-tmin)*x + (a-((b-a)/(tmax-tmin)*tmin)) 
        % Keep the below formulas to go back to original domain
        
        tmin = min(feature_matrix_X(:,i));
        tmax = max(feature_matrix_X(:,i));
        a_bar_local = (b-a)/(tmax-tmin);
        b_bar_local = a-(a_bar_local*tmin);
        a_bar = [a_bar a_bar_local];
        b_bar = [b_bar b_bar_local];
        
        % If A is a matrix, table, or timetable, then normalize operates on each column of data separately.
        feature_matrix_X(:,i) = normalize(feature_matrix_X(:,i), 'range', [a b]);
    end



end
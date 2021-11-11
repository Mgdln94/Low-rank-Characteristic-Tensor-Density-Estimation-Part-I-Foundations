clc
close all
clear all
 
rng(0) % For reproducibility

folder = fileparts(which(mfilename)); 
addpath(genpath(folder)); % Add current folder plus all subfolders to the path.


data = csvread('winequality_red.csv');
% data = csvread('Winequality_white.csv');
% data = csvread('FOTP.csv');
% data = swap_ends(data);
% data = csvread('PCB.csv');
% data = swap_ends(data);
% data = importdata('Superconduct.data');
% data = csvread('Gas_Sensor.csv');

% Tunable Parameters - Smoothing parameter/ Tensor rank
W_sel = [10 15 20];
F_vec = [10 20 30 50];

M = size(data,1);
N = size(data,2);

test_sample_size =  floor(0.2*M);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Data normalization %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

a = 0;
b = 1;
[data_norm, minVal_vec, maxVal_vec] = normalize_in_hypercube(data);
%% Going back to the original domain
% your_original_data = [];
%     for i = 1:N
%         your_original_data(:,i) = minVal_vec(i) + data_norm(:,i).*(maxVal_vec(i) - minVal_vec(i));
%     end

%csvwrite('norm_data_wine_red.txt',data_norm)
%csvwrite('norm_data_wine_white.txt',data_norm)
%csvwrite('norm_FOTP.txt',data_norm)

stp = 0.01;
t = a:stp:b;
      
marg_third = combnk(1:N,3);
marg = num2cell(marg_third,2);
opts.marg   = marg;

marg_third = [];
marg_third = combnk(1:N,3);

% Uncomment the following lines for higher dimensional datasets (other the red/white wine datasets)
% marg_third = marg_third(randperm(size(marg_third, 1), 1000), :);
% for i=1:2:N-2
%     marg_third = [marg_third;[i i+1 i+2]];
% end

marg_third = unique(marg_third,'rows');
marg = num2cell(marg_third,2);
opts.marg   = marg;

sims = 5;
log_lik_proposed_approach = zeros(1,sims);
MSE_proposed_approach = zeros(1,sims);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
k = 5; % Number of partitions 
     
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
        opts.rho  = 1e-2;
        opts.tol_impr = 1e-3;
        
        % Random Shuffling
        data_norm = data_norm(randperm(size(data_norm, 1)), :);
        %%% Training/ Testing split
        test_data = data_norm(1:test_sample_size,:);
        train_data = data_norm(test_sample_size+1:M,:);
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% Choosing the number of coefficients - K-fold cross validation %%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
        [v,F] = K_cross_validation(train_data,k,W_sel,F_vec,marg,opts,minVal_vec, maxVal_vec);

      
        smooth_param = (length(v)-1)/2;
   
        I = ones(1,N)*length(v);
                            
        [opts.A0,~] = gen_CF_factors(I,F);
        opts.l0 = 1/F*ones(F,1);
              
        Y = form_datasets_from_3s(train_data,marg,I,v);       
        [A_3,l_3,~] = C3_CPD(Y,I,F,opts);
        
        %% Calculate the density of samples  %%
        estimated_density_samples = zeros(size(test_data,1),1);           
        for s = 1:size(test_data,1)           
            estimated_density_samples(s) = PDF_point_eval(A_3,l_3,test_data(s,:),smooth_param);          
        end       
        log_lik_proposed_approach(sim) = (1/length(test_data))*sum(log(estimated_density_samples));
        
        %% Calculate the MAE of samples  %%
        estimated_value = zeros(size(test_data,1),1); 
        
        for s = 1:size(test_data,1)           
            estimated_value(s) = PDF_predict(A_3,l_3,test_data(s,:),smooth_param,minVal_vec, maxVal_vec);          
        end
        
        test_original_dom =  minVal_vec(N) + test_data(:,N).*(maxVal_vec(N) - minVal_vec(N));
        MSE_proposed_approach(sim) = (1/length(test_data))*sum(abs(estimated_value-test_original_dom)); 
        
    end
%end

 function B = swap_ends(A)     % function definition
A(:,[1 end])=A(:,[end 1]);
B=A;
end

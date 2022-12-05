% The script is designed to evaluate the performance of the subspace
% merging procedure.
%
% Based on work of Grammenos et al.: https://arxiv.org/abs/1907.08059
%
% Author: Andreas Grammenos (ag926@cl.cam.ac.uk)
%
% Last touched date: 02/06/2020
% 
% License: GPLv3
%

%% Initliasation
clc; clear; close all;

% for reproducibility
rng(300);

% the type used
params.type = "merge";
% enable printing
params.pflag = 0;

params = setup_vars(params);

% put everything in one plot
one_plot = 1;

fprintf("\n -- Merging Test suite starting\n\n");

%% Test execution (Basic)

fprintf("\n >> Running basic merging tests...\n");

% configuration for basic test
feats = 500;  % number of features
alpha1 = 0.1; % alpha for first distribution
alpha2 = 1;   % alpha for second distribution
T1 = 600;     % columns for first distribution
T2 = 400;     % columns for second distribution
T3 = 300;     % columns for third distribution

r = feats;  % common rank
r1 = 10;    % rank 1 test
r2 = 10;    % rank 2 test
r3 = 5;     % rank 3 test

% use basic merge to test
% 1 - naive svd merge
% 2 - naive qr merge
% 3 - block matrix merge using econ. svd and qr
algo_type = 1; 
lambda1 = 1;
lambda2 = 1;

% generate the synthetic dataset
Y1 = rand(feats, T1);
Y2 = rand(feats, T2);
Y3 = rand(feats, T3);

% perform the svds with same rank (at first)

yy = [Y1, Y3, Y2];
[Uff, Sff, ~] = svds(yy, r);

% try to merge


% diff
fprintf(" ** Using equal ranks (r: %d)\n", r);
fprintf(" ** Subspace (abs) diff: %d\n", norm(abs(Uf)-abs(Uff), 'fro'));
fprintf(" ** Subspace diff: %d\n", norm(Uf-Uff, 'fro'));
fprintf(" ** Singular value diff: %d\n", norm(Sf-Sff(1:r, 1:r), 'fro'));

% perform the svds with same rank (at first)
[U1, S1, ~] = svds(Y1, r1);
[U2, S2, ~] = svds(Y2, r2);
[Uff, Sff, ~] = svds([Y1, Y2], max(r1, r2));

% try to merge


% diff
fprintf("\n ** Using unequal ranks (r1: %d, r2: %d)\n", r1, r2);
fprintf(" ** Subspace (abs) diff: %d\n", norm(abs(Uf)-abs(Uff), 'fro'));
fprintf(" ** Subspace diff: %d\n", norm(Uf-Uff, 'fro'));
fprintf(" ** Singular value diff: %d\n", norm(Sf-Sff, 'fro'));


fprintf("\n >> Finished basic merging tests...\n");


%% Running over variable sizes to evaluate error scaling against SVD

fprintf("\n >> Running over variable sizes...\n\n");

% number of features (ambient dimension)
feats = 800;
% number of vectors (columns)
T = [feats, 2*feats, 3*feats, 4*feats, 5*feats];
% T = [200, 400, 600, 800, 1000];
% target rank
r = 100;

% synthetic dataset parameter for Power Law
synth_params.spectrum_type = "pl";
synth_params.alpha = 1;
synth_params.lambda = .01;

% preallocation of error arrays
errf_fast_u = zeros(1, size(T, 2));
errf_fast_g = zeros(1, size(T, 2));
errf_svd_u = zeros(1, size(T, 2));
errf_svd_g = zeros(1, size(T, 2));

errf_svd_u1 = zeros(1, size(T, 2));
errf_svd_g1 = zeros(1, size(T, 2));
errf_svd_v1 = zeros(1, size(T, 2));

errf_svd_u2 = zeros(1, size(T, 2));
errf_svd_g2 = zeros(1, size(T, 2));
errf_svd_v2 = zeros(1, size(T, 2));


f_times = zeros(1, size(T, 2));
s_times = zeros(1, size(T, 2));


errf_svd_L = zeros(1, size(T, 2));
errf_svd_S = zeros(1, size(T, 2));
% run for T
for i = 1:size(T, 2)
  fprintf("\n == Running for T: %d\n", T(i));
  % define the chunk size for this particular instance
  chunkSize = T(i)/2;
  % generate the data
  Y = synthetic_data_gen(feats, T(i), synth_params);
  % perform the offline r-SVD on the full dataset
  % True value of RPCA
  [L_u,L_sigma,L_v,S_u,S_sigma,S_v] = RobustPCA(Y); 

  % use halves
  
  % first half
  [L_u1,L_sigma1,L_v1,S_u1,S_sigma1,S_v1] = RobustPCA(Y(:, 1:chunkSize));

  % second half
  [L_u2,L_sigma2,L_v2,S_u2,S_sigma2,S_v2] = RobustPCA(Y(:, chunkSize+1:end));
  % svd merge
  s_tic = tic;
  
  [Ul_f_svd, Sl_f_svd,Vl_f_svd] = frpca_subspace_merge2(L_u1, L_sigma1, L_u2, L_sigma2);
  [Us_f_svd, Ss_f_svd,Vs_f_svd] = frpca_subspace_merge2(S_u1, S_sigma1, S_u2, S_sigma2);
  s_times(i) = toc(s_tic);
  
  % fast merge
  f_tic = tic;
  f_times(i) = toc(f_tic);
  fprintf('\n -- One epoch...');
  % check the errors using fast method


  old_L = L_u*L_sigma*L_v;
  old_S = S_u*S_sigma*S_v;
  new_L = Ul_f_svd * Sl_f_svd * Vl_f_svd;
  new_S = Us_f_svd * Ss_f_svd * Vs_f_svd;
  errf_svd_L(i) = (1/T(i)) * immse(old_L, new_L);
  errf_svd_S(i) = (1/T(i)) * immse(old_S, new_S);

end
%% 

fprintf("\n >> Finishd running over variable sizes...");
fprintf("\n >> Plotting results.");

my_ticks = size(T, 2);

plottools("on");
figure;
subplot(1, 3, 1)
plot(1:my_ticks, errf_svd_L, '*-', 'LineWidth', 2);
hold on;
plot(1:my_ticks, errf_svd_S, '+-', 'LineWidth', 2);
hold off;
title("L AND S ERROR RPCA VS FRPCA");
legend("fast", "svd");
xticks(1:my_ticks);
xticklabels(num2cell(T));
xlabel("T");
ylabel("error (mse)");



% finally set the fonts to be larger
set(findall(gcf,'-property','FontSize'),'FontSize',14)
% make the figure larger from the get go
set(gcf, 'Units', 'Normalized', 'Position',  [.4, .1, .3, .6])

fprintf("\n -- Merging Test suite finished\n");



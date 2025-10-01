%% Geometric Spatio-Spectral Total Variation for Hyperspectral Image Denoising and Destriping
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: Shingo Takemoto (takemoto.s.e908@m.isct.ac.jp)
% Last version: Oct. 1, 2025
% Article: S. Takemoto, S. Ono, 
%   ``Geometric Spatio-Spectral Total Variation for Hyperspectral Image Denoising and Destriping''
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear
close all
addpath(genpath('sub_functions'))
fprintf('******* initium *******\n');
rng('default')

%% Generating observation
%%%%%%%%%%%%%%%%%%%%% User settings of experiment %%%%%%%%%%%%%%%%%%%%%%%%%%%%
deg.Gaussian_sigma      = 0.1; % Standard derivation of Gaussian noise
deg.sparse_rate         = 0.05; % Rate of sparse noise
deg.stripe_rate         = 0.05; % Rate of stripe noise
deg.stripe_intensity    = 0.5; % Range of intensity for stripe noise
deg.deadline_rate       = 0.01; % Rate of deadline noise

% image = 'JasperRidge';
image = 'PaviaUniversity';

noise_seed = 'default';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[HSI_clean] = Load_HSI(image);
[HSI_noisy] = Generate_obsv(HSI_clean, deg, noise_seed);

HSI_clean = single(HSI_clean);
HSI_noisy = single(HSI_noisy);

[n1, n2, n3] = size(HSI_noisy);


%% Setting parameters
%%%%%%%%%%%%%%%%%%%%% User Settings %%%%%%%%%%%%%%%%%%%%%%%%%%%%
params.rho = 0.95; % Parameter for the radii of the noise constraints

params.omega = 0.03; % Balancing parameter between first and second-order differences

params.stopcri = 1e-5; % Stopping criterion
params.maxiter = 20000; % Maximum number of iterations

use_GPU = 1; % 1, if you use GPU, 0, if you use CPU
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Radius of v-centered l2-ball constraint serving data-fidelity
params.epsilon = params.rho * deg.Gaussian_sigma * sqrt(n1*n2*n3 * (1 - deg.sparse_rate));

% Radius of l1-ball constraint for sparse noise
params.alpha = params.rho * (0.5 * n1*n2*n3 * deg.sparse_rate);

% Radius of l1-ball constraint for stripe noise
params.beta = params.rho * n1*n2*n3 * (1 - deg.sparse_rate) ...
    * deg.stripe_rate * deg.stripe_intensity / 2; 


%% Showing settings
fprintf('~~~ SETTINGS ~~~\n');
fprintf('Image: %s Size: (%d, %d, %d)\n', image, n1, n2, n3);
fprintf('Gaussian sigma: %g\n', deg.Gaussian_sigma);
fprintf('Sparse rate: %g\n', deg.sparse_rate);
fprintf('Stripe rate: %g\n', deg.stripe_rate);
fprintf('Stripe intensity: %g\n', deg.stripe_intensity);
fprintf('Rho: %g\n', params.rho);
fprintf('Omega: %g\n', params.omega);
fprintf('Stopping criterion: %g\n', params.stopcri);


%% Denoising and destriping
if use_GPU == 1
    [HSI_restored, ~, ~, ~] = GeoSSTV_GPU(HSI_noisy, params); % for GPU user
elseif use_GPU == 0
    [HSI_restored, ~, ~, ~] = GeoSSTV_CPU(HSI_noisy, params); % for CPU user
else
end


%% Plotting results
val_mpsnr  = calc_MPSNR(HSI_restored, HSI_clean);
val_mssim  = calc_MSSIM(HSI_restored, HSI_clean);

fprintf('~~~ RESULTS ~~~\n');
fprintf('MPSNR: %#.4g\n', val_mpsnr);
fprintf('MSSIM: %#.4g\n', val_mssim);


%% Showing result images
visband = 59; % band to visualize

figure(1);
subplot(1,3,1), imshow(HSI_clean(:,:,visband)), title('Ground-truth');
subplot(1,3,2), imshow(HSI_noisy(:,:,visband)), title('Noisy HSI');
subplot(1,3,3), imshow(HSI_restored(:,:,visband)), title('Restored HSI');

fprintf('******* finis *******\n');

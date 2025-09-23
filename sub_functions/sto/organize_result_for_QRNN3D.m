clear
close all;

addpath(genpath('sub_functions'))
addpath('func_metrics')


%% Selecting conditions
noise_conditions = { ...
    %g      ps     pt     t_int pd
    {0.1,   0,     0,     0,    0   }, ...  
    {0.1,   0.05,  0,     0,    0   }, ... 
    {0.1,   0,     0.05,  0.5,  0   }, ... 
    {0.1,   0,     0,     0,    0.01}, ... 
    {0.1,   0.05,  0.05,  0.5,  0.01}, ... 
};

% idc_noise_conditions = 1:size(noise_conditions, 2);
idc_noise_conditions = 1:5;


images = {...
    "JasperRidge", ...
    "PaviaU", ...
};

idc_images = 1:numel(images);


ckpts = {
    "paviaft", ...
    "complex", ...
};

idc_ckpts = 1:numel(ckpts);


name_method = "QRNN3D";

dir_result_root = "H:\マイドライブ\MATLAB_Share\Data_QRNN3D\results";

load("dir_save_folder.mat", "dir_save_folder");
load("dir_save_comp_folder.mat", "dir_save_comp_folder")


idx_data = 1;
total_data = length(idc_noise_conditions) * length(idc_images) * length(ckpts);

for idx_noise_condition = idc_noise_conditions
for idx_image = idc_images
for idx_ckpt = idc_ckpts
%% Determining conditions
deg.gaussian_sigma      = noise_conditions{idx_noise_condition}{1};
deg.sparse_rate         = noise_conditions{idx_noise_condition}{2};
deg.stripe_rate         = noise_conditions{idx_noise_condition}{3};
deg.stripe_intensity    = noise_conditions{idx_noise_condition}{4};
deg.deadline_rate       = noise_conditions{idx_noise_condition}{5};
image = images{idx_image};
ckpt = ckpts{idx_ckpt};

params.ckpt = ckpt;
name_params_savetext = ckpt;


%% Loading results
dir_result_file_name = fullfile(dir_result_root, ...
    ckpt, ...
    image, ...
    append("qrnn3d_", image, "_Case", num2str(idx_noise_condition), ".mat") ...
    );

load(dir_result_file_name, "HSI_restored", "HSI_noisy")
other_result.HSI_input = HSI_noisy;

[HSI_clean, hsi] = Load_HSI(image);
noise_seed = "default";
[HSI_noisy, deg] = Generate_obsv_for_denoising(HSI_clean, deg, noise_seed);

HSI_clean = single(HSI_clean);
HSI_noisy = single(HSI_noisy);


removed_noise.all = HSI_noisy - HSI_restored;


%% Plotting results
val_mpsnr  = calc_MPSNR(HSI_restored, HSI_clean);
val_mssim  = calc_MSSIM(HSI_restored, HSI_clean);
val_sam    = calc_SAM(HSI_restored, HSI_clean);

[vals_psnr_per_band, vals_ssim_per_band] ...
    = calc_PSNR_SSIM_per_band(HSI_restored, HSI_clean);


fprintf("\n~~~ SETTINGS ~~~\n");
fprintf("Method: %s\n", name_method);
fprintf("Image: %s Size: (%d, %d, %d)\n", image, hsi.n1, hsi.n2, hsi.n3);
fprintf("Gaussian sigma: %g\n", deg.gaussian_sigma);
fprintf("Sparse rate: %g\n", deg.sparse_rate);
fprintf("Stripe rate: %g\n", deg.stripe_rate);
fprintf("Stripe intensity: %g\n", deg.stripe_intensity);
fprintf("Deadline rate: %g\n", deg.deadline_rate)

fprintf("~~~ RESULTS ~~~\n");
fprintf("MPSNR: %#.4g\n", val_mpsnr);
fprintf("MSSIM: %#.4g\n", val_mssim);
fprintf("SAM  : %#.4g\n", val_sam);


%% Saving results
dir_save_method_folder = fullfile(...
    dir_save_folder, ...
    append("denoising_", image), ...
    append("g", num2str(deg.gaussian_sigma), "_ps", num2str(deg.sparse_rate), ...
        "_pt", num2str(deg.stripe_rate), "_pd", num2str(deg.deadline_rate)), ...
    name_method ...
);

dir_save_comp_method_folder = fullfile(...
    dir_save_comp_folder, ...
    append("denoising_", image), ...
    append("g", num2str(deg.gaussian_sigma), "_ps", num2str(deg.sparse_rate), ...
        "_pt", num2str(deg.stripe_rate), "_pd", num2str(deg.deadline_rate)), ...
    name_method ...
);

mkdir(dir_save_method_folder);
mkdir(dir_save_comp_method_folder);


save(fullfile(dir_save_method_folder, append(name_params_savetext, ".mat")), ...
    "HSI_clean", "HSI_noisy", "hsi", "deg", "image", ...
    "HSI_restored", "removed_noise", ...
    "val_mpsnr", "val_mssim", "val_sam", ...
    "vals_psnr_per_band", "vals_ssim_per_band", ...
    "params", "other_result", ...
    "-v7.3", "-nocompression" ...
);
save(fullfile(dir_save_comp_method_folder, append(name_params_savetext, ".mat")), ...
    "HSI_restored", "params", "val_mpsnr", "val_mssim", "val_sam", ...
    "-v7.3", "-nocompression" ...
);

close all


fprintf("Move (%d/%d)\n", idx_data, total_data);
idx_data = idx_data + 1;


end
end
end
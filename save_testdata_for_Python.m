clear
close all;

addpath(genpath('sub_functions'))

%% Selecting conditions
noise_conditions = { ...
    %g      ps     pt     t_int pd
    {0.1,   0,     0,     0,    0   }, ...
    {0.1,   0.05,  0,     0,    0   }, ...
    {0.1,   0,     0.05,  0.5,  0   }, ...
    {0.1,   0,     0,     0,    0.01}, ...
    {0.1,   0.05,  0.05,  0.5,  0.01}, ...
    {0.2,   0.05,  0.05,  0.5,  0.01}, ...
    };

% idc_noise_conditions = 1:size(noise_conditions, 2);
idc_noise_conditions = 6;


images = {...
    "JasperRidge", ...
    "PaviaU", ...
    };

idc_images = 1:numel(images);

dir_save_root = "H:\マイドライブ\MATLAB_Share\Dataset_GeoSSTV_OJSP";


idx_data = 1;
total_data = length(idc_noise_conditions) * length(idc_images);

for idx_noise_condition = idc_noise_conditions
    for idx_image = idc_images
        %% Generating observation
        deg.gaussian_sigma      = noise_conditions{idx_noise_condition}{1};
        deg.sparse_rate         = noise_conditions{idx_noise_condition}{2};
        deg.stripe_rate         = noise_conditions{idx_noise_condition}{3};
        deg.stripe_intensity    = noise_conditions{idx_noise_condition}{4};
        deg.deadline_rate       = noise_conditions{idx_noise_condition}{5};
        image = images{idx_image};

        [HSI_clean, hsi] = Load_HSI(image);
        noise_seed = "default";
        [HSI_noisy, deg] = Generate_obsv_for_denoising(HSI_clean, deg, noise_seed);


        gt = HSI_clean;
        input = HSI_noisy;

        % input(input>1) = 1;
        % input(input<0) = 0;


        dir_save_folder = fullfile(dir_save_root, ...
            image, ...
            append("Case", num2str(idx_noise_condition)) ...
            );
        mkdir(dir_save_folder)

        save(fullfile(dir_save_folder, "data.mat"), "input", "gt")

        fprintf("Save (%d/%d)\n", idx_data, total_data);
        idx_data = idx_data + 1;


    end
end
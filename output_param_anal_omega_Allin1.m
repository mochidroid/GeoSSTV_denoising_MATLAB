clear;
close all;

addpath('./func_metrics')
addpath(genpath('sub_functions'))


%% Switching operator
is_save = 1;


%% Setting conditions
noise_conditions = { ...
    %g      ps     pt     t_int pd
    {0.1,   0,     0,     0,    0   }, ...  
    {0.1,   0.05,  0,     0,    0   }, ... 
    {0.1,   0,     0.05,  0.5,  0   }, ... 
    {0.1,   0,     0,     0,    0.01}, ... 
    {0.1,   0.05,  0.05,  0.5,  0.01}, ... 
};

idx_noise_condition = 1;

images = {...
    "JasperRidge", ...
    "PaviaU", ...
};

images_legend = { ...
    'Jasper Ridge', ...
    'Pavia University', ...
};

idc_images = 1:numel(images);
% idc_images = 1;
num_images = numel(idc_images);


%% Selecting common parameters
rho = 0.98;
stopcri_idx = 5;

name_method = "GeoSSTV";
get_params_savetext = @(lambda) ...
    sprintf("o%g_r%.2f_stop1e-%d", lambda, rho, stopcri_idx);


lambda_list = [0:0.005:0.07];
num_params_lambda = numel(lambda_list);


vals_mpsnr_list = NaN([num_params_lambda, num_images]);
vals_mssim_list = NaN([num_params_lambda, num_images]);

load("dir_save_comp_folder.mat", "dir_save_comp_folder")


%% Loading results
for idx_image = idc_images
    %% Generating observation
    deg.gaussian_sigma      = noise_conditions{idx_noise_condition}{1};
    deg.sparse_rate         = noise_conditions{idx_noise_condition}{2};
    deg.stripe_rate         = noise_conditions{idx_noise_condition}{3};
    deg.stripe_intensity    = noise_conditions{idx_noise_condition}{4};
    deg.deadline_rate       = noise_conditions{idx_noise_condition}{5};
    image = images{idx_image};
    
    [HSI_clean, hsi] = Load_HSI(image);
    HSI_clean = single(HSI_clean);

    dir_result_folder = fullfile(...
        dir_save_comp_folder, ...
        append("denoising_", image), ...
        append("g", num2str(deg.gaussian_sigma), "_ps", num2str(deg.sparse_rate), ...
                "_pt", num2str(deg.stripe_rate), "_pd", num2str(deg.deadline_rate)), ...
        name_method ...
    );
    
    
    %% Loading MPSNRs and MSSIMs
    for idx_lambda = 1:num_params_lambda
        lambda = lambda_list(idx_lambda);
        params_savetext = get_params_savetext(lambda);

        
        load(fullfile(dir_result_folder, append(params_savetext, ".mat")), ...
            "val_mpsnr", "val_mssim");

        vals_mpsnr_list(idx_lambda, idx_image) = val_mpsnr;
        vals_mssim_list(idx_lambda, idx_image) = val_mssim;
    end
end

%% Plotting and saving result
% plot_style = {'LineWidth', 2.0};
plot_style = {'LineWidth', 5.0};
label_style = {'FontSize', 30};
% label_style = {'FontSize', 300};
% figure_style = {'Position', [100, 100, 800, 600]};
figure_style = {'Position', [100, 100, 600, 500]};
ax_fontsize = 17;


figure(figure_style{:});
hold on
for idx_image = idc_images
    image_legend = images_legend{idx_image};
    plot(lambda_list, vals_mpsnr_list(:, idx_image), ...
        'DisplayName', image_legend, plot_style{:});
end


xlabel('\omega', label_style{:});
ylabel('MPSNR', label_style{:});
xlim([0, 0.07]);
% ylim([33, 37]);
xticks(0:0.01:0.07);
ax = gca;
ax.FontSize = ax_fontsize;
grid on
legend('show')

hold off

if exist('is_save', 'var') && is_save == 1
    dir_output_folder = fullfile(...
        dir_save_comp_folder, ...
        "GeoSSTV_OJSP", ...
        "param_anal");
    mkdir(dir_output_folder)

    % saveas(gcf, append(dir_output_folder, 'omega_mpsnr.png'))
    exportgraphics(gcf, fullfile(dir_output_folder, 'omega_mpsnr.eps'), 'ContentType', 'image')
end


figure(figure_style{:});
hold on
for idx_image = idc_images
    image_legend = images_legend{idx_image};
    plot(lambda_list, vals_mssim_list(:, idx_image), ...
        'DisplayName', image_legend, plot_style{:});
end
xlabel('\omega', label_style{:});
ylabel('MSSIM', label_style{:});
xlim([0, 0.07]);
ylim([0.915, 0.95]) 
xticks(0:0.01:0.07);
ax = gca;
ax.FontSize = ax_fontsize;
grid on
legend('show')

if exist('is_save', 'var') && is_save == 1
    % saveas(gcf, append(dir_output_folder, 'omega_mssim.png'))
    exportgraphics(gcf, fullfile(dir_output_folder, 'omega_mssim.eps'), 'ContentType', 'image')
end

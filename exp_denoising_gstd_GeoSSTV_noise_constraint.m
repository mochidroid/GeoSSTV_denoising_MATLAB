clear
close all;

addpath(genpath("sub_functions"))
addpath("func_metrics")
addpath("methods")

fprintf("******* initium *******\n");

%% Selecting conditions
% noise_conditions = { ...
%     %g      ps     pt     t_int pd
%     {0.05,  0.05,  0.05,  0.5,  0.01}, ... % g0.05 ps0.05 pt0.05 pd0.001
%     {0.1,   0.05,  0.05,  0.5,  0.01}, ... % g0.1 ps0.05 pt0.05 pd0.001
%     {0.05,  0.05,  0.05,  0.5,  0.03}, ... % g0.1 ps0.05 pt0.05 pd0.001
%     {0.1,   0.05,  0.05,  0.5,  0.03}, ... % g0.1 ps0.05 pt0.05 pd0.001
%     {0.05,  0.05,  0,     0,    0   }, ...
%     {0.1,   0.05,  0,     0,    0   }, ...
%     {0.05,  0.05,  0.05,  0.5,  0   }, ...
%     {0.1,   0.05,  0.05,  0.5,  0   }, ...
%     {0.2,   0.05,  0,     0,    0   }, ...
%     {0.2,   0.05,  0.05,  0.5,  0   }, ...
% };

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
idc_noise_conditions = 5;


images = {...
    "JasperRidge", ...
    "PaviaU", ...
    };

% images = {...
%     "JasperRidge64", ...
%     "PaviaU64", ...
% };

% is_flipped = {...
%     0, ...
%     1, ...
%     0, ...
%     1, ...
% };

% idc_images = 1:numel(images);
idc_images = 1;


idx_exp = 0;
total_exp = length(idc_noise_conditions) * length(idc_images);


%% Setting common parameters
% rhos = {0.93, 0.95, 0.98};
rhos = {0.90};
% rhos = {0.98};

% epsilon_rho = 0.01;

% % epsilon = epsilon_rho * sqrt(hsi.N * (1 - deg.sparse_rate));
% epsilon = rhos * deg.gaussian_sigma * sqrt(hsi.N * (1 - deg.sparse_rate));
% alpha = rhos * (0.5 * hsi.N * deg.sparse_rate);
% beta = rhos * hsi.N * (1 - deg.sparse_rate) ...
%     * deg.stripe_rate * deg.stripe_intensity / 2;

stopcri_idx = 5;
stopcri = 10 ^ -stopcri_idx;

maxiter = 20000;
% maxiter = 5;

load("dir_save_folder.mat", "dir_save_folder");
load("dir_save_comp_folder.mat", "dir_save_comp_folder")


%% Setting each methods info
% GeoSSTV
% GeoSSTV.lambda = {0.001, 0.005, 0.01, 0.03};
GeoSSTV.lambda = {0.01, 0.02, 0.03, 0.05};

% Multipliers for Ablation Study
mults_epsilon = {0.5, 1, 2, 5};
mults_alpha = {0.5, 1, 2, 5};
mults_beta = {0.5, 1, 2, 5};

addpath("./methods/GeoSSTV")

clear methods_info; % Initialize methods_info completely

methods_info(1) = struct( ...
    "name", "GeoSSTV_g", ...
    "func", @(HSI_clean, HSI_noisy, params, deg) ...
    func_GeoSSTV_g_for_denoising_v2(HSI_clean, HSI_noisy, params), ...
    "param_names", {{"lambda", "maxiter", "stopcri", "rho_radius", "mult_epsilon"}}, ...
    "params", {{GeoSSTV.lambda, maxiter, stopcri, rhos, mults_epsilon}}, ...
    "get_params_savetext", @(params) ...
    sprintf("o%g_r%.2f_mE%g_stop1e-%d", params.lambda, params.rho_radius, params.mult_epsilon, stopcri_idx), ...
    "enable", true ...
    );

methods_info(end+1) = struct( ...
    "name", "GeoSSTV_gs", ...
    "func", @(HSI_clean, HSI_noisy, params, deg) ...
    func_GeoSSTV_gs_for_denoising_v2(HSI_clean, HSI_noisy, params), ...
    "param_names", {{"lambda", "maxiter", "stopcri", "rho_radius", "mult_epsilon", "mult_alpha"}}, ...
    "params", {{GeoSSTV.lambda, maxiter, stopcri, rhos, mults_epsilon, mults_alpha}}, ...
    "get_params_savetext", @(params) ...
    sprintf("o%g_r%.2f_mE%g_mA%g_stop1e-%d", params.lambda, params.rho_radius, params.mult_epsilon, params.mult_alpha, stopcri_idx), ...
    "enable", true ...
    );

methods_info(end+1) = struct( ...
    "name", "GeoSSTV_gt", ...
    "func", @(HSI_clean, HSI_noisy, params, deg) ...
    func_GeoSSTV_gt_for_denoising_v2(HSI_clean, HSI_noisy, params), ...
    "param_names", {{"lambda", "maxiter", "stopcri", "rho_radius", "mult_epsilon", "mult_beta"}}, ...
    "params", {{GeoSSTV.lambda, maxiter, stopcri, rhos, mults_epsilon, mults_beta}}, ...
    "get_params_savetext", @(params) ...
    sprintf("o%g_r%.2f_mE%g_mB%g_stop1e-%d", params.lambda, params.rho_radius, params.mult_epsilon, params.mult_beta, stopcri_idx), ...
    "enable", true ...
    );

methods_info = methods_info([methods_info.enable]); % removing false methods
num_methods = numel(methods_info);


%% Running Expt.
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

        HSI_clean = single(HSI_clean);
        HSI_noisy = single(HSI_noisy);


        % if is_flipped{idx_image}
        %     HSI_clean = flip(HSI_clean, 2);
        %     HSI_noisy = flip(HSI_noisy, 2);
        %     image = image + "f";
        % end


        idx_exp = idx_exp + 1;

        i_method = 0;


        %% Running methods
        for idx_method = 1:num_methods
            name_method = methods_info(idx_method).name;
            func_method = methods_info(idx_method).func;
            params_name = methods_info(idx_method).param_names;
            params_cell = methods_info(idx_method).params;


            i_method = i_method + 1;


            [params_comb, num_params_comb] = ParamsList2Comb(params_cell);

            for idx_params_comb = 1:num_params_comb

                params = struct();
                for idx_params = 1:numel(params_name)
                    % Assigning parameters to the structure
                    params.(params_name{idx_params}) = params_comb{idx_params_comb}{idx_params};
                end

                % If rho_radius exists, calculate epsilon, alpha, and beta
                if isfield(params, "rho_radius")
                    % Calculating the rate except sparse noise and deadline noise
                    r_sp = deg.sparse_rate;
                    r_dl = deg.deadline_rate;
                    dl_mean_width = 2;
                    cr_dl = 1 - exp(-dl_mean_width * r_dl); % cover_rate of deadline

                    mask_valid = (HSI_noisy > 0) & (HSI_noisy < 1);
                    mu = mean(HSI_noisy(mask_valid), 'all');

                    params.epsilon = params.rho_radius * deg.gaussian_sigma * sqrt(hsi.N * (1 - r_sp) * (1 - cr_dl));
                    params.alpha = params.rho_radius * hsi.N * (0.5*r_sp + mu*cr_dl);
                    params.beta    = params.rho_radius * hsi.N * (1 - r_sp) * (1 - cr_dl) ...
                        * (deg.stripe_rate * deg.stripe_intensity / 2);

                    % Applying ablation multipliers
                    if isfield(params, "mult_epsilon")
                        params.epsilon = params.epsilon * params.mult_epsilon;
                    end
                    if isfield(params, "mult_alpha")
                        params.alpha = params.alpha * params.mult_alpha;
                    end
                    if isfield(params, "mult_beta")
                        params.beta = params.beta * params.mult_beta;
                    end
                end

                name_params_savetext = methods_info(idx_method).get_params_savetext(params);


                fprintf("\n~~~ SETTINGS ~~~\n");
                fprintf("Method: %s\n", name_method);
                fprintf("Image: %s Size: (%d, %d, %d)\n", image, hsi.n1, hsi.n2, hsi.n3);
                fprintf("Gaussian sigma: %g\n", deg.gaussian_sigma);
                fprintf("Sparse rate: %g\n", deg.sparse_rate);
                fprintf("Stripe rate: %g\n", deg.stripe_rate);
                fprintf("Stripe intensity: %g\n", deg.stripe_intensity);
                fprintf("Deadline rate: %g\n", deg.deadline_rate);
                fprintf("Parameter settings: %s\n", name_params_savetext)
                fprintf("Methods: (%d/%d), Cases: (%d/%d), Params:(%d/%d)\n", ...
                    i_method, num_methods, idx_exp, total_exp, idx_params_comb, num_params_comb);

                [HSI_restored, removed_noise, other_result]...
                    = func_method(HSI_clean, HSI_noisy, params, deg);


                % Plotting results
                val_mpsnr  = calc_MPSNR(HSI_restored, HSI_clean);
                val_mssim  = calc_MSSIM(HSI_restored, HSI_clean);
                val_sam    = calc_SAM(HSI_restored, HSI_clean);

                fprintf("~~~ RESULTS ~~~\n");
                fprintf("MPSNR: %#.4g\n", val_mpsnr);
                fprintf("MSSIM: %#.4g\n", val_mssim);
                fprintf("SAM  : %#.4g\n", val_sam);

                [vals_psnr_per_band, vals_ssim_per_band] = calc_PSNR_SSIM_per_band(HSI_restored, HSI_clean);


                % Saving each result
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

            end

        end
    end
end

fprintf("******* finis *******\n");

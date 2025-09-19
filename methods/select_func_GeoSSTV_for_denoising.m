function [HSI_restored, removed_noise, other_result] ...
     = select_func_GeoSSTV_for_denoising(HSI_clean, HSI_noisy, params, deg)
addpath('./methods/GeoSSTV/')

% Selecting RISSTV function based on noise conditions
if deg.sparse_rate == 0 && deg.stripe_rate == 0 && deg.deadline_rate == 0
    % [HSI_restored, removed_noise, other_result] = ...
    %     func_GeoSSTV_g_for_denoising(HSI_clean, HSI_noisy, params);
    [HSI_restored, removed_noise, other_result] = ...
        func_GeoSSTV_g_for_denoising_v2(HSI_clean, HSI_noisy, params);

elseif deg.stripe_rate == 0
    [HSI_restored, removed_noise, other_result] = ...
        func_GeoSSTV_gs_for_denoising(HSI_clean, HSI_noisy, params);

elseif deg.sparse_rate == 0 && deg.deadline_rate == 0
    % [HSI_restored, removed_noise, other_result] = ...
    %     func_GeoSSTV_gt_for_denoising(HSI_clean, HSI_noisy, params);
    [HSI_restored, removed_noise, other_result] = ...
        func_GeoSSTV_gt_for_denoising_v2(HSI_clean, HSI_noisy, params);

else
    [HSI_restored, removed_noise, other_result] = ...
        func_GeoSSTV_gst_for_denoising(HSI_clean, HSI_noisy, params);
end
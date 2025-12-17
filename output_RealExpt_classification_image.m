clear;
close all;

addpath(genpath("sub_functions"))
addpath("func_metrics")


%% Switching operator
% is_show_cropped_image = 1;
% is_show_HSI = 1;  
is_output_image = 1;

%% Setting common parameters
images = {...
    "IndianPines", ...
    "Suwannee", ...
};

% Radius parameters; 1st: IndianPines, 2nd: Suwannee
epsilon = {30, 100};
alpha = {200, 800};
beta = {100, 5000};

stopcri_idx = 5;
stopcri = 10 ^ -stopcri_idx;

maxiter = 20000;


load("dir_save_comp_folder.mat", "dir_save_comp_folder");



%% Setting output parameters
idx_output = 1;

switch idx_output
    case 1
        idx_image = 1;

        % gain_restoration = 1;
        gain_restoration = 1.5; 
        
        gain_diff = 5;

        save_band = 32; % IndianPines
        % save_band = 40; % IndianPines
        % save_band = 88; % IndianPines

        crop_start_pos = [2, 30];
        crop_size = [20, 20];
        crop_expansion_rate = 2;
        crop_embed_tblr = 'br';
        
        arrow_head_pos = [20, 8];
        % arrow_head_pos = [38, 74];
        % arrow_head_pos = [53, 64];
        arrow_length = 15;
        arrow_handle_width = 3;
        arrow_head_width = 3;
        arrow_dir_tblr = "l";
        arrow_methods_idc = [];

    case 2
        idx_image = 2;

        gain_restoration = 1;
        % gain_restoration = 2;
        
        gain_diff = 8;

        save_band = 196; % Suwannee
        
        % crop_start_pos = [27, 80];
        crop_start_pos = [30, 80];
        crop_size = [20, 20];
        crop_expansion_rate = 2;
        crop_embed_tblr = 'bl';
        
        arrow_head_pos = [20, 8];
        % arrow_head_pos = [38, 74];
        % arrow_head_pos = [53, 64];
        arrow_length = 15;
        arrow_handle_width = 3;
        arrow_head_width = 3;
        arrow_dir_tblr = "l";
        arrow_methods_idc = [];

end

 
%% Setting each methods info
% SSTV
methods_info(1) = struct( ...
    "name", "SSTV", ...
    "output_name", "SSTV", ...
    "get_params_savetext", ...
        sprintf("e%g_a%g_b%g_stop1e-%d", ...
            epsilon{idx_image}, alpha{idx_image}, beta{idx_image}, ...
            stopcri_idx), ...
    "line_style", "--", ...
    "enable", true ...
);

% l0l1HTV
l0l1HTV.stepsize_reduction = 0.999; % 0.999
l0l1HTV.L10ball_th = 0.03; % 0.02, 0.03, 0.04

methods_info(end+1) = struct( ...
    "name", "l0l1HTV", ...
    "output_name", "$\llHTV$", ...
    "get_params_savetext", ...
        sprintf("e%g_a%g_b%g_sr%.5g_th%.2f_maxiter%d", ...
            epsilon{idx_image}, alpha{idx_image}, beta{idx_image}, ...
            l0l1HTV.stepsize_reduction, l0l1HTV.L10ball_th, maxiter), ...
    "line_style", "--", ...
    "enable", true ...
);

% HSSTV_L1
HSSTV.omega = 0.05;

methods_info(end+1) = struct( ...
    "name", "HSSTV_L1", ...
    "output_name", "HSSTV1", ...
    "get_params_savetext", ...
        sprintf("e%g_a%g_b%g_o%.2f_stop1e-%d", ...
            epsilon{idx_image}, alpha{idx_image}, beta{idx_image}, ...
            HSSTV.omega, stopcri_idx), ...
    "line_style", "--", ...
    "enable", true ...
);

% HSSTV_L12
methods_info(end+1) = struct( ...
    "name", "HSSTV_L12", ...
    "output_name", "HSSTV2", ...
    "get_params_savetext", ...
        sprintf("e%g_a%g_b%g_o%.2f_stop1e-%d", ...
            epsilon{idx_image}, alpha{idx_image}, beta{idx_image}, ...
            HSSTV.omega, stopcri_idx), ...
    "line_style", "--", ...
    "enable", true ...
);

% TPTV
TPTV.Rank = [7,7,5];
TPTV.initial_rank = 2;
TPTV.maxIter = 50; % 50, 100
TPTV.lambdas = 1e-3;% 5e-4, 1e-4, 1e-3, 1e-2, 1.5e-2

methods_info(end+1) = struct( ...
    "name", "TPTV", ...
    "output_name", "TPTV", ...
    "get_params_savetext", ...
        sprintf("maxiter%d_l%.4g", TPTV.maxIter, TPTV.lambdas), ...
    "line_style", "--", ...
    "enable", true ...
);

% QRNN3D
QRNN3D.ckpt = "complex"; % "complex", "paviaft"

methods_info(end+1) = struct( ...
    "name", "QRNN3D", ...
    "output_name", "QRNN3D", ...
    "get_params_savetext", ...
        sprintf("p_%s_0normalized", QRNN3D.ckpt), ...
    "line_style", "--", ...
    "enable", true ...
);

% GeoSSTV
GeoSSTV.lambda = 0.01; % 0.005, 0.01, 0.03, 0.05

% methods_info(end+1) = struct( ...
%     "name", "GeoSSTV", ...
%     "output_name", "GeoSSTV", ...
%     "get_params_savetext", ...
%         sprintf("e%g_a%g_b%g_o%.2f_stop1e-%d", ...
%             epsilon{idx_image}, alpha{idx_image}, beta{idx_image}, ...
%             GeoSSTV.lambda, stopcri_idx), ...
%     "line_style", "-", ...
%     "enable", true ...
% );

methods_info(end+1) = struct( ...
    "name", "GeoSSTV", ...
    "output_name", "GeoSSTV", ...
    "get_params_savetext", ...
        sprintf("e%g_a%g_b%g_o%.2f_stop1e-%d", ...
            20, alpha{idx_image}, beta{idx_image}, ...
            GeoSSTV.lambda, stopcri_idx), ...
    "line_style", "-", ...
    "enable", true ...
);

% methods_info(end+1) = struct( ...
%     "name", "GeoSSTV", ...
%     "output_name", "GeoSSTV_v2", ...
%     "get_params_savetext", ...
%         sprintf("e%g_a%g_b%g_o%.2f_stop1e-%d", ...
%             110, 850, beta{idx_image}, ...
%             GeoSSTV.lambda, stopcri_idx), ...
%     "line_style", "-", ...
%     "enable", true ...
% );

methods_info = methods_info([methods_info.enable]); % removing true methods
num_methods = numel(methods_info);


%% Loading clean and noisy HS images
image = images{idx_image};

[HSI_noisy, hsi] = Load_real_HSI(image);
HSI_noisy = single(HSI_noisy);


fprintf("Image: %s Size: (%d, %d, %d)\n", image, hsi.n1, hsi.n2, hsi.n3);


%% Cropping clean and noisy HS images
image_noisy = HSI_noisy(:,:,save_band)*gain_restoration;
image_noisy = Crop_Embed_image(image_noisy, ...
            crop_start_pos, crop_size, crop_expansion_rate, crop_embed_tblr);

if ~isempty(arrow_methods_idc)
    image_noisy = Embed_arrow(image_noisy, ...
                arrow_head_pos, arrow_length, ...
                arrow_handle_width, arrow_head_width, arrow_dir_tblr);
end


%% Loading and cropping result images
for idx_method = 1:num_methods 
    name_method = methods_info(idx_method).name;

    dir_result_folder = fullfile(...
        dir_save_comp_folder, ...
        append("denoising_", image), ...
        name_method ...
    );

    params_savetext = methods_info(idx_method).get_params_savetext;

    load(fullfile(dir_result_folder, append(params_savetext, ".mat")), ...
        "HSI_restored");

    methods_info(idx_method).HSI_restored = HSI_restored;


    % Cropping result images
    image_restored = HSI_restored(:,:,save_band)*gain_restoration;
    image_restored = Crop_Embed_image(image_restored, ...
                crop_start_pos, crop_size, crop_expansion_rate, crop_embed_tblr);

    
    if find(arrow_methods_idc == idx_method)
        image_restored = Embed_arrow(image_restored, ...
                    arrow_head_pos, arrow_length, ...
                    arrow_handle_width, arrow_head_width, arrow_dir_tblr);
    end

    methods_info(idx_method).image_restored         = image_restored;

end


%% Showing cropped images
if exist("is_show_cropped_image", "var") && is_show_cropped_image == 1
    cat_restored_image = image_noisy;

    for idx_method = 1:num_methods 
        cat_restored_image = cat(2, cat_restored_image, methods_info(idx_method).image_restored);
    end

    figure;
    imshow(cat_restored_image);
end


%% Showing HSI
if exist("is_show_HSI", "var") && is_show_HSI == 1
    cat_HSI = HSI_noisy;

    for idx_method = 1:num_methods
        cat_HSI = cat(2, cat_HSI, methods_info(idx_method).HSI_restored);
    end

    cat_diff_GeoSSTV = abs(repmat(methods_info(num_methods).HSI_restored, [1, num_methods+1, 1]) - cat_HSI) * gain_diff;

    implay(cat(1, cat_HSI, cat_diff_GeoSSTV));
end


%% Output restored image
if exist("is_output_image", "var") && is_output_image == 1
    dir_output_folder = fullfile( ...
        dir_save_comp_folder, ...
        "GeoSSTV_OJSP", ...
        sprintf("%s_b%d", image, save_band));
    mkdir(dir_output_folder)


    output_file_name = fullfile(...
            dir_output_folder, ...
            "image_noisy.png");

    imwrite(image_noisy, output_file_name, 'BitDepth', 8);
    

    for idx_method = 1:num_methods 
        name_method = methods_info(idx_method).name;
        image_restored = methods_info(idx_method).image_restored;

        output_file_name = fullfile(...
            dir_output_folder, ...
            sprintf("image_%s.png", name_method));

        imwrite(image_restored, output_file_name, 'BitDepth', 8);
    end

    fprintf("save dir: %s\n", dir_output_folder)
end



%% =========================================================================
%% Classification Experiment (Corrected for 16 Classes)
%% =========================================================================
fprintf('\n--- Starting Classification Experiment (16 Classes) ---\n');

% 1. Ground Truthの読み込みとクロップ
% -------------------------------------------------------------------------
% ユーザー指定のパス
gt_path = 'H:\マイドライブ\MATLAB_Share\HSIData\IndianPine\Indian_pines_gt.mat';

if isfile(gt_path)
    load(gt_path, 'indian_pines_gt');
    GT_original = indian_pines_gt;
else
    error('GT file not found. Please check the path.');
end

% クロップ処理 (画像データと合わせる)
GT_cropped = GT_original(1:end-25, 26:end);
[rows, cols] = size(GT_cropped);

% 2. データ準備
% -------------------------------------------------------------------------
train_ratio = 0.10;  % 10%を学習に使用
rng(42);             % シード固定

% 分類対象リストの作成
class_targets = struct('name', "Noisy", 'HSI', HSI_noisy);
for i = 1:numel(methods_info)
    class_targets(end+1).name = methods_info(i).name;
    class_targets(end).HSI = methods_info(i).HSI_restored;
end

metrics_table = table();
n_targets = numel(class_targets);

% カラーマップの設定（16クラスを見やすく表示するため）
% 背景(0)を黒、1-16をジェットカラーで表現
cmap_base = jet(16);
cmap_final = [0 0 0; cmap_base]; % 0番目を黒にする

figure('Name', 'Classification Maps 16 Classes', 'Position', [100, 100, 1200, 600]);

%% =========================================================================
%% 3. 分類ループ (Classification Map & Error Map & Output)
%% =========================================================================
% Figureの準備
fig_class = figure('Name', 'Classification Maps', 'Position', [50, 100, 1200, 500]);
fig_error = figure('Name', 'Error Maps (Red=Wrong)', 'Position', [50, 650, 1200, 500]);

% カラーマップ定義
cmap_class = [0 0 0; jet(16)];       % クラス用: 0番目(背景)=黒, 1-16=Jet
cmap_error = [0 0 0; 1 0 0];         % エラー用: 0(正解)=黒, 1(誤り)=赤

% 保存用フォルダの作成 (is_output_image = 1 の場合)
if exist("is_output_image", "var") && is_output_image == 1
    dir_out_class = fullfile(dir_save_comp_folder, "GeoSSTV_OJSP", sprintf("Classification_%s", image));
    if ~exist(dir_out_class, 'dir')
        mkdir(dir_out_class);
    end
    fprintf("Classification results will be saved to: %s\n", dir_out_class);
end

% Ground Truthのベクトル化
gt_vector = double(GT_cropped(:));
idx_labeled = find(gt_vector > 0);
y_labeled = gt_vector(idx_labeled);

% 学習・テスト分割
cv = cvpartition(y_labeled, 'HoldOut', 1 - train_ratio);
idx_train = cv.training;
idx_test = cv.test;

for i = 1:n_targets
    method_name = class_targets(i).name;
    fprintf('Classifying: %s ... ', method_name);
    
    % --- データ準備 ---
    img = class_targets(i).HSI;
    [nr, nc, nb] = size(img);
    X_vector = reshape(img, nr*nc, nb);
    X_labeled = X_vector(idx_labeled, :);
    
    % 正規化 (Z-score)
    X_labeled = zscore(X_labeled); 
    
    % 学習
    X_train = X_labeled(idx_train, :);
    y_train = y_labeled(idx_train);
    t = templateSVM('Standardize', true, 'KernelFunction', 'linear');
    Mdl = fitcecoc(X_train, y_train, 'Learners', t, 'Coding', 'onevsone');
    
    % テスト予測 & 精度算出
    X_test = X_labeled(idx_test, :);
    y_test = y_labeled(idx_test);
    pred_test = predict(Mdl, X_test);
    acc = calc_accuracies(y_test, pred_test);
    
    % テーブル保存
    new_row = table({method_name}, acc.OA, acc.AA, acc.Kappa, ...
        'VariableNames', {'Method', 'OA', 'AA', 'Kappa'});
    metrics_table = [metrics_table; new_row];
    
    % --- マップ作成 ---
    pred_labeled_all = predict(Mdl, X_labeled);
    
    % Classification Map 構築
    class_map = zeros(rows, cols);
    class_map(idx_labeled) = pred_labeled_all; 
    
    % Error Map 構築
    error_map = zeros(rows, cols);
    is_error = (GT_cropped > 0) & (GT_cropped ~= class_map);
    error_map(is_error) = 1; 
    
    % --- 画像の保存 (ここを追加) ---
    if exist("is_output_image", "var") && is_output_image == 1
        % 1. Classification Map 保存
        % uint8に変換することで、0がcmapの1行目(黒)、1が2行目...に対応します
        filename_c = fullfile(dir_out_class, sprintf("ClassMap_%s.png", method_name));
        imwrite(uint8(class_map), cmap_class, filename_c);
        
        % 2. Error Map 保存
        filename_e = fullfile(dir_out_class, sprintf("ErrorMap_%s.png", method_name));
        imwrite(uint8(error_map), cmap_error, filename_e);
    end
    
    % --- 画面表示 (Classification) ---
    figure(fig_class);
    subplot(2, ceil(n_targets/2), i);
    imagesc(class_map); 
    axis image off;
    title(sprintf('%s\nOA: %.2f%%', method_name, acc.OA));
    colormap(fig_class, cmap_class);
    clim([0 16]);
    
    % --- 画面表示 (Error) ---
    figure(fig_error);
    subplot(2, ceil(n_targets/2), i);
    imagesc(error_map);
    axis image off;
    error_rate = sum(error_map(:)) / length(idx_labeled) * 100;
    title(sprintf('%s\nError Rate: %.2f%%', method_name, error_rate));
    colormap(fig_error, cmap_error);
    clim([0 1]);
    
    fprintf('Done. OA: %.2f%%\n', acc.OA);
end

disp(metrics_table);

%% =========================================================================
%% 4. クラス別詳細精度の算出とテーブル出力
%% =========================================================================

% Indian Pinesのクラス名定義 (16クラス)
class_names = [
    "Alfalfa"; "Corn-notill"; "Corn-mintill"; "Corn"; ...
    "Grass-pasture"; "Grass-trees"; "Grass-pasture-mowed"; "Hay-windrowed"; ...
    "Oats"; "Soybean-notill"; "Soybean-mintill"; "Soybean-clean"; ...
    "Wheat"; "Woods"; "Bldgs-Grass-Trees-Drives"; "Stone-Steel-Towers"
];

% 結果格納用の行列初期化
num_classes = 16;
results_acc_per_class = zeros(num_classes, n_targets); % 各クラスの精度
results_summary = zeros(3, n_targets);                 % OA, AA, Kappa

% Train/Test数のカウント用文字列配列
train_test_str = strings(num_classes, 1);

% --- ループ処理 ---
for i = 1:n_targets
    fprintf('Evaluating: %s ... ', class_targets(i).name);
    
    % データ準備
    img = class_targets(i).HSI;
    [nr, nc, nb] = size(img);
    X_vector = reshape(img, nr*nc, nb);
    X_labeled = X_vector(idx_labeled, :);
    
    % 正規化 (Z-score)
    X_labeled = zscore(X_labeled); 
    
    % 学習・テストデータ分割
    X_train = X_labeled(idx_train, :);
    y_train = y_labeled(idx_train);
    X_test  = X_labeled(idx_test, :);
    y_test  = y_labeled(idx_test);
    
    % SVM学習 (線形カーネル)
    t = templateSVM('Standardize', true, 'KernelFunction', 'linear');
    Mdl = fitcecoc(X_train, y_train, 'Learners', t, 'Coding', 'onevsone');
    
    % テストデータ予測
    pred_test = predict(Mdl, X_test);
    
    % --- 詳細精度計算 ---
    % 混同行列 (必ず16x16になるようにOrderを指定)
    conf_mat = confusionmat(y_test, pred_test, 'Order', 1:num_classes);
    
    % クラスごとの精度 (対角成分 / 行の合計)
    class_acc = diag(conf_mat) ./ sum(conf_mat, 2);
    results_acc_per_class(:, i) = class_acc; % 0.0 ~ 1.0 の値
    
    % OA, AA, Kappa
    acc_summary = calc_accuracies_detailed(conf_mat);
    results_summary(1, i) = acc_summary.Kappa;
    results_summary(2, i) = acc_summary.OA; % 0.0 ~ 1.0
    results_summary(3, i) = acc_summary.AA; % 0.0 ~ 1.0
    
    % --- Train/Test数のカウント (初回ループ時のみ実行) ---
    if i == 1
        for c = 1:num_classes
            n_train = sum(y_train == c);
            n_test  = sum(y_test == c);
            train_test_str(c) = sprintf("%d/%d", n_train, n_test);
        end
    end
    
    fprintf('Done. (OA: %.4f)\n', acc_summary.OA);
end

%% =========================================================================
%% 5. テーブルの整形と表示
%% =========================================================================

% メソッド名を取得
method_names = [class_targets.name];

% テーブル作成
% 1. クラスごとの精度部分
T_class = array2table(results_acc_per_class, 'RowNames', class_names, 'VariableNames', method_names);

% 2. Train/Test列を追加 (先頭に)
T_class = addvars(T_class, train_test_str, 'Before', 1, 'NewVariableNames', 'Train_Test');

% 3. Summary部分 (Kappa, OA, AA)
summary_names = ["Kappa"; "OA"; "AA"];
% Train/Test列の場所には "-" を入れる
dummy_str = ["-"; "-"; "-"];
T_summary = array2table(results_summary, 'RowNames', summary_names, 'VariableNames', method_names);
T_summary = addvars(T_summary, dummy_str, 'Before', 1, 'NewVariableNames', 'Train_Test');

% 4. 結合して表示
T_final = [T_class; T_summary];

fprintf('\n==================================================================\n');
fprintf('CLASSIFICATION RESULTS ON INDIAN PINES DATA BY DIFFERENT METHODS\n');
fprintf('==================================================================\n');
disp(T_final);

% (オプション) CSVとして保存
% writetable(T_final, 'classification_results.csv', 'WriteRowNames', true);


%% 関数: 詳細精度計算 (比率 0-1 で出力)
function res = calc_accuracies_detailed(conf_mat)
    % OA
    oa = trace(conf_mat) / sum(conf_mat(:));
    
    % AA (各クラス精度の平均)
    class_acc = diag(conf_mat) ./ sum(conf_mat, 2);
    % NaN (テストデータがないクラス) を除外して平均
    valid_mask = ~isnan(class_acc);
    aa = mean(class_acc(valid_mask));
    
    % Kappa
    po = oa;
    pe = (sum(conf_mat, 1) * sum(conf_mat, 2)) / (sum(conf_mat(:))^2);
    kappa = (po - pe) / (1 - pe);
    
    res.OA = oa;
    res.AA = aa;
    res.Kappa = kappa;
end

%% 関数: 精度計算
function res = calc_accuracies(y_true, y_pred)
    conf_mat = confusionmat(y_true, y_pred);
    oa = trace(conf_mat) / sum(conf_mat(:)) * 100;
    
    % Class-wise accuracy (AA)
    class_counts = sum(conf_mat, 2);
    class_acc = diag(conf_mat) ./ class_counts;
    % サンプル数0のクラスがあれば除外
    valid_classes = class_counts > 0;
    aa = mean(class_acc(valid_classes)) * 100;
    
    % Kappa
    po = oa / 100;
    pe = (sum(conf_mat, 1) * sum(conf_mat, 2)) / (sum(conf_mat(:))^2);
    kappa = (po - pe) / (1 - pe);
    
    res.OA = oa;
    res.AA = aa;
    res.Kappa = kappa;
end
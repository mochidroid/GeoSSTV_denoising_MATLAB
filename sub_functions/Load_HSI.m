function  [HSI_clean] = Load_HSI(set_image)
rng('default')
%% Loading HSI
switch set_image
    case 'JasperRidge'
        load(fullfile("dataset", "jasperRidge2_R198.mat"));
        U_tmp = reshape(Y, [198, 100, 100]);
        HSI_clean = normalize01(permute(U_tmp, [2,3,1]));


    case 'PaviaUniversity'
        load(fullfile("dataset", "PaviaU.mat"))
        hsi.start_pos = [170, 200, 5];
        hsi_size = [140, 140, size(paviaU, 3)-hsi.start_pos(3)+1];
        hsi.end_pos = hsi.start_pos + hsi_size - 1;
        HSI_clean = normalize01(paviaU(hsi.start_pos(1):hsi.end_pos(1), hsi.start_pos(2):hsi.end_pos(2), hsi.start_pos(3):hsi.end_pos(3)));
end


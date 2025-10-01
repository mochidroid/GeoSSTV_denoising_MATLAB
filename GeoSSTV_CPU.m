%% Geometric Spatio-Spectral Total Variation for Hyperspectral Image Denoising and Destriping
%% =========================== First part notes===========================
% Author: Shingo Takemoto (takemoto.s.e908@m.isct.ac.jp)
% Last version: Oct. 1, 2025
% Article: S. Takemoto, S. Ono, 
%   ``Geometric Spatio-Spectral Total Variation for Hyperspectral Image Denoising and Destriping''
% -------------------------------------------------------------------------
%% =========================== Second part notes =========================== 
% INPUT:
%   HSI_noisy: noisy hyperspectral image of size n1*n2*n3 normalized to [0,1]
%   params: an option structure whose fields are as follows:           
%       alpha: radius of l_1 ball for sparse noise
%       beta: radius of l_1 ball for stripe noise
%       epsilon: radius of l_2 ball serving data-fidelity
%       omega: balancing parameter between 
%           first- and second-order differences
%       max_iter: maximum number of iterations
%       stop_cri: stopping criterion of P-PDS
%       dispiter: Period to display intermediate results
% OUTPUT:
%   restored_HSI: denoised hyperspectral image
%   removed_noise: removed noise
%   iteration: number of P-PDS iteration
%  ========================================================================

function [HSI_restored, removed_noise, iteration, converge_rate_U] ...
     = GeoSSTV_CPU(HSI_noisy, params)
fprintf('** Running GeoSSTV_CPU **\n');
HSI_noisy = single(HSI_noisy);
[n1, n2, n3] = size(HSI_noisy);

alpha       = single(params.alpha);
beta        = single(params.beta);
epsilon     = single(params.epsilon);
omage      = single(params.omega);
maxiter     = single(params.maxiter);
stopcri     = single(params.stopcri);

%% Setting params
dispiter    = unique([1:10, 1000:1000:maxiter]);
dispband    = round(n3/2);


%% Initializing primal and dual variables

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% primal variables
% U: clean HSI
% W1: interpolation (first-order difference)
% W2: interpolation (second-order difference)
% S: sparse noise(salt-and-pepper noise)
% T: stripe noise
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

U = zeros([n1, n2, n3], 'single');
W1 = zeros([n1, n2, n3, 2, 3], 'single');
W2 = zeros([n1, n2, n3, 2, 3], 'single');
S = zeros([n1, n2, n3], 'single');
T = zeros([n1, n2, n3], 'single');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% dual variables
% Y1: term of connection with Du and L*w
% Y2: term of connection with DDsu and L*w 
% Y2: term of l2ball
% Y3: term of stripe noise
%
% Y1 = D(U) - L*(W1)
% Y2 = D(Ds(U)) - L*(W2) 
% Y3 = U + S + T
% Y4 = Dv(T)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Y1 = zeros([n1, n2, n3, 2], 'single');
Y2 = zeros([n1, n2, n3, 2], 'single');
Y3 = zeros([n1, n2, n3], 'single');
Y4 = zeros([n1, n2, n3], 'single');


%% Setting operators
% Difference operators with Neumann boundary
D       = @(z) cat(4, z([2:end, end],:,:) - z, z(:,[2:end, end],:) - z);
Dt      = @(z) cat(1, -z(1, :, :, 1), -z(2:end-1, :, :, 1) + z(1:end-2, :, :, 1), z(end-1, :, :, 1)) ...
        + cat(2, -z(:, 1, :, 2), -z(:, 2:end-1, :, 2) + z(:, 1:end-2, :, 2), z(:, end-1, :, 2));

Dl      = @(z) z(:, :, [2:end, end]) - z;
Dlt     = @(z) cat(3, -z(:, :, 1), -z(:, :, 2:end-1) + z(:, :, 1:end-2), z(:, :, end-1));

Dv      = @(z) z([2:end, end], :, :) - z;
Dvt     = @(z) cat(1, -z(1, :, :), -z(2:(end-1), :, :) + z(1:(end-2), :, :), z(end-1, :, :));


%% Setting step size for P-PDS
gamma1_U    = single(1/(8 + 8*4 + 1)));
gamma1_W1   = single(1/4));
gamma1_W2   = single(1/4));
gamma1_S    = single(1));
gamma1_T    = single(1/(4 + 1)));
gamma2      = single(1/5));


%% main loop (P-PDS)
converge_rate_U = zeros([1, maxiter], 'single');
move_function_value = zeros([1, maxiter], 'single');
fprintf('~~~ P-PDS STARTS ~~~\n');

for i = 1:maxiter   
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Updating U
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    U_tmp   = U - gamma1_U*(Dt(Y1) + Dlt(Dt(Y2)) + Y3);
    U_next  = ProjBox(U_tmp, 0, 1);
    U_res   = 2*U_next - U;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Updating S
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    S_tmp   = S - gamma1_S*Y3;
    S_next  = ProjFastL1Ball(S_tmp, alpha);
    S_res = 2*S_next - S;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Updating T
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    T_tmp   = T - gamma1_T*(Y3 + Dvt(Y4));
    T_next  = ProjFastL1Ball(T_tmp, beta);
    T_res = 2*T_next - T;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Updating W1
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    W1_tmp   = W1 + gamma1_W1*(opL(Y1));
    W1_next  = Prox_l12norm_d4(W1_tmp, omage*gamma1_W1);
    W1_res   = 2*W1_next - W1;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Updating W2
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    W2_tmp   = W2 + gamma1_W2*(opL(Y2));
    W2_next  = Prox_l12norm_d4(W2_tmp, gamma1_W2);
    W2_res   = 2*W2_next - W2;

    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Updating Y1
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    Y1_next  = Y1 + gamma2*(D(U_res)-opLadj(W1_res));

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Updating Y2
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    Y2_next  = Y2 + gamma2*(D(Dl(U_res))-opLadj(W2_res));

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Updating Y3
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    Y3_tmp  = Y3 + gamma2*(U_res + S_res + T_res);    
    Y3_next = Y3_tmp - gamma2*ProjL2ball(Y3_tmp/gamma2, HSI_noisy, epsilon);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Updating Y4
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    Y4_next = Y4 + gamma2*Dv(T_res);

    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Calculating error
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    converge_rate_U(i) = norm(U_next(:) - U(:),2)/norm(U(:),2);
    move_function_value(i) = omage*sum(sqrt(sum(W1.^2,4)), 'all') + sum(sqrt(sum(W2.^2,4)), 'all');
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Updating all variables
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    U   = U_next;
    W1  = W1_next;
    W2  = W2_next;
    S   = S_next;
    T   = T_next;

    
    Y1  = Y1_next;
    Y2  = Y2_next;
    Y3  = Y3_next;
    Y4  = Y4_next;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Convergence checking
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if i>=2 && converge_rate_U(i) < stopcri % Checking convergence condition
        fprintf('Iter: %d, Error: %0.6f, Function value: %#.4g\n', ...
            i, converge_rate_U(i), move_function_value(i));
        break
    end
    if ismember(i, dispiter) % Displaying intermediate results
        fprintf('Iter: %d, Error: %0.6f, Function value: %#.4g\n', ...
            i, converge_rate_U(i), move_function_value(i));
        figure(1)
            subplot(2,2,1)
            imshow(HSI_noisy(:,:,dispband));
            title('Observed')
            
            subplot(2,2,2)
            imshow(gather(U(:,:,dispband)));
            title('Restored')
    
            subplot(2,2,3)
            semilogy(converge_rate_U(1:i));
            title('Converge rate U TV')
            
            subplot(2,2,4)
            semilogy(move_function_value(1:i));
            title('Function value')
            drawnow
    end
end

fprintf('~~~ P-PDS ENDS ~~~\n');

%% Organizing results for output
HSI_noisy                       = gather(HSI_noisy);
HSI_restored                    = gather(U);
iteration                       = gather(i);
removed_noise.sparse_noise      = gather(S);
removed_noise.stripe_noise      = gather(T);
removed_noise.gaussian_noise    = HSI_noisy - HSI_restored - ...
                                    removed_noise.sparse_noise - removed_noise.stripe_noise;
converge_rate_U                 = gather(converge_rate_U(1:iteration));

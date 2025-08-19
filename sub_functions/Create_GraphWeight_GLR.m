function [L_delta, info] = Create_GraphWeight_GLR(X, num_segments, k_lap, sigma_L)
% X: [N1 x N2 x N3]  HSI（ガイド用 or 観測）
% num_segments: セグメント数（ガイド画像のk-means数）
% k_lap: （任意）GLR用の相互kNN数（既定: min(10, K-1)）
% sigma_L:（任意）GLR用のガウシアン幅（既定: 距離分布の中央値から自動設定）
%
% 返り値:
%  L_delta: (K x K) のラプラシアン（K=N3-1; ノード=隣接バンド差分）
%  info: 構築に用いた中間情報（代表スペクトル等）

[n1, n2, n3] = size(X);
K = n3 - 1;
if nargin < 3 || isempty(k_lap), k_lap = min(10, K-1); end

% --- 1) ガイド画像をセグメント化 ---
guide_image = mean(X, 3);
labels = imsegkmeans(guide_image, num_segments);
segID = reshape(labels, n1, n2);

% --- 2) セグメント代表スペクトル r(s,:) を算出 ---
r = zeros(num_segments, n3);
counts = zeros(num_segments,1);
for s = 1:num_segments
    mask = (segID == s);
    counts(s) = nnz(mask);
    if counts(s) > 0
        tmp = X .* repmat(mask, [1,1,n3]);
        r(s,:) = sum(reshape(tmp, [], n3), 1) / counts(s);
    end
end

% --- 3) 隣接バンド差分（ノードの値） ---
Rdiff = r(:,2:end) - r(:,1:end-1);    % [num_segments x K]

% --- 4) セグメントサイズで重み付けした距離^2行列 D2（K x K）---
Wseg = diag(counts + eps);             % 画素数で重み付け
Gmat = sqrt(Wseg) * Rdiff;             % [num_segments x K]
GTG  = Gmat.' * Gmat;                  % [K x K]
nrm2 = diag(GTG);
D2   = (nrm2 + nrm2.') - 2*GTG;        % 距離^2
D2(1:K+1:end) = 0;

% --- 5) ガウシアンでノード間重み W_band ---
d2_vec = D2(triu(true(K),1));
d2_vec = d2_vec(d2_vec > 0);
if nargin < 4 || isempty(sigma_L)
    if isempty(d2_vec), sigma_L = 1;
    else, sigma_L = sqrt(median(d2_vec));
end
W_band = exp(-D2./(2*sigma_L^2));
W_band(1:K+1:end) = 0;

% 相互kNNで疎化（任意）
if k_lap > 0 && K > 1
    Wk = zeros(K);
    for i = 1:K
        [~, idx] = maxk(W_band(i,:), min(k_lap, K-1));
        Wk(i, idx) = W_band(i, idx);
    end
    W_band = max(Wk, Wk.');            % 対称化
end

% --- 6) ラプラシアン L = D - W ---
d = sum(W_band, 2);
L_delta = diag(d) - W_band;

% --- 7) 情報 ---
info.labels = segID;
info.representative_spectra = r;   % [num_segments x N3]
info.Rdiff = Rdiff;                % [num_segments x K]
info.segment_sizes = counts;
info.W_band = W_band;
info.sigma_L = sigma_L;
info.k_lap = k_lap;
end

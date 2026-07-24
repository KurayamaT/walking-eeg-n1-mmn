function stats_table3a_gmd_tanova(matpath)
% STATS_TABLE3A_GMD_TANOVA
%   Table 3A (spatial characterisation) — direct scalp-map comparison of the
%   MMN-window (130-200 ms) deviant-minus-standard maps, Walk-Water vs Walk-Clay:
%   global map dissimilarity (GMD) of the grand-average maps and a label-swap
%   TANOVA (topographic ANOVA). Each participant's 32-channel map is average-
%   referenced and GFP-normalised, then averaged; GMD is the RMS channel difference
%   of the two grand-average normalised maps (Murray et al., 2008).
%
%   The TANOVA p is the EXACT label-swap permutation: all 2^n within-participant
%   Water/Clay swaps are enumerated, so the value is deterministic and platform-
%   independent (the published .131 was a 10,000-sample Monte-Carlo estimate).
%
%   Usage: stats_table3a_gmd_tanova('path/to/erp_arrays_n22.mat')

if nargin < 1 || isempty(matpath)
    matpath = fullfile('..', 'data', 'erp_arrays_n22.mat');
end
L = load(matpath);
d = L.diff_array; t = L.times_ms(:); conds = cellstr(L.conditions);   % diff_array = deviant - standard
cw = find(strcmp(conds, 'walk_water'), 1); cc = find(strcmp(conds, 'walk_clay'), 1);
mmn = t >= 130 & t <= 200;

topo = @(c) squeeze(mean(d(mmn, :, c, :), 1)).';        % (subj, ch) window-mean map
A = gfp_normalize(topo(cw)); B = gfp_normalize(topo(cc));
v = ~(any(isnan(A), 2) | any(isnan(B), 2)); A = A(v, :); B = B(v, :);
n = size(A, 1);

obs = gmd(mean(A, 1), mean(B, 1));
p = tanova_exact(A, B);

fprintf('\n=== Table 3A (spatial): MMN-window direct scalp-map GMD / TANOVA, Water vs Clay ===\n');
fprintf('  n=%d  GMD=%.3f  TANOVA_p(exact)=%.3f   [target GMD=0.310  TANOVA=.134 (exact; .131 was Monte-Carlo)]\n', ...
        n, obs, p);
end


function M = gfp_normalize(M)     % average-reference then GFP-normalise each row (subject)
Mc = M - mean(M, 2);
M = Mc ./ sqrt(mean(Mc.^2, 2));
end


function g = gmd(a, b)            % RMS channel difference of two 1xch maps
g = sqrt(mean((a - b).^2));
end


function p = tanova_exact(A, B)
% Exact label-swap TANOVA. With per-subject-normalised maps fixed, a swap flips the
% sign of that subject's map difference D_i = A_i - B_i, so the permuted statistic is
%   sqrt(mean_ch( (1/n) sum_i s_i D_i )^2),  s in {-1,+1}^n.
% Enumerate all 2^n sign vectors (deterministic).
D = A - B; n = size(D, 1); N = 2^n;
obs = sqrt(mean(mean(D, 1).^2)); tol = 1e-12;
cnt = 0; batch = 2^16;
for st = 0:batch:N-1
    idx = (st:min(st + batch, N) - 1)';
    Sg = zeros(numel(idx), n);
    for b = 1:n; Sg(:, b) = 1 - 2 * bitget(idx, b); end
    stat = sqrt(mean(((Sg * D) / n).^2, 2));
    cnt = cnt + sum(stat >= obs - tol);
end
p = cnt / N;
end

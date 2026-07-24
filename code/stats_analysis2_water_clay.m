function stats_analysis2_water_clay(matpath)
% STATS_ANALYSIS2_WATER_CLAY
%   Analysis 2 (Walk-Water vs Walk-Clay). Reproduces the amplitude rows of Table 3
%   Panel A — the N1-window Stimulus x Condition interaction, the frontocentral-ROI
%   MMN amplitude (with its distribution-free corroboration and the equivalence/TOST),
%   from the deposited ERP arrays. Frontocentral ROI {Fz, FC1, FC2, Cz}. The
%   topographic-consistency rows of Panel A are reproduced by check_topo_consistency.m.
%
%   Permutation p is the EXACT sign-flip test: all 2^n sign reassignments are
%   enumerated, so the value is deterministic and platform-independent (no RNG seed).
%   Wilcoxon p uses the exact signed-rank distribution (MATLAB's default signrank
%   switches to a normal approximation for n > 15).
%
%   Usage: stats_analysis2_water_clay('path/to/erp_arrays_n22.mat')

if nargin < 1 || isempty(matpath)
    matpath = fullfile('..', 'data', 'erp_arrays_n22.mat');
end
S = load(matpath);
dev = S.dev_array; sd = S.std_array; t = S.times_ms(:);
chs = cellstr(S.ch_names); conds = cellstr(S.conditions);
roi = {'Fz', 'FC1', 'FC2', 'Cz'};
ridx = cellfun(@(c) find(strcmp(chs, c), 1), roi);
n1w = t >= 80 & t <= 130; mmnw = t >= 130 & t <= 200;
cw = find(strcmp(conds, 'walk_water'), 1);
cc = find(strcmp(conds, 'walk_clay'), 1);
amp = @(A, w) squeeze(mean(mean(A(w, ridx, :, :), 1), 2)).';   % (subj, cond)
devN1 = amp(dev, n1w); stdN1 = amp(sd, n1w);
devM  = amp(dev, mmnw); stdM  = amp(sd, mmnw);
n = size(devN1, 1);

fprintf('\n=== Table 3A (amplitude rows): Analysis 2 (Walk-Water vs Walk-Clay), n=%d ===\n', n);

% N1-window Stimulus x Condition interaction
interc = (devN1(:, cw) - devN1(:, cc)) - (stdN1(:, cw) - stdN1(:, cc));
FI = (mean(interc) / (std(interc) / sqrt(n)))^2; pI = 1 - fcdf(FI, 1, n - 1);
fprintf('N1-window Stimulus x Condition interaction: F(1,%d)=%.2f p=%.3f   [target 0.32, .577]\n', n - 1, FI, pI);

% MMN amplitude (deviant-evoked minus standard-evoked), Water vs Clay
mmnW = devM(:, cw) - stdM(:, cw); mmnC = devM(:, cc) - stdM(:, cc);
d = mmnW - mmnC;
[tA, pA, dzA] = ptest(d);
fprintf('MMN amplitude: Water %+.3f+/-%.3f  Clay %+.3f+/-%.3f  diff(W-C)=%+.3f\n', ...
        mean(mmnW), std(mmnW), mean(mmnC), std(mmnC), mean(d));
fprintf('  t(%d)=%+.2f p=%.3f d_z=%+.2f  perm=%.3f  Wilcoxon=%.3f   [target t=0.35 p=.727 dz=0.08 perm=.725 W=.824]\n', ...
        n - 1, tA, pA, dzA, exact_signflip(d), signrank(d, 0, 'method', 'exact'));

% TOST equivalence on the MMN-amplitude Water-Clay difference (SESOI = 0.5 d_z)
md = mean(d); sdd = std(d); se = sdd / sqrt(n); sesoi = 0.5 * sdd;
tl = (md + sesoi) / se; tu = (md - sesoi) / se;
pl = 1 - tcdf(tl, n - 1); pu = tcdf(tu, n - 1);
tc = tinv(0.95, n - 1); ci_lo = md - tc * se; ci_hi = md + tc * se;
fprintf('MMN-amplitude TOST (SESOI +/-%.3f uV): t_lower=%+.2f p=%.3f  t_upper=%+.2f p=%.3f  90%%CI[%+.3f,%+.3f]\n', ...
        sesoi, tl, pl, tu, pu, ci_lo, ci_hi);
fprintf('  [target t_lower=+2.70 p=.007  t_upper=-1.99 p=.030  90%%CI[-0.179,+0.271]  SESOI +/-0.307]\n');
end


function [t, p, dz] = ptest(d)
n = numel(d); t = mean(d) / (std(d) / sqrt(n)); p = 2 * (1 - tcdf(abs(t), n - 1)); dz = mean(d) / std(d);
end


function p = exact_signflip(d)
% Exact two-sided sign-flip permutation p: enumerate all 2^n sign reassignments of the
% difference scores and count |sum| >= |observed sum|. Deterministic, platform-
% independent (no RNG). For n = 22 this is 2^22 patterns (fast, exhaustive).
d = d(:); n = numel(d); N = 2^n; obs = abs(sum(d)); tol = 1e-9 * max(1, obs);
cnt = 0; batch = 2^16;
for st = 0:batch:N-1
    idx = (st:min(st + batch, N) - 1)';
    B = zeros(numel(idx), n);
    for b = 1:n; B(:, b) = bitget(idx, b); end
    cnt = cnt + sum(abs((1 - 2 * B) * d) >= obs - tol);
end
p = cnt / N;
end

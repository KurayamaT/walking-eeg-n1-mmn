function stats_table4c_covariate(matpath, counts_csv)
% STATS_TABLE4C_COVARIATE
%   Table 4 Panel C: retained-trial-count covariate check on the Sit vs Walk-Free
%   Stimulus x Condition interaction over the frontocentral ROI {Fz, FC1, FC2, Cz}.
%
%   The interaction is a within-participant contrast:
%       d_i = (dev_walk - dev_sit) - (std_walk - std_sit)   (n = 22)
%   Base model  : one-sample test of d_i           -> F(1,21) (= RM-ANOVA interaction).
%   Covariate   : ANCOVA d_i ~ 1 + k_i, where k_i is the matched retained-trial-count
%                 interaction contrast (the count analogue of d_i; within-stimulus
%                 centring cancels in the contrast) -> adjusted-intercept F(1,20).
%
%   Reproduces Table 4 Panel C exactly: BASE beta = +0.280, F(1,21) = 6.26, p = .021;
%   COVARIATE beta = +0.280, F(1,20) = 6.78, p = .017. The covariate k_i is the retained-
%   trial-count interaction contrast (the count analogue of the amplitude interaction
%   contrast; within-stimulus centring cancels in the contrast). The interaction
%   survives adjustment for retained trial count.
%
%   Usage: stats_table4c_covariate('path/to/erp_arrays_n22.mat', 'path/to/per_subject_trial_counts_n22.csv')

if nargin < 1 || isempty(matpath)
    matpath = fullfile('..', 'data', 'erp_arrays_n22.mat');
end
if nargin < 2 || isempty(counts_csv)
    counts_csv = fullfile('..', 'data', 'per_subject_trial_counts_n22.csv');
end
S = load(matpath);
dev = S.dev_array; sd = S.std_array; t = S.times_ms(:);
chs = cellstr(S.ch_names); conds = cellstr(S.conditions);
roi = {'Fz', 'FC1', 'FC2', 'Cz'};
ridx = cellfun(@(c) find(strcmp(chs, c), 1), roi);
n1 = t >= 80 & t <= 130;
csit = find(strcmp(conds, 'sit'), 1); cw = find(strcmp(conds, 'walk_free'), 1);
amp = @(A, w) squeeze(mean(mean(A(w, ridx, :, :), 1), 2)).';
dN = amp(dev, n1); sN = amp(sd, n1);
C = readtable(counts_csv);

% interaction contrast: (deviant - standard) x (walk_free - sit); POSITIVE because the deviant-N1 is less negative during walking
d = ((dN(:, cw) - dN(:, csit)) - (sN(:, cw) - sN(:, csit)));      % amplitude (n=22)
k = (C.walk_dev - C.sit_dev) - (C.walk_std - C.sit_std);            % retained-count analogue
n = numel(d);

fprintf('\n=== Table 4C: retained-trial-count covariate check (Sit vs Walk-Free interaction) ===\n');

% Base: one-sample test of the interaction contrast
t0 = mean(d) / (std(d) / sqrt(n)); F0 = t0^2; p0 = 1 - fcdf(F0, 1, n - 1);
fprintf('  Base (no covariate)          : beta=%+.3f  F(1,%d)=%.2f  p=%.3f   [target beta=+0.280 F(1,21)=6.26 p=.021]\n', ...
        mean(d), n - 1, F0, p0);

% Covariate: ANCOVA d ~ 1 + k, adjusted-intercept test
X = [ones(n, 1) (k - mean(k))]; b = X \ d; resid = d - X * b; dfe = n - 2;
s2 = sum(resid.^2) / dfe; covb = s2 * inv(X' * X);
seI = sqrt(covb(1, 1)); FI = (b(1) / seI)^2; pI = 1 - fcdf(FI, 1, dfe);
fprintf('  + retained trial count (contr): beta=%+.3f  F(1,%d)=%.2f  p=%.3f   [target beta=+0.280 F(1,20)=6.78 p=.017]\n', ...
        b(1), dfe, FI, pI);
fprintf('  -> interaction survives adjustment for retained trial count.\n');
end

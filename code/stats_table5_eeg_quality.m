function stats_table5_eeg_quality(matpath)
% STATS_TABLE5_EEG_QUALITY
%   Table 5: condition-wise EEG data-quality measures and the Sit vs Walk-Free
%   contrast, reproduced from the deposited canonical ERP arrays.
%
%   Four measures per participant and condition:
%     1. retained standard trials          (n_std_mat)
%     2. retained deviant trials           (n_dev_mat)
%     3. prestimulus averaged-waveform RMS (-50..0 ms)
%     4. late-window averaged-waveform RMS (350..450 ms)
%
%   The RMS measures are taken on the participant's ROI-averaged waveform
%   {Fz, FC1, FC2, Cz}, computed separately for the standard- and deviant-evoked
%   waveforms and then averaged over the two stimulus types. They are descriptive
%   proxies for data quality — not noise, signal-to-noise or single-trial estimates.
%   The late window lies outside every analysis window and is unaffected by
%   baseline correction.
%
%   The focal Sit vs Walk-Free contrast uses the manuscript's normality-first rule:
%   a Shapiro-Wilk gate at alpha = .05 on the within-participant difference selects
%   a paired t-test (d_z) or an exact Wilcoxon signed-rank test (matched-pairs
%   rank-biserial r), and Holm-Bonferroni is applied across the four measures.
%
%   Usage: stats_table5_eeg_quality('path/to/erp_arrays_n22.mat')

if nargin < 1 || isempty(matpath)
    matpath = fullfile('..', 'data', 'erp_arrays_n22.mat');
end
S = load(matpath);
t    = S.times_ms(:);
chs  = cellstr(S.ch_names);
conds = cellstr(S.conditions);
roi  = {'Fz', 'FC1', 'FC2', 'Cz'};
ridx = cellfun(@(c) find(strcmp(chs, c), 1), roi);

BASELINE = [-50 0];
LATE_WIN = [350 450];
bl   = t >= BASELINE(1) & t <= BASELINE(2);
late = t >= LATE_WIN(1) & t <= LATE_WIN(2);

% (time, channel, condition, subject) -> ROI-averaged waveform (time, cond, subj)
wstd = squeeze(mean(S.std_array(:, ridx, :, :), 2));
wdev = squeeze(mean(S.dev_array(:, ridx, :, :), 2));

nstd = S.n_std_mat; ndev = S.n_dev_mat;
if size(nstd, 1) ~= size(wstd, 3); nstd = nstd.'; ndev = ndev.'; end   % -> (subj, cond)

rmsw = @(W, m) squeeze(sqrt(mean(W(m, :, :).^2, 1))).';                % -> (subj, cond)
q_pre  = 0.5 * (rmsw(wstd, bl)   + rmsw(wdev, bl));
q_late = 0.5 * (rmsw(wstd, late) + rmsw(wdev, late));

names = {'Retained standard trials      ', ...
         'Retained deviant trials       ', ...
         'Prestimulus avg-waveform RMS  ', ...
         'Late-window avg-waveform RMS  '};
M = {nstd, ndev, q_pre, q_late};
targets = {'[target Sit 1028+/-18  WF 1022+/-22  t=-1.21 p=.240 pHolm=.719 dz=-0.26]', ...
           '[target Sit  206+/-18  WF  204+/-17  t=-0.48 p=.637 pHolm=1.000 dz=-0.10]', ...
           '[target Sit 0.091+/-0.053 WF 0.086+/-0.053 t=-0.53 p=.605 pHolm=1.000 dz=-0.11]', ...
           '[target Sit 0.359+/-0.204 WF 0.283+/-0.171 W=61 p=.033 pHolm=.132 r=-0.52]'};

ci_sit = find(strcmp(conds, 'sit'), 1);
ci_wf  = find(strcmp(conds, 'walk_free'), 1);
n = size(nstd, 1);

fprintf('\n=== Table 5: condition-wise EEG data quality (n=%d) ===\n', n);
fprintf('  ROI = %s; baseline %d..%d ms; late window %d..%d ms\n\n', ...
    strjoin(roi, ', '), BASELINE(1), BASELINE(2), LATE_WIN(1), LATE_WIN(2));

fprintf('  -- condition means (mean +/- SD) --\n');
fprintf('  %-30s %12s %12s %12s %12s\n', 'Measure', conds{:});
for k = 1:4
    A = M{k};
    fprintf('  %-30s', names{k});
    for c = 1:numel(conds)
        fprintf(' %6.3f+/-%-5.3f', mean(A(:, c)), std(A(:, c)));
    end
    fprintf('\n');
end

fprintf('\n  -- Sit vs Walk-Free (normality-gated, Holm across the four measures) --\n');
p_raw = nan(4, 1); rows = cell(4, 1);
for k = 1:4
    A = M{k};
    d = A(:, ci_wf) - A(:, ci_sit);
    [~, p_sw] = shapiro_wilk_norm(d);
    if p_sw >= 0.05
        dz = mean(d) / std(d);
        tval = mean(d) / (std(d) / sqrt(n));
        p_raw(k) = 2 * (1 - tcdf(abs(tval), n - 1));
        rows{k} = sprintf('t(%d) = %+.2f, dz = %+.2f', n - 1, tval, dz);
    else
        [Wstat, p_raw(k)] = wilcoxon_signed_rank_p(d);
        rk = tiedrank(abs(d));
        Wp = sum(rk(d > 0)); Wm = sum(rk(d < 0));
        r = -(Wm - Wp) / (Wp + Wm);
        rows{k} = sprintf('W = %g, r = %+.2f', Wstat, r);
    end
    rows{k} = sprintf('%s  [Shapiro-Wilk p = %.3f]', rows{k}, p_sw);
end
p_holm = holm_correct(p_raw);
for k = 1:4
    fprintf('  %-30s %-44s p = %.3f  pHolm = %.3f\n', names{k}, rows{k}, p_raw(k), p_holm(k));
    fprintf('  %-30s %s\n', '', targets{k});
end
fprintf(['\n  None of the four measures shows a statistically detectable difference between\n' ...
         '  Sit and Walk-Free after Holm correction. A non-significant difference is not\n' ...
         '  evidence of equivalence.\n\n']);
end

% ---------------------------------------------------------------- helpers

function padj = holm_correct(p)
p = p(:); m = numel(p);
[ps, idx] = sort(p);
adj = zeros(m, 1); run = 0;
for i = 1:m
    run = max(run, (m - i + 1) * ps(i));
    adj(i) = min(run, 1);
end
padj = zeros(m, 1); padj(idx) = adj;
end

function [W, p] = shapiro_wilk_norm(x)
% Shapiro-Wilk test, Royston (1992) AS R94; matches scipy.stats.shapiro to ~4 decimals.
x = x(~isnan(x)); x = sort(x(:)); n = length(x);
if n < 3 || n > 5000; W = NaN; p = NaN; return; end
i = (1:n)';
m = norminv((i - 0.375) / (n + 0.25));
ssm = m' * m; c = m / sqrt(ssm); u = 1 / sqrt(n);
an  = c(n)   + 0.221157*u - 0.147981*u^2 - 2.071190*u^3 + 4.434685*u^4 - 2.706056*u^5;
an1 = c(n-1) + 0.042981*u - 0.293762*u^2 - 1.752461*u^3 + 5.682633*u^4 - 3.582633*u^5;
if n > 5
    phi = (ssm - 2*m(n)^2 - 2*m(n-1)^2) / (1 - 2*an^2 - 2*an1^2);
    a = m / sqrt(phi); a(n) = an; a(n-1) = an1; a(1) = -an; a(2) = -an1;
else
    phi = (ssm - 2*m(n)^2) / (1 - 2*an^2);
    a = m / sqrt(phi); a(n) = an; a(1) = -an;
end
W = (a' * x)^2 / sum((x - mean(x)).^2);
ln = log(n);
if n <= 11
    gamma_v = -2.273 + 0.459 * n;
    w1 = -log(gamma_v - log(1 - W));
    mu_v = 0.5440 - 0.39978 * n + 0.025054 * n^2 - 0.0006714 * n^3;
    sigma_v = exp(1.3822 - 0.77857 * n + 0.062767 * n^2 - 0.0020322 * n^3);
else
    w1 = log(1 - W);
    mu_v = -1.5861 - 0.31082 * ln - 0.083751 * ln^2 + 0.0038915 * ln^3;
    sigma_v = exp(-0.4803 - 0.082676 * ln + 0.0030302 * ln^2);
end
p = 1 - normcdf((w1 - mu_v) / sigma_v);
end

function [W_stat, p] = wilcoxon_signed_rank_p(d)
% Exact Wilcoxon signed-rank p. MATLAB's signrank switches to the normal
% approximation above n = 15, which does not match scipy; force the exact method.
% Trial counts carry ties and zeros, for which the exact distribution is undefined;
% there the tie- and continuity-corrected normal approximation is used, as in scipy.
d = d(~isnan(d)); n0 = numel(d); d = d(d ~= 0); n = numel(d);
if n < 1; W_stat = NaN; p = NaN; return; end
rk = tiedrank(abs(d));
W_stat = min(sum(rk(d > 0)), sum(rk(d < 0)));
tied = numel(unique(abs(d))) < numel(d) || numel(d) < n0;
if tied
    mu = n * (n + 1) / 4;
    [~, ~, ic] = unique(abs(d));
    tie_sum = 0;
    for k = 1:max(ic); tk = sum(ic == k); tie_sum = tie_sum + tk^3 - tk; end
    sigma = sqrt(n * (n + 1) * (2 * n + 1) / 24 - tie_sum / 48);
    z = (abs(W_stat - mu) - 0.5) / sigma;
    p = 2 * (1 - normcdf(z));
else
    p = signrank(d, 0, 'method', 'exact');
end
end

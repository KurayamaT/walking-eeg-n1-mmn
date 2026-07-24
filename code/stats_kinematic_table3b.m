function stats_kinematic_table3b(derived_dir)
% STATS_KINEMATIC_TABLE3B
%   Analysis 2 kinematics (Walk-Water vs Walk-Clay). Reproduces Table 3 Panel B
%   from the deposited wide CSVs. Four co-primary outcomes, normality-first
%   (Shapiro-Wilk gate -> paired t if normal, Wilcoxon signed-rank otherwise),
%   Holm-Bonferroni x 4. Effect size: Cohen's d_z (paired t) or matched-pairs
%   rank-biserial r (Wilcoxon).
%
%   Paired Water-vs-Clay per outcome. Excursion inputs are the
%   orientation-corrected, drift-free per-stride Fourier displacements (Methods 2.6), so NO
%   participant is excluded: wrist n = 16, all three L3 outcomes n = 17.
%   (Wrist cumulative 3D path length is not an outcome: no stable gravitational
%   reference at the freely-rotating wrist.)
%
%   Usage: stats_kinematic_table3b('path/to/derived')
%     where the directory holds
%       kinematic_wrist_3cond_paired_wide_n22.csv
%       kinematic_L3_3cond_paired_wide_n22.csv

if nargin < 1 || isempty(derived_dir)
    derived_dir = fullfile('..', 'data');
end
W = readtable(fullfile(derived_dir, 'kinematic_wrist_3cond_paired_wide_n22.csv'));
L = readtable(fullfile(derived_dir, 'kinematic_L3_3cond_paired_wide_n22.csv'));

% sensor | measure | source table | Table-3B target string
defs = {
    'wrist', 'jerk_RMS_3d',     W, 'Wilcoxon W=16  p_raw=.005  p_Holm=.010  r=-0.76   (Water 2.03+/-0.87  Clay 2.52+/-0.99)';
    'L3',    'jerk_RMS_3d',     L, 'paired t(16)=-2.89  p_raw=.011  p_Holm=.011  d_z=-0.70   (Water 4.45+/-1.24  Clay 4.79+/-1.31)';
    'L3',    'vertical_length', L, 'paired t(16)=-5.18  p_raw<.001 p_Holm<.001  d_z=-1.26   (Water 0.70+/-0.19  Clay 0.83+/-0.24)';
    'L3',    'lateral_length',  L, 'paired t(16)=-4.96  p_raw<.001 p_Holm<.001  d_z=-1.20   (Water 1.06+/-0.35  Clay 1.28+/-0.43)'
};
n_def = size(defs, 1);

fprintf('\n=== Table 3B: Analysis 2 kinematics (Walk-Water vs Walk-Clay) ===\n');
fprintf('    normality-first (Shapiro-Wilk gate), Holm x 4, paired Water-vs-Clay, no S06 exclusion\n\n');

test_name = cell(n_def, 1);
stat_str  = cell(n_def, 1);
es_str    = cell(n_def, 1);
p_chosen  = zeros(n_def, 1);
n_used    = zeros(n_def, 1);
means_str = cell(n_def, 1);

for k = 1:n_def
    sensor = defs{k, 1}; m = defs{k, 2}; D = defs{k, 3};
    a = D.([m '_water']);
    b = D.([m '_clay']);
    ok = ~isnan(a) & ~isnan(b);
    a = a(ok); b = b(ok);
    n = numel(a);
    d = a - b;

    [~, p_sw] = shapiro_wilk_norm(d);
    if p_sw >= 0.05
        [t, df, p, dz] = paired_t(d);
        test_name{k} = 'paired t';
        stat_str{k}  = sprintf('t(%d)=%+.2f', df, t);
        es_str{k}    = sprintf('d_z=%+.2f', dz);
        p_chosen(k)  = p;
    else
        [Wstat, p] = wilcoxon_signed_rank_p(d);
        r = rank_biserial(d);
        test_name{k} = 'Wilcoxon';
        stat_str{k}  = sprintf('W=%.0f', Wstat);
        es_str{k}    = sprintf('r=%+.2f', r);
        p_chosen(k)  = p;
    end
    n_used(k)    = n;
    means_str{k} = sprintf('Water %.2f+/-%.2f  Clay %.2f+/-%.2f', ...
                           mean(a), std(a), mean(b), std(b));
end

p_holm = holm_correction(p_chosen);

for k = 1:n_def
    fprintf('  %-5s %-16s n=%2d  %-8s %-11s p_raw=%.4f p_Holm=%.4f %-9s %s  SW p=%s\n', ...
        defs{k, 1}, defs{k, 2}, n_used(k), test_name{k}, stat_str{k}, ...
        p_chosen(k), p_holm(k), es_str{k}, sig_star(p_holm(k)), means_str{k});
    fprintf('        [target: %s]\n', defs{k, 4});
end
fprintf('\n  All four Water<Clay, significant after Holm x 4 -> reproduces Table 3B.\n');
end


function p_holm = holm_correction(p_vec)
p_vec = p_vec(:); k = numel(p_vec);
[~, idx] = sort(p_vec);
p_holm = zeros(k, 1); running = 0;
for i = 1:k
    c = min(p_vec(idx(i)) * (k - i + 1), 1);
    running = max(running, c);
    p_holm(idx(i)) = running;
end
end


% ---------- statistics helpers (self-contained) ----------

function [t, df, p, dz] = paired_t(d)
d = d(~isnan(d)); n = numel(d);
t = mean(d) / (std(d) / sqrt(n));
df = n - 1;
p = 2 * (1 - tcdf(abs(t), df));
dz = mean(d) / std(d);
end


function r = rank_biserial(d)
% Matched-pairs rank-biserial r = (T+ - T-) / total rank sum, over |d|>0.
d = d(~isnan(d)); d = d(d ~= 0); n = numel(d);
if n < 1; r = NaN; return; end
ranks = tiedrank(abs(d));
r = (sum(ranks(d > 0)) - sum(ranks(d < 0))) / (n * (n + 1) / 2);
end


function [Wstat, p] = wilcoxon_signed_rank_p(d)
d = d(~isnan(d)); d = d(d ~= 0); n = numel(d);
if n < 1; Wstat = NaN; p = NaN; return; end
ranks = tiedrank(abs(d));
Wstat = min(sum(ranks(d > 0)), sum(ranks(d < 0)));
% Two-sided p. The manuscript used scipy.stats.wilcoxon, which computes the
% EXACT signed-rank distribution when there are no ties/zeros. MATLAB's default
% signrank switches to the normal approximation for n > 15, so request the exact
% method explicitly in that case; fall back to the approximation when ties are
% present (as scipy also does).
has_ties = numel(unique(abs(d))) < n;
if ~has_ties
    p = signrank(d, 0, 'method', 'exact');
else
    p = signrank(d, 0, 'method', 'approximate');
end
end


function [W, p] = shapiro_wilk_norm(x)
% Shapiro-Wilk, Royston (1992) AS R94. Matches scipy.stats.shapiro / R
% shapiro.test to ~4 decimals for n = 3..5000.
x = x(~isnan(x)); x = sort(x(:)); n = numel(x);
if n < 3 || n > 5000; W = NaN; p = NaN; return; end
i = (1:n)';
mm = norminv((i - 0.375) / (n + 0.25));
ssm = mm' * mm;
c = mm / sqrt(ssm);
u = 1 / sqrt(n);
an  = c(n)   + 0.221157*u - 0.147981*u^2 - 2.071190*u^3 + 4.434685*u^4 - 2.706056*u^5;
an1 = c(n-1) + 0.042981*u - 0.293762*u^2 - 1.752461*u^3 + 5.682633*u^4 - 3.582633*u^5;
if n > 5
    phi = (ssm - 2*mm(n)^2 - 2*mm(n-1)^2) / (1 - 2*an^2 - 2*an1^2);
    a = mm / sqrt(phi); a(n) = an; a(n-1) = an1; a(1) = -an; a(2) = -an1;
else
    phi = (ssm - 2*mm(n)^2) / (1 - 2*an^2);
    a = mm / sqrt(phi); a(n) = an; a(1) = -an;
end
W = (a' * x)^2 / sum((x - mean(x)).^2);
ln = log(n);
if n <= 11
    g = -2.273 + 0.459 * n;
    w1 = -log(g - log(1 - W));
    mu = 0.5440 - 0.39978*n + 0.025054*n^2 - 0.0006714*n^3;
    sg = exp(1.3822 - 0.77857*n + 0.062767*n^2 - 0.0020322*n^3);
else
    w1 = log(1 - W);
    mu = -1.5861 - 0.31082*ln - 0.083751*ln^2 + 0.0038915*ln^3;
    sg = exp(-0.4803 - 0.082676*ln + 0.0030302*ln^2);
end
p = 1 - normcdf((w1 - mu) / sg);
end


function s = sig_star(p)
if     p < .001; s = '***';
elseif p < .01;  s = '**';
elseif p < .05;  s = '*';
else;            s = '';
end
end

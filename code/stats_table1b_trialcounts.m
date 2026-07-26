function stats_table1b_trialcounts(derived_dir)
% STATS_TABLE1B_TRIALCOUNTS
%   Table 1 Panel B: retained trial counts per condition (after artifact rejection)
%   and the four-condition (Sit, Walk-Free, Walk-Water, Walk-Clay) omnibus test on
%   the counts, separately for Standard and Deviant stimuli. Reproduced from the
%   deposited per-participant trial-count table.
%
%   Test selection follows the manuscript's normality-first rule (Methods): the
%   Shapiro-Wilk test is applied to the additive-model residuals of the participant
%   x condition matrix, and the omnibus is a Friedman test when those residuals
%   depart from normality (p < .05) and a repeated-measures ANOVA otherwise. The
%   Standard counts fail the gate and are therefore tested with Friedman; the
%   Deviant counts pass and are tested with RM-ANOVA, Greenhouse-Geisser corrected.
%   Both tests are hugely non-significant, so the conclusion — retained trial counts
%   did not differ across conditions — does not depend on the choice.
%
%   Usage: stats_table1b_trialcounts('path/to/derived')

if nargin < 1 || isempty(derived_dir)
    derived_dir = fullfile('..', 'data');
end
T = readtable(fullfile(derived_dir, 'per_subject_trial_counts_n22.csv'));

fprintf('\n=== Table 1B: retained trial counts, 4-condition omnibus (n=%d) ===\n', height(T));
fprintf('    normality-first (Shapiro-Wilk gate on the additive-model residuals)\n\n');
labels = {'Standard', 'std', '[target means 1,014-1,028  SD 18-31  Friedman chi2(3)=4.34  p=.227  W=.066]';
          'Deviant',  'dev', '[target means 204-209      SD 15-18  RM-ANOVA F(3,63)=0.27  p_GG=.811]'};
for i = 1:size(labels, 1)
    s = labels{i, 2};
    Y = [T.(['sit_' s]), T.(['walk_' s]), T.(['water_' s]), T.(['clay_' s])];
    cm = mean(Y, 1); cs = std(Y, 0, 1);
    resid = Y - mean(Y, 1) - mean(Y, 2) + mean(Y(:));
    [~, sw] = shapiro_wilk_norm(resid(:));
    fprintf('  %-8s: means [%.0f, %.0f]  SD [%.0f, %.0f]\n', ...
        labels{i, 1}, min(cm), max(cm), min(cs), max(cs));
    if sw < 0.05
        [chi2, W, p] = friedman_1way(Y);
        fprintf('            Friedman chi2(3)=%.2f  p=%.3f  W=%.3f          [SW p=%.4f -> non-normal]\n', ...
            chi2, p, W, sw);
    else
        [F, df1, df2, p_unc, p_GG, epsGG] = rmanova_gg(Y);
        fprintf('            RM-ANOVA F(%d,%d)=%.2f  p_unc=%.3f  eps_GG=%.3f  p_GG=%.3f   [SW p=%.4f -> normal]\n', ...
            df1, df2, F, p_unc, epsGG, p_GG, sw);
    end
    fprintf('            %s\n', labels{i, 3});
end
end


function [chi2, W, p] = friedman_1way(Y)
% Friedman test with the tie correction used by scipy.stats.friedmanchisquare;
% W is Kendall's coefficient of concordance.
[n, k] = size(Y);
R = zeros(n, k);
for i = 1:n; R(i, :) = tiedrank(Y(i, :)); end
Rj = sum(R, 1);
chi2 = 12 / (n * k * (k + 1)) * sum(Rj.^2) - 3 * n * (k + 1);
tie_sum = 0;
for i = 1:n
    [~, ~, ic] = unique(Y(i, :));
    for c = 1:max(ic); tc = sum(ic == c); tie_sum = tie_sum + tc^3 - tc; end
end
chi2 = chi2 / (1 - tie_sum / (n * (k^3 - k)));
W = chi2 / (n * (k - 1));
p = 1 - chi2cdf(chi2, k - 1);
end


function [F, df1, df2, p_unc, p_GG, eps] = rmanova_gg(Y)
[n, k] = size(Y); grand = mean(Y(:)); cm = mean(Y, 1); sm = mean(Y, 2);
SSc = n * sum((cm - grand).^2); SSs = k * sum((sm - grand).^2);
SSt = sum((Y(:) - grand).^2); SSe = SSt - SSc - SSs;
df1 = k - 1; df2 = (n - 1) * (k - 1);
F = (SSc / df1) / (SSe / df2);
p_unc = 1 - fcdf(F, df1, df2);
S = cov(Y); M = zeros(k, k - 1);
for j = 1:k - 1; M(1:j, j) = 1 / j; M(j + 1, j) = -1; end
M = M ./ sqrt(sum(M.^2, 1));
Sp = M' * S * M; lam = real(eig(Sp)); lam(lam < 1e-12) = 1e-12;
eps = (sum(lam))^2 / ((k - 1) * sum(lam.^2));
eps = max(min(eps, 1), 1 / (k - 1));
p_GG = 1 - fcdf(F, df1 * eps, df2 * eps);
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

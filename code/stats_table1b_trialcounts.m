function stats_table1b_trialcounts(derived_dir)
% STATS_TABLE1B_TRIALCOUNTS
%   Table 1 Panel B: retained trial counts per condition (after artifact rejection)
%   and the four-condition (Sit, Walk-Free, Walk-Water, Walk-Clay) repeated-measures
%   ANOVA on the counts, separately for Standard and Deviant stimuli. Reproduced from
%   the deposited per-participant trial-count table.
%
%   NOTE on the ANOVA p-values: the Standard row's published p (.16) is the
%   Greenhouse-Geisser-corrected p (uncorrected would be .13); the Deviant row's
%   published .86 equals the UNCORRECTED p (GG-corrected = .83, HF = .85). The two
%   rows are hugely non-significant either way; this script prints both so the value
%   is unambiguous.
%
%   Usage: stats_table1b_trialcounts('path/to/derived')

if nargin < 1 || isempty(derived_dir)
    derived_dir = fullfile('..', 'data');
end
T = readtable(fullfile(derived_dir, 'per_subject_trial_counts_n22.csv'));

fprintf('\n=== Table 1B: retained trial counts, 4-condition RM-ANOVA (n=%d) ===\n', height(T));
labels = {'Standard', 'std', '[target means 1014-1028  SD 18-31  F(3,63)=1.95  p_GG=.16]';
          'Deviant',  'dev', '[target means 204-209    SD 15-18  F(3,63)=0.25  p=.86 (uncorrected)]'};
for i = 1:size(labels, 1)
    s = labels{i, 2};
    Y = [T.(['sit_' s]), T.(['walk_' s]), T.(['water_' s]), T.(['clay_' s])];
    cm = mean(Y, 1); cs = std(Y, 0, 1);
    [F, df1, df2, p_unc, p_GG, epsGG] = rmanova_gg(Y);
    fprintf('  %-8s: means [%.0f, %.0f]  SD [%.0f, %.0f]  F(%d,%d)=%.2f  p_unc=%.3f  eps_GG=%.3f  p_GG=%.3f\n', ...
        labels{i, 1}, min(cm), max(cm), min(cs), max(cs), df1, df2, F, p_unc, epsGG, p_GG);
    fprintf('           %s\n', labels{i, 3});
end
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

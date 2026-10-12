function r = rv_one_sample(x)
% RV_ONE_SAMPLE  One-sample summary used throughout the revision analyses.
%   Mean, Student-t 95% CI, d_z, t, F = t^2, parametric p, Shapiro-Wilk p,
%   exact Wilcoxon signed-rank p, and the test selected by the normality rule
%   (t test when Shapiro-Wilk p >= .05, otherwise exact Wilcoxon).
x = x(:); n = numel(x);
[~, p, ci, st] = ttest(x, 0);
r.n = n; r.mean = mean(x); r.sd = std(x);
r.ci_low = ci(1); r.ci_high = ci(2);
r.dz = mean(x) / std(x); r.t = st.tstat; r.F = st.tstat^2; r.p_t = p;
if n >= 12, [~, r.shapiro_p] = rv_shapiro_wilk(x); else, r.shapiro_p = NaN; end
r.wilcoxon_p = signrank(x, 0, 'method', 'exact');
if r.shapiro_p < .05
    r.selected_test = "exact Wilcoxon"; r.selected_p = r.wilcoxon_p;
else
    r.selected_test = "t"; r.selected_p = r.p_t;
end
end

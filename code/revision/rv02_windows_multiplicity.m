function R = rv02_windows_multiplicity(erp)
% RV02_WINDOWS_MULTIPLICITY  N1- and MMN-window tests for both contrasts with
%   Holm-Bonferroni sensitivity (Supplementary Table S4), the descriptive simple
%   effects (Results 3.1), and the mean N1 interaction contrast (Fig. 2e).
E = rv_load_erp(erp);
win = [80 130; 130 200]; wname = ["N1", "MMN"];
pairs = [1 2; 4 3]; cname = ["Sit_WalkFree", "Clay_Water"];
T = table();
for c = 1:2
    for k = 1:2
        d = rv_roi_mean(E.dev, E, win(k, 1), win(k, 2)); s = rv_roi_mean(E.std, E, win(k, 1), win(k, 2));
        x = (d(:, pairs(c, 2)) - d(:, pairs(c, 1))) - (s(:, pairs(c, 2)) - s(:, pairs(c, 1)));
        r = rv_one_sample(x);
        T = [T; table(cname(c), wname(k), r.n, r.mean, r.ci_low, r.ci_high, r.dz, r.F, r.p_t, ...
            r.shapiro_p, r.wilcoxon_p, r.selected_test, r.selected_p, 'VariableNames', ...
            {'contrast', 'window', 'n', 'mean_uv', 'ci_low_uv', 'ci_high_uv', 'dz', 'F', 'p_parametric', ...
             'shapiro_p', 'p_wilcoxon', 'selected_test', 'p_selected'})]; %#ok<AGROW>
        if c == 1 && k == 1, contrast = x; dN = d; sN = s; end
    end
end
for f = ["p_parametric", "p_selected"]
    p = T.(f);
    T.(f + "_holm2") = [rv_holm(p(1:2)); rv_holm(p(3:4))];
    T.(f + "_holm4") = rv_holm(p);
end
assert(abs(T.F(1) - 6.26) < .005);
S = table(); v = {dN, sN}; nm = ["deviant", "standard"];
for j = 1:2
    x = v{j}(:, 2) - v{j}(:, 1); r = rv_one_sample(x);
    S = [S; table(nm(j), mean(v{j}(:, 1)), mean(v{j}(:, 2)), r.mean, r.ci_low, r.ci_high, r.dz, r.t, r.p_t, ...
        r.wilcoxon_p, 'VariableNames', {'stimulus', 'sit_mean_uv', 'walk_free_mean_uv', 'change_uv', ...
        'ci_low_uv', 'ci_high_uv', 'dz', 't', 'p_t', 'p_wilcoxon'})]; %#ok<AGROW>
end
R.tests = T; R.simple = S;
R.contrasts = table(string(E.subjects(:)), contrast, 'VariableNames', {'subject', 'interaction_uv'});
end

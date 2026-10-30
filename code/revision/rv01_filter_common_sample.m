function R = rv01_filter_common_sample(paths)
% RV01_FILTER_COMMON_SAMPLE  High-pass sensitivity of the N1 interaction
%   (main-text Table 3 Panel B; Supplementary Table S2).
%   paths: ERP arrays reprocessed at 1.0, 0.5 and 0.1 Hz (cell array of three .mat files).
%   Available = participants with Sit and Walk-Free ERPs at that setting.
%   Common 19 = participants with both ERPs at all three settings.
%   Common 18 = the common sample without S11 (historical output, not reported in the current manuscript).
cut = [1 .5 .1];
C = nan(22, 3); ok = false(22, 3);
for f = 1:3
    E = rv_load_erp(paths{f});
    d = rv_roi_mean(E.dev, E, 80, 130); s = rv_roi_mean(E.std, E, 80, 130);
    C(:, f) = (d(:, 2) - d(:, 1)) - (s(:, 2) - s(:, 1));
    fin = squeeze(all(all(isfinite(E.dev(:, :, 1:2, :)) & isfinite(E.std(:, :, 1:2, :)), 1), 2));
    ok(:, f) = all(fin, 1)' & isfinite(C(:, f));
end
subj = E.subjects(:);
common19 = all(ok, 2); common18 = common19 & ~strcmp(subj, 'S11');
assert(sum(common19) == 19 && sum(common18) == 18);
S = table();
for f = 1:3
    sets = {ok(:, f), common19, common18}; names = ["available", "common19", "common18"];
    for g = 1:3
        r = rv_one_sample(C(sets{g}, f));
        S = [S; table(cut(f), names(g), r.n, r.mean, r.ci_low, r.ci_high, r.dz, r.F, r.p_t, ...
            r.shapiro_p, r.wilcoxon_p, r.selected_test, r.selected_p, 'VariableNames', ...
            {'cutoff_hz', 'sample', 'n', 'mean_uv', 'ci_low_uv', 'ci_high_uv', 'dz', 'F', 'p_F', ...
             'shapiro_p', 'wilcoxon_p', 'selected_test', 'selected_p'})]; %#ok<AGROW>
    end
end
P = table(); pr = [3 1; 3 2; 2 1]; sets = {common19, common18};
for g = 1:2
    for k = 1:3
        x = C(sets{g}, pr(k, 1)) - C(sets{g}, pr(k, 2)); m = rv_mean_ci(x);
        P = [P; table(sum(sets{g}), cut(pr(k, 1)), cut(pr(k, 2)), m.mean, m.ci_low, m.ci_high, ...
            'VariableNames', {'n', 'cutoff_a_hz', 'cutoff_b_hz', 'mean_a_minus_b_uv', 'ci_low_uv', 'ci_high_uv'})]; %#ok<AGROW>
    end
end
R.summary = S; R.paired = P;
R.subjects = table(string(subj), C(:, 1), C(:, 2), C(:, 3), common19, common18, ...
    'VariableNames', {'subject', 'C_1p0_uv', 'C_0p5_uv', 'C_0p1_uv', 'common19', 'common18'});
end

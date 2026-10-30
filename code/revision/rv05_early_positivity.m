function R = rv05_early_positivity(erp, hpf_paths)
% RV05_EARLY_POSITIVITY  Descriptive checks of the positivity preceding N1
%   (Supplementary Table S7; Supplementary Results S-R1; values underlying Fig. S6).
%   erp: canonical arrays; hpf_paths: arrays reprocessed at 1.0, 0.5 and 0.1 Hz.
E = rv_load_erp(erp);
early = E.t >= 40 & E.t <= 70;
stim = {"standard", E.std; "pooled_deviant", E.dev}; cond = ["Sit", "WalkFree"];
S = table(); P = table();
for a = 1:2
    x = squeeze(mean(stim{a, 2}(:, E.roi, :, :), 2));
    for c = 1:2
        w = squeeze(x(:, c, :)); z = mean(w(early, :), 1)'; gm = mean(w, 2);
        loo = (sum(z) - z) / 21; h = tinv(.975, 21) * std(z) / sqrt(22);
        ix = find(early); [pk, j] = max(gm(early));
        S = [S; table(stim{a, 1}, cond(c), 22, mean(z), median(z), std(z), mean(z) - h, mean(z) + h, sum(z > 0), ...
            min(loo), max(loo), E.t(ix(j)), pk, 'VariableNames', {'stimulus', 'condition', 'n', 'mean_uv', ...
            'median_uv', 'sd_uv', 'ci_low_uv', 'ci_high_uv', 'positive_subjects', 'loo_min_uv', 'loo_max_uv', ...
            'grand_peak_within_40_70_ms', 'grand_peak_uv'})]; %#ok<AGROW>
        P = [P; table(string(E.subjects(:)), repmat(stim{a, 1}, 22, 1), repmat(cond(c), 22, 1), z, ...
            'VariableNames', {'subject', 'stimulus', 'condition', 'early_40_70_mean_uv'})]; %#ok<AGROW>
    end
end
R.summary = S; R.participants = P;
names = ["canonical_1Hz", "rerun_1p0", "rerun_0p5", "rerun_0p1"]; files = [{erp}, hpf_paths(:)'];
X = cell(4, 2); common = true(22, 1);
for f = 1:4
    F = rv_load_erp(files{f});
    X{f, 1} = squeeze(mean(F.std(:, F.roi, :, :), 2)); X{f, 2} = squeeze(mean(F.dev(:, F.roi, :, :), 2));
    for a = 1:2, common = common & squeeze(all(all(isfinite(X{f, a}(:, 1:2, :)), 1), 2)); end
end
assert(sum(common) == 19);
win = {"early_descriptive", E.t >= 0 & E.t < 80; "N1", E.t >= 80 & E.t <= 130; "prestimulus", E.t >= -50 & E.t < 0};
D = table(); lab = ["standard", "pooled_deviant"];
for f = 1:4
    for a = 1:2
        for c = 1:2
            w = squeeze(X{f, a}(:, c, common)); y = mean(w, 2);
            for k = 1:3
                m = win{k, 2}; z = mean(w(m, :), 1)'; h = tinv(.975, 18) * std(z) / sqrt(19);
                ix = find(m); [~, j] = max(abs(y(m)));
                D = [D; table(names(f), cond(c), lab(a), win{k, 1}, 19, mean(z), mean(z) - h, mean(z) + h, ...
                    y(ix(j)), E.t(ix(j)), 'VariableNames', {'array', 'condition', 'stimulus', 'window', 'n', ...
                    'mean_uv', 'ci_low_uv', 'ci_high_uv', 'grand_average_abs_peak_uv', 'grand_average_abs_peak_ms'})]; %#ok<AGROW>
            end
        end
    end
end
R.common19 = D; R.common_subjects = string(E.subjects(common));
end

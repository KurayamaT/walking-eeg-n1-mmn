function [T, I] = rv04_latency_windows(erp)
% RV04_LATENCY_WINDOWS  Additional latency-window sensitivity for both contrasts
%   (Supplementary Table S6; Results 3.2). Window boundaries are inclusive.
E = rv_load_erp(erp);
W = {'N1', 80, 130; 'N1', 80, 120; 'N1', 80, 150; 'MMN', 130, 200; 'MMN', 100, 200; 'MMN', 100, 250; 'MMN', 150, 250};
con = {"Walk-Free minus Sit", 2, 1; "Water minus Clay", 3, 4};
T = table(); I = table();
for k = 1:size(W, 1)
    d = rv_roi_mean(E.dev, E, W{k, 2}, W{k, 3}); s = rv_roi_mean(E.std, E, W{k, 2}, W{k, 3});
    for c = 1:2
        a = con{c, 2}; b = con{c, 3};
        x = (d(:, a) - s(:, a)) - (d(:, b) - s(:, b)); r = rv_one_sample(x);
        T = [T; table(string(W{k, 1}), W{k, 2}, W{k, 3}, con{c, 1}, r.n, r.mean, r.ci_low, r.ci_high, r.dz, ...
            r.t, r.F, r.p_t, r.shapiro_p, r.wilcoxon_p, r.selected_test, r.selected_p, 'VariableNames', ...
            {'component', 'lo_ms', 'hi_ms', 'contrast', 'n', 'estimate_uv', 'ci95_low_uv', 'ci95_high_uv', 'dz', ...
             't', 'F', 'parametric_p', 'shapiro_p', 'wilcoxon_p', 'selected_test', 'selected_p'})]; %#ok<AGROW>
        I = [I; table(repmat(string(W{k, 1}), 22, 1), repmat(W{k, 2}, 22, 1), repmat(W{k, 3}, 22, 1), ...
            repmat(con{c, 1}, 22, 1), string(E.subjects(:)), x, 'VariableNames', ...
            {'component', 'lo_ms', 'hi_ms', 'contrast', 'subject', 'contrast_uv'})]; %#ok<AGROW>
    end
end
assert(abs(T.F(1) - 6.2574797646) < 1e-7);
end

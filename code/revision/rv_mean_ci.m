function r = rv_mean_ci(v)
% RV_MEAN_CI  Mean and Student-t 95% CI of the finite values in v.
v = v(isfinite(v)); n = numel(v);
r.n = n; r.mean = mean(v);
h = tinv(.975, n - 1) * std(v) / sqrt(n);
r.ci_low = r.mean - h; r.ci_high = r.mean + h; r.positive = sum(v > 0);
end

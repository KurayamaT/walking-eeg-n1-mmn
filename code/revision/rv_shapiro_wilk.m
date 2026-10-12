function [W, p] = rv_shapiro_wilk(x)
% RV_SHAPIRO_WILK  Shapiro-Wilk test, Royston (1992) approximation, 12 <= n <= 5000.
x = sort(x(:)); n = numel(x);
assert(n > 11 && n <= 5000);
m = norminv(((1:n)' - .375) / (n + .25));
ssm = m' * m; c = m / sqrt(ssm); u = 1 / sqrt(n);
an  = c(n)     + .221157*u - .147981*u^2 - 2.071190*u^3 + 4.434685*u^4 - 2.706056*u^5;
an1 = c(n - 1) + .042981*u - .293762*u^2 - 1.752461*u^3 + 5.682633*u^4 - 3.582633*u^5;
phi = (ssm - 2*m(n)^2 - 2*m(n - 1)^2) / (1 - 2*an^2 - 2*an1^2);
a = m / sqrt(phi); a(n) = an; a(n - 1) = an1; a(1) = -an; a(2) = -an1;
W = (a' * x)^2 / sum((x - mean(x)).^2);
ln = log(n);
mu = -1.5861 - .31082*ln - .083751*ln^2 + .0038915*ln^3;
sigma = exp(-.4803 - .082676*ln + .0030302*ln^2);
p = normcdf((log(1 - W) - mu) / sigma, 'upper');
end

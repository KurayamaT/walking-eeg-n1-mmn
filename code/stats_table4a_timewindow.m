function stats_table4a_timewindow(matpath)
% STATS_TABLE4A_TIMEWINDOW
%   Table 4 Panel A: time-window sensitivity of the Sit vs Walk-Free
%   Stimulus x Condition interaction over the frontocentral ROI {Fz, FC1, FC2, Cz}.
%   For each window the interaction contrast per participant is
%       [(dev_walk - dev_sit) - (std_walk - std_sit)]
%   and is tested with a one-sample t-test (t^2 = F(1,n-1)); d_z = mean/SD.
%   Reproduced from the deposited canonical ERP arrays.
%
%   Usage: stats_table4a_timewindow('path/to/erp_arrays_n22.mat')

if nargin < 1 || isempty(matpath)
    matpath = fullfile('..', 'data', 'erp_arrays_n22.mat');
end
S = load(matpath);
dev = S.dev_array; sd = S.std_array; t = S.times_ms(:);
chs = cellstr(S.ch_names); conds = cellstr(S.conditions);
roi  = {'Fz', 'FC1', 'FC2', 'Cz'};
ridx = cellfun(@(c) find(strcmp(chs, c), 1), roi);
csit = find(strcmp(conds, 'sit'), 1);
cw   = find(strcmp(conds, 'walk_free'), 1);
amp  = @(A, w) squeeze(mean(mean(A(w, ridx, :, :), 1), 2)).';   % (subj, cond)

wins = { 'N1 (80-130 ms) ', 80, 130,  '[target t=+2.50 p=.021 dz=+0.53]';
         '100-200 ms     ', 100, 200, '[target t=+1.86 p=.076 dz=+0.40]';
         'MMN (130-200ms)', 130, 200, '[target t=+1.13 p=.270 dz=+0.24]';
         '150-250 ms     ', 150, 250, '[target t=+0.28 p=.781 dz=+0.06]';
         '200-300 ms     ', 200, 300, '[target t=-0.01 p=.992 dz=-0.00]' };

fprintf('\n=== Table 4A: time-window sensitivity of Sit vs Walk-Free interaction (ROI) ===\n');
for k = 1:size(wins, 1)
    w = t >= wins{k, 2} & t <= wins{k, 3};
    dN = amp(dev, w); sN = amp(sd, w);
    ic = (dN(:, cw) - dN(:, csit)) - (sN(:, cw) - sN(:, csit));
    n = numel(ic);
    tt = mean(ic) / (std(ic) / sqrt(n));
    pp = 2 * (1 - tcdf(abs(tt), n - 1));
    dz = mean(ic) / std(ic);
    fprintf('  %s  t(%d)=%+.2f  p=%.3f  d_z=%+.2f   %s\n', wins{k, 1}, n - 1, tt, pp, dz, wins{k, 4});
end
end

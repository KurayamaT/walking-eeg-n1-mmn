function stats_analysis1_n1_interaction(matpath)
% STATS_ANALYSIS1_N1_INTERACTION
%   Analysis 1 headline: the N1-window (80-130 ms) Stimulus x Condition
%   interaction, Sit vs Walk-Free, over the frontocentral ROI {Fz, FC1, FC2, Cz}.
%   Reads the deposited subject-/condition-averaged ERP arrays and reproduces the
%   reported statistic deterministically on any machine.
%
%   Target (machine of record): F(1,21) = 6.26, p = .021, d_z = 0.53,
%   deviant-N1 Sit = -0.764 uV, Walk-Free = -0.453 uV.
%
%   Usage:
%     stats_analysis1_n1_interaction                       % default deposited path
%     stats_analysis1_n1_interaction('path/to/erp_arrays_n22.mat')

if nargin < 1 || isempty(matpath)
    matpath = fullfile('..', 'data', 'erp_arrays_n22.mat');
end
S = load(matpath);

dev   = S.dev_array;      % (time, ch, cond, subj)
sd    = S.std_array;      % (time, ch, cond, subj)
t     = S.times_ms(:);
chs   = cellstr(S.ch_names);
conds = cellstr(S.conditions);

roi   = {'Fz', 'FC1', 'FC2', 'Cz'};
ridx  = cellfun(@(c) find(strcmp(chs, c), 1), roi);
win   = t >= 80 & t <= 130;                 % N1 window
csit  = find(strcmp(conds, 'sit'), 1);
cfree = find(strcmp(conds, 'walk_free'), 1);

% Mean amplitude over the N1 window (dim 1) and the ROI channels (dim 2) -> (subj, cond)
devN1 = squeeze(mean(mean(dev(win, ridx, :, :), 1), 2)).';
stdN1 = squeeze(mean(mean(sd(win, ridx, :, :), 1), 2)).';

% Stimulus x Condition interaction contrast (paired over subjects)
inter = (devN1(:, cfree) - devN1(:, csit)) - (stdN1(:, cfree) - stdN1(:, csit));
n     = numel(inter);
tval  = mean(inter) / (std(inter) / sqrt(n));   % std = sample SD (N-1)
F     = tval ^ 2;
dz    = mean(inter) / std(inter);
p     = 1 - fcdf(F, 1, n - 1);

fprintf('\n=== Analysis 1: N1-window Stimulus x Condition interaction (Sit vs Walk-Free) ===\n');
fprintf('  ROI {Fz, FC1, FC2, Cz}, 80-130 ms, n = %d\n', n);
fprintf('  F(1,%d) = %.2f    p = %.3f    d_z = %.2f\n', n - 1, F, p, dz);
fprintf('  deviant-N1: Sit = %+.3f uV    Walk-Free = %+.3f uV\n', ...
        mean(devN1(:, csit)), mean(devN1(:, cfree)));
fprintf('  [target: F=6.26, p=.021, d_z=0.53, Sit=-0.764, WF=-0.453]\n');
end

function R = rv08_stimulus_sequences(eventdir)
% RV08_STIMULUS_SEQUENCES  Audit of the recorded standard/deviant sequences
%   (Supplementary Methods S-M9; Methods 2.3). eventdir holds one event list
%   per recording (Sxx_<condition>_events.csv; columns sample_1based, is_deviant).
cond = ["sit", "walk_clay", "walk_free", "walk_water"];
T = table(); seqs = cell(22, 4);
for s = 1:22
    for c = 1:4
        E = readtable(fullfile(eventdir, sprintf('S%02d_%s_events.csv', s, cond(c))));
        q = E.is_deviant(:); smp = E.sample_1based(:); a = q(1:end-1); b = q(2:end);
        seqs{s, c} = q;
        T = [T; table(sprintf("S%02d", s), cond(c), numel(q), sum(q == 0), sum(q == 1), sum(a == 0 & b == 0), ...
            sum(a == 1 & b == 0), sum(a == 0 & b == 1), sum(a == 1 & b == 1), min(diff(smp)) * 1000 / 1024, ...
            max(diff(smp)) * 1000 / 1024, 'VariableNames', {'subject', 'condition', 'events', 'standards', ...
            'deviants', 'S_to_S', 'D_to_S', 'S_to_D', 'D_to_D', 'interval_min_ms', 'interval_max_ms'})]; %#ok<AGROW>
    end
end
identical = 0; npairs = 0;
for s = 1:22
    for a = 1:3
        for b = a+1:4
            n = min(numel(seqs{s, a}), numel(seqs{s, b})); npairs = npairs + 1;
            identical = identical + all(seqs{s, a}(1:n) == seqs{s, b}(1:n));
        end
    end
end
pairs = T.S_to_S + T.D_to_S + T.S_to_D + T.D_to_D;
assert(all(pairs == T.events - 1));
R.records = T;
R.totals = table(height(T), sum(pairs), sum(T.D_to_D), 100 * sum(T.D_to_D) / sum(pairs), sum(T.D_to_D > 0), ...
    min(T.D_to_D), max(T.D_to_D), npairs, identical, sum(T.standards), sum(T.deviants), ...
    sum(T.standards) / sum(T.deviants), 'VariableNames', {'records', 'adjacent_pairs', 'deviant_to_deviant_pairs', ...
    'deviant_to_deviant_pct', 'records_with_consecutive_deviants', 'DD_min', 'DD_max', 'condition_pairs', ...
    'identical_condition_pairs', 'standards', 'deviants', 'standard_to_deviant_ratio'});
end

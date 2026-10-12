function R = rv07_speed_order(csvpath)
% RV07_SPEED_ORDER  Self-selected walking speed and allocation of the four
%   conditions to block positions (Methods 2.2; opening of Results).
T = readtable(csvpath, 'TextType', 'string');
assert(height(T) == 22);
grp = {"all", true(22, 1); "male", T.sex == "male"; "female", T.sex == "female"};
S = table();
for k = 1:3
    v = T.speed_m_min(grp{k, 2});
    S = [S; table(grp{k, 1}, numel(v), mean(v), std(v), min(v), max(v), 'VariableNames', ...
        {'group', 'n', 'mean_m_min', 'sd_m_min', 'min_m_min', 'max_m_min'})]; %#ok<AGROW>
end
O = [T.Sit_position, T.Walk_Free_position, T.Walk_Water_position, T.Walk_Clay_position];
assert(all(all(sort(O, 2) == 1:4)));
N = zeros(4);
for j = 1:4, for k = 1:4, N(j, k) = sum(O(:, j) == k); end, end
R.speed = S;
R.positions = array2table(N, 'VariableNames', {'position1', 'position2', 'position3', 'position4'}, ...
    'RowNames', {'Sit', 'Walk-Free', 'Walk-Water', 'Walk-Clay'});
R.unique_orders = size(unique(O, 'rows'), 1);
end

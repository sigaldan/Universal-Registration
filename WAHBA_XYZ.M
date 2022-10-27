function [data_b_new,C_b2a] = Wahba_XYZ(data_a, data_b, sns_pos_a, sns_pos_b)
% find rotation from "b" to "a"
% and correct "b" measurements (in "b" axis)

data_a(:,1) = data_a(:,1) + sns_pos_a(1) - sns_pos_b(1);
data_a(:,2) = data_a(:,2) + sns_pos_a(2) - sns_pos_b(2);
data_a(:,3) = data_a(:,3) + sns_pos_a(3) - sns_pos_b(3);

%% b-->a (BAT_b)

[C_b2a, ~] = solve_wahba_svd(data_a, data_b, 0);

data_b_new = C_b2a*data_b';
data_b_new = data_b_new';

% 

end
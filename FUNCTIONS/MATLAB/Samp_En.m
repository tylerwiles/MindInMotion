function ent = Samp_En(dat, pattern_length, tolerance, standard_dev)

% INPUTS:
% dat = Vector time series.
% pattern_length = Length of a vector to be compared.
% tolerance = Criteria of similarity.
% standard_dev = Standard Deviation of dat

% OUTPUT:
% ent = Sample entropy.

%%
% TITLE: Samp_En.m
% DATE: December 1, 2023
% AUTHOR: Aaron D. Likens, PHD
% EMAIL: alikens@unomaha.edu

% DESCRIPTION:
% Function to calculate sample entropy.

% Reference: Richman, J. S., & Moorman, J. R. (2000). Physiological
% time-series analysis using approximate entropy and sample entropy.
% American Journal of Physiology-Heart and Circulatory Physiology,
% 278(6), H2039-H2049.

% Copyright 2023, Aaron D. Likens

% Redistribution and use of this script, with or without
% modification, is permitted provided this copyright notice,
% the original authors name and the following disclaimer remains.

% DISCLAIMER: It is the user's responsibility to check the code is returning
% the appropriate results before reporting any scientific findings.
% The author will not accept any responsibility for incorrect results
% directly relating to user error.

tolerance = tolerance*standard_dev;

x_mat = buffer(dat, pattern_length + 1, pattern_length, 'nodelay')';
match_m = sum(pdist(x_mat, 'chebychev') < tolerance);
match_m1 = sum(pdist(x_mat(:, 1:pattern_length), 'chebychev') < tolerance);

if match_m == 0 || match_m1 == 0
    ent = NaN;
    return
end

ent = log(match_m1 / match_m);

end

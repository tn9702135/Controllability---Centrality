function ICC = icc_3_1(X)
% ICC(3,1): two-way mixed, consistency, single measure
% Rows = regions, Columns = subjects

[n, k] = size(X);   % n regions, k subjects

% Means
mean_region = mean(X, 2);
mean_subject = mean(X, 1);
grand_mean = mean(X(:));

% Sum of squares
SSR = k * sum((mean_region - grand_mean).^2);   % regions
SSC = n * sum((mean_subject - grand_mean).^2);  % subjects
SSE = sum(sum((X - mean_region - mean_subject + grand_mean).^2));

% Degrees of freedom
dfR = n - 1;
dfE = (n - 1) * (k - 1);

% Mean squares
MSR = SSR / dfR;
MSE = SSE / dfE;

% ICC(3,1)
ICC = (MSR - MSE) / (MSR + (k - 1) * MSE);
end

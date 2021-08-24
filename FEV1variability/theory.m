%% confidence interval coefficients for the sample's standard deviation
% CI = [coef_down * SD, coef_up * SD]
% there is 95% chance that this interval contains the true sample's std
% CIs look at the reliability of the estimation procedure, not 

% assumption: normal distribution

alpha = 0.05;
n = [2,3,5,10,25,50,70,100,200,500,1000]';
std_coef_down = sqrt((n-1)./chi2inv(1-alpha/2,n-1));
std_coef_up = sqrt((n-1)./chi2inv(alpha/2,n-1));
diff_std_coef_down = [diff(std_coef_down).*100; nan];
diff_std_ci_up = [diff(std_coef_up).*100; nan];
disp(table(n,std_coef_down,diff_std_coef_down,std_coef_up,diff_std_ci_up));

% after 25, the diff is quite low

%% for sample's mean

% samples come from a normal distribution

% clarify for sample mean when we don't know the sample's std
% known std dev: CI = mean +/- std(x)/sqrt(n)*1.96 
% unknown std dev use the Student's distribution

% how is it possible that the CI is based on the sample's population and
% independant from the global population?

% why in stats 30 is significative to make statistics? works only for mean

% chebyshev inequality assumptions
% CI assumptions
% t test, z test
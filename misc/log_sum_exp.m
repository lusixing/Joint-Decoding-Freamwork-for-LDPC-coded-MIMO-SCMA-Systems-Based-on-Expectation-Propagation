function y = log_sum_exp(x)
xm = max(x);
x = x-repmat(xm, size(x,1), 1);
y = xm + log(sum(exp(x)));
function bv = bv_fun_test(x)
bv = sin(pi*x(:, 1)) .* exp(pi*x(:, 2));
end
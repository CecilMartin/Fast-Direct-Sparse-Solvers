function bv = bv_fun_test(x)
bv = sin(pi*x(:, 1)) .* cos(pi*x(:, 2));
end
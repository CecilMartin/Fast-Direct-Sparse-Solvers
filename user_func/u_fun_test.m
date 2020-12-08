function u = u_fun_test(x)
u = sin(pi*x(:, 1)) .* exp(pi*x(:, 2));
end
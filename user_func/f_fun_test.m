function f = f_fun_test(x)
f = 2*pi^2 * sin(pi*x(:, 1)) .* cos(pi*x(:, 2));
end
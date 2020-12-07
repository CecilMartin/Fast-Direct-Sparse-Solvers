function f = get_f_orig(n,bv_fun,f_fun)
% load vector for Laplace
% Dirichlet Boundary Value
% The region is set as [0, 1] * [0, 1]
% With (n-1)*(n-1) interior points
% f = -ones((n-1)^2,1);
h = 1/n;
[x1, x2] = ndgrid((1:(n-1))*h); x = [x2(:) 1-x1(:)];
f = h^2 * f_fun(x);
Il = 1:n; f(Il) = f(Il) + bv_fun(x(Il, :));  % The left boundary, {0}*[0, 1]
Ir = (n-1)^2 - (1:n); f(Ir) = f(Ir) + bv_fun(x(Ir, :));  % The right boundary, {1}*[0, 1]
Iu = 1:(n-1):(n-1)^2; f(Iu) = f(Iu) + bv_fun(x(Iu, :));  % The upper boundary, [0, 1]*{1}
Id = (n-1):(n-1):(n-1)^2; f(Id) = f(Id) + bv_fun(x(Id, :));  % The down boundary, [0, 1]*{0}
end